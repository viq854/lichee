package lineage;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.ConsoleHandler;
import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import lineage.AAFClusterer.Cluster;
import lineage.AAFClusterer.ClusteringAlgorithms;
import lineage.PHYTree;
import util.*;
import util.Configs.format;

/**
 * Main cell lineage construction pipeline
 * 
 * @autor viq
 */
public class LineageEngine {

	private static final String NET_FILE_EXTENSION = ".net";
	private static final String TREE_FILE_EXTENSION = ".tree";
	private static final String TREE_JPG_FILE_EXTENSION = ".tree.jpg";
	private static final String SAMPLE_LIN_FILE_EXTENSION = ".lin";
	protected static final Logger logger = Logger.getLogger("lineage.engine");
	
	/**
	 * The main pipeline for constructing the cell lineage
	 */
	public static void buildLineage(Args args) {
				
		// 1. load validation/VCF data
		SNVDataStore db = new SNVDataStore(args.inputFileName, args.normalSampleId);
		db.annotateSNVs(args.cnvFileName, args.annFileName, args.cosmicFileName, args.tcgaFileName);
		
		// [2. normal cell contamination, CNVs (+ any additional filtering, pre-processing)]
		
		// 3. get the SNVs partitioned by group tag and create the appropriate SNV group objects
		HashMap<String, ArrayList<SNVEntry>> snvsByTag = db.getSomaticGroups();
		ArrayList<SNVGroup> groups = new ArrayList<SNVGroup>();
		for(String groupTag : snvsByTag.keySet()) {
			groups.add(new SNVGroup(groupTag, snvsByTag.get(groupTag), db.isRobustGroup(groupTag)));
		}
		System.out.println("Total number of somatic SNV groups: " + groups.size());
		if(groups.size() == 0) {
			logger.log(Level.WARNING, "All SNV groups have been filtered out.");
			return;
		}
		
		// 4. cluster SNVs in each group
		AAFClusterer clusterer = new AAFClusterer();
		for(SNVGroup group : groups) {
			Cluster[] clusters = clusterer.clusterSubPopulations(group, ClusteringAlgorithms.EM, 1);
			logger.fine("Clustering results for group: " + group.getTag());
			for(Cluster c : clusters) {
				logger.log(Level.FINE, c.toString());
			}
			group.setSubPopulations(clusters);
		}
		
		// 5. construct the constraint network
		PHYNetwork constrNetwork = new PHYNetwork(groups, db.getNumSamples());
		logger.log(Level.FINE, constrNetwork.toString());
		logger.log(Level.INFO, "Nodes:\n" + constrNetwork.getNodesAsString());
		
		// -- temporary -- evaluating results
		String nodesFileName = args.outputFilePrefix + ".nodes";
		try {
			FileWriter fw = new FileWriter(nodesFileName);
			fw.write(constrNetwork.getNodesAsString());
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + nodesFileName);
			System.exit(-1);
		}
		
		// 6. find all the lineage trees that pass the AAF constraints
		ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees();  
		logger.log(Level.INFO, "Found " + spanningTrees.size() + " valid trees");
		
		if(spanningTrees.size() == 0) {
			logger.log(Level.INFO, "No valid trees found. Adjusting the constraints...");	
			// if no valid trees were found, fix the network (e.g. remove group nodes that are not robust)
			int delta = 0;
			do {
				int numNodes = constrNetwork.numNodes;
				constrNetwork = constrNetwork.fixNetwork();
				spanningTrees = constrNetwork.getLineageTrees();  
				logger.log(Level.INFO, "Found " + spanningTrees.size() + " valid trees after constraint adjustment");	
				delta = numNodes - constrNetwork.numNodes; 
			} while((delta != 0) && (spanningTrees.size() <= 0));
		}
		
		// 7. evaluate/rank the trees
		constrNetwork.evaluateLineageTrees();
		if(spanningTrees.size() > 0) {
			logger.log(Level.INFO, "Best tree error score: " + spanningTrees.get(0).getErrorScore());
		}
		
		// 8. result visualization
		String[] sampleNames = new String[db.getNumSamples()];
		logger.log(Level.INFO, "Samples: ");
		for(int i = 0; i < db.getNumSamples(); i++) {
			sampleNames[i] = db.getSampleName(i);
			logger.log(Level.INFO, i + ": " + sampleNames[i]);
		}
		
		if(args.showNetwork) {
			constrNetwork.displayNetwork();
		}
		if(args.showBest && (spanningTrees.size() > 0)) {
			if(args.persist) {
				constrNetwork.displayTree(spanningTrees.get(0), sampleNames, null, args.outputFilePrefix + TREE_JPG_FILE_EXTENSION);
			} else {
				constrNetwork.displayTree(spanningTrees.get(0), sampleNames, null, null);
			}
			for(int i = 1; i < args.numShow; i++) {
				if(spanningTrees.size() > i) {
					constrNetwork.displayTree(spanningTrees.get(i), sampleNames, null, null);
				} else {
					break;
				}
			}
		}
		
		// 9. persistent storage
		if(args.persist) {
			writeNetworkToFile(constrNetwork, args.outputFilePrefix);
			if(spanningTrees.size() > 0) {
				writeTreesToFile(spanningTrees, sampleNames, args.outputFilePrefix);
				writeSampleLineageToFile(spanningTrees, sampleNames,args.outputFilePrefix);
			}
		}
		
	}
	
	/** 
	 * The main pipeline for tree result visualization
	 */
	public static void showTrees(Args args) {
		String netFileName = args.showFileNamePrefix + NET_FILE_EXTENSION;
		String treeFileName = args.showFileNamePrefix + TREE_FILE_EXTENSION;
		
		SNVDatabase db = new SNVDatabase(args.showFileNamePrefix, args.normalSampleId);
		db.resolveNonRobustconflicts();
		HashMap<String, ArrayList<SNVEntry>> snvsByTag = db.generateFilteredTAG2SNVsMap(null);
		
		PHYNetwork net = readNetworkFromFile(netFileName);
		String[] sampleNames = new String[net.numSamples];
		ArrayList<PHYTree> trees = readTreesFromFile(treeFileName, sampleNames);
		if(args.numShow == 1) {
			net.displayTree(trees.get(0), sampleNames, snvsByTag, null);
		} else if(args.numShow > 1) {
			for(int i = 0; i < args.numShow; i++) {
				if(trees.size() > i) {
					net.displayTree(trees.get(i), sampleNames, snvsByTag, null);
				} else {
					break;
				}
			}
		}
	}
	
	// I/O
	
	private static void writeNetworkToFile(PHYNetwork net, String outputNamePrefix) {
		String netFileName = outputNamePrefix + NET_FILE_EXTENSION;
		try {
			FileOutputStream fos = new FileOutputStream(netFileName);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(net);
			oos.close();
			fos.close(); // move to finally
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + netFileName);
			System.exit(-1);
		}
	}
	
	private static void writeTreesToFile(ArrayList<PHYTree> trees, String[] sampleNames, String outputNamePrefix) {
		String treeFileName = outputNamePrefix + TREE_FILE_EXTENSION;
		try {
			FileOutputStream fos = new FileOutputStream(treeFileName);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(sampleNames);
			oos.writeObject(trees);
			oos.close();
			fos.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + treeFileName);
			System.exit(-1);
		}
	}
	
	private static void writeSampleLineageToFile(ArrayList<PHYTree> trees, String[] sampleNames, String outputNamePrefix) {
		String linFileName = outputNamePrefix + SAMPLE_LIN_FILE_EXTENSION;
		try {
			FileWriter fw = new FileWriter(linFileName);
			for(int i = 0; i < sampleNames.length; i++) {
				fw.write(trees.get(0).getLineage(i, sampleNames[i]));
				fw.write("\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + linFileName);
			System.exit(-1);
		}
	}
	
	private static PHYNetwork readNetworkFromFile(String netFileName) {
		try {
			FileInputStream fis = new FileInputStream(netFileName);
			ObjectInputStream ois = new ObjectInputStream(fis);
			PHYNetwork net = (PHYNetwork) ois.readObject();
			ois.close();
			fis.close();
			return net;
		} catch (IOException e) {
			//e.printStackTrace();
			System.err.println("Failed to read to the file: " + netFileName + "\nNote: The -build command (with -s option) has to be run before the -show command for the given file.");
			System.exit(1);
		} catch (ClassNotFoundException e) {
			//e.printStackTrace();
			System.err.println("Failed to load network from file: " + netFileName + "\nNote: The -build command (with -s option) has to be run before the -show command for the given file.");
			System.exit(1);
		}
		return null;
	}
	
	private static ArrayList<PHYTree> readTreesFromFile(String treeFileName, String[] outSampleNames) {
		try {
			FileInputStream fis = new FileInputStream(treeFileName);
			ObjectInputStream ois = new ObjectInputStream(fis);
			String[] sampleNames = (String[]) ois.readObject();
			@SuppressWarnings("unchecked")
			ArrayList<PHYTree> trees = (ArrayList<PHYTree>) ois.readObject();
			ois.close();
			fis.close();
			for(int i = 0; i < sampleNames.length; i++) {
				outSampleNames[i] = sampleNames[i];
			}
			
			return trees;
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to read to the file: " + treeFileName);
			System.exit(1);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			System.err.println("Failed to read to the file: " + treeFileName);
			System.exit(1);
		}
		return null;
	}

	// ---- LAUNCH ----
	
	public static void main(String[] args) {
		Options options = new Options(); 
		options.addOption("build", false, "Construct the sample cell lineage trees");
		options.addOption("show", false, "Display the saved cell lineage tree(s)");
		options.addOption("i", true, "Input file path");
		options.addOption("o", true, "Output file path (default: input file path)");
		options.addOption("cnv", true, "File path to CNV regions used to mark SNVs");
		options.addOption("ann", true, "File path to ANNOVAR SNV annotations");
		options.addOption("cosmic", true, "File path to COSMIC SNV annotations");
		options.addOption("tcga", true, "File path to TCGA SNV annotations");
		options.addOption("n", "normal", true, "Normal sample id (default: 0)");
		options.addOption("closestParentOnly", false, "Only add edges to the ancestor in the closest possible level");
		options.addOption("maxAAF", true, "Maximum allowed AAF in a sample (default: 0.6)");
		options.addOption("minAAFPresent", true, "Minimum AAF to robustly call an SNV in a sample (default: 0.04)");
		options.addOption("maxAAFAbsent", true, "Maximum AAF to robustly consider an SNV as absent from a sample (default: 0.015)");
		options.addOption("minGroupSize", true, "Minimum SNV group size (default: 1)");
		options.addOption("minRobustGroup", true, "Minimum number of robust SNVs in group to be robust (default: 2)");
		options.addOption("minRobustGroupKeep", true, "Minimum number of robust SNVs per group to keep the group in the network initially (default: 0 - keeps all the groups)");
		options.addOption("minTargetDistRatio", true, "Minimum non-robust SNV to target group's SNV ratio per sample (default: 0.5)");
		options.addOption("e", true, "AAF error margin (default: 0.08)");
		options.addOption("minClusterSize", true, "Minimum size a cluster must have to be a considered a node in the network (default: 2)");
		options.addOption("maxClusterDist", true, "Maximum mean AAF difference up to which two clusters can be collapsed (default: 0.2)");
		options.addOption("vcf", false, "Input data type is vcf (default: validation)");
		options.addOption("s", "save", false, "Save the output to file");
		options.addOption("tree", "showTree", false, "Display the top ranking lineage tree");
		options.addOption("net", "showNetwork", false, "Display the constraint network");
		options.addOption("top", true, "Number of top ranking lineage trees to display");
		options.addOption("v", "verbose", false, "Verbose mode");
		options.addOption("h", "help", false, "Print usage");
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmdLine = null;
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			new HelpFormatter().printHelp("smuth", options);
			System.exit(-1);
		}
		
		// Set-up input args
		Args params = new Args();	
		
		// required options
		if(!cmdLine.hasOption("i")) {
			System.out.println("Required parameter: input file path [-i]");
			new HelpFormatter().printHelp("smuth", options);
			System.exit(-1);
		}
		if(!cmdLine.hasOption("n")) {
			System.out.println("Required parameter: normal sample id [-n]");
			new HelpFormatter().printHelp("smuth", options);
			System.exit(-1);
		}
		
		if(cmdLine.hasOption("n")) {
			params.normalSampleId = Integer.parseInt(cmdLine.getOptionValue("n"));
		}
		if(cmdLine.hasOption("minAAFPresent")) {
			Configs.VALIDATION_THR = Double.parseDouble(cmdLine.getOptionValue("minAAFPresent"));
		}
		if(cmdLine.hasOption("maxAAFAbsent")) {
			Configs.VALIDATION_SOFT_THR = Double.parseDouble(cmdLine.getOptionValue("maxAAFAbsent"));
		}
		if(cmdLine.hasOption("minGroupSize")) {
			Configs.GROUP_SIZE_THR = Integer.parseInt(cmdLine.getOptionValue("minGroupSize"));
		}
		if(cmdLine.hasOption("maxAAF")) {
			Configs.MAX_ALLOWED_VAF = Integer.parseInt(cmdLine.getOptionValue("maxAAF"));
		}
		if(cmdLine.hasOption("minRobustGroup")) {
			Configs.ROBUSTGROUP_SIZE_THR = Integer.parseInt(cmdLine.getOptionValue("minRobustGroup"));
		}
		if(cmdLine.hasOption("minRobustGroupKeep")) {
			Configs.GROUP_ROBUST_NUM_THR = Integer.parseInt(cmdLine.getOptionValue("minRobustGroupKeep"));
		}
		if(cmdLine.hasOption("minTargetDistRatio")) {
			Configs.MIN_VAF_TARGET_RATIO_PER_SAMPLE = Double.parseDouble(cmdLine.getOptionValue("minTargetDistRatio"));
		}
		if(cmdLine.hasOption("e")) {
			Parameters.AAF_ERROR_MARGIN = Double.parseDouble(cmdLine.getOptionValue("e"));
		}
		if(cmdLine.hasOption("minClusterSize")) {
			Parameters.MIN_CLUSTER_SIZE = Integer.parseInt(cmdLine.getOptionValue("minClusterSize"));
		}
		if(cmdLine.hasOption("maxClusterDist")) {
			Parameters.MAX_COLLAPSE_CLUSTER_DIFF = Double.parseDouble(cmdLine.getOptionValue("maxClusterDist"));
		}
		if(cmdLine.hasOption("vcf")) {
			Configs.INFORMAT = format.VCF;
		}
		if(cmdLine.hasOption("s")) {
			params.persist = true;
		}
		if(cmdLine.hasOption("tree")) {
			params.showBest = true;
		}
		if(cmdLine.hasOption("net")) {
			params.showNetwork = true;
		}
		if(cmdLine.hasOption("closestParentOnly")) {
			Parameters.ALL_EDGES = false;
		}
		if(cmdLine.hasOption("top")) {
			params.numShow = Integer.parseInt(cmdLine.getOptionValue("top"));
		}
		if(cmdLine.hasOption("h")) {
			new HelpFormatter().printHelp(" ", options);
		}
		// logger
		ConsoleHandler h = new ConsoleHandler();
		h.setFormatter(new LogFormatter());
		h.setLevel(Level.INFO);
		logger.setLevel(Level.INFO);
		if(cmdLine.hasOption("v")) {
			h.setLevel(Level.FINEST);
			logger.setLevel(Level.FINEST);
		}
		logger.addHandler(h);
		logger.setUseParentHandlers(false);
		
		if(cmdLine.hasOption("build")) {
			// input file
			params.inputFileName = cmdLine.getOptionValue("i");
			// output file	
			if(cmdLine.hasOption("o")) {
				params.outputFilePrefix = cmdLine.getOptionValue("o");	
			} else {
				params.outputFilePrefix = params.inputFileName;	
			}
			if(cmdLine.hasOption("cnv")) {
				params.cnvFileName = cmdLine.getOptionValue("cnv");	
			}
			if(cmdLine.hasOption("ann")) {
				params.annFileName = cmdLine.getOptionValue("ann");	
			}
			if(cmdLine.hasOption("cosmic")) {
				params.cosmicFileName = cmdLine.getOptionValue("cosmic");	
			}
			if(cmdLine.hasOption("tcga")) {
				params.tcgaFileName = cmdLine.getOptionValue("tcga");	
			}
			buildLineage(params);
		} else if (cmdLine.hasOption("show")){
			params.showFileNamePrefix = cmdLine.getOptionValue("i");
			showTrees(params);
		} else {
			System.out.println("Required command: build OR show");
			new HelpFormatter().printHelp("smuth", options);
			System.exit(-1);
		}
	}
	
	protected static class Args {
		// --- 'build' command ---
		String inputFileName;
		String cnvFileName;
		String annFileName;
		String cosmicFileName;
		String tcgaFileName;
		String outputFilePrefix;	
		int normalSampleId = 0;
		
		// --- 'show' command ---
		String showFileNamePrefix;
		int numShow = 1;
		
		// flags
		boolean showBest = false;
		boolean showNetwork = false;
		boolean verbose = false;
		boolean persist = false;
	}

	protected static class LogFormatter extends Formatter {
		public String format(LogRecord rec) {
			return rec.getMessage() + "\r\n";
		}
	}
}
