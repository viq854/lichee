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

/**
 * Main cell lineage construction pipeline
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
		SNVDatabase db = new SNVDatabase(args.inputFileName, args.normalSampleId);
		
		// 2. handle normal cell contamination, CNVs, 
		//    determine the minimum number of clusters using LOH
		//    (+ any additional filtering, pre-processing)
		
		// 3. get the SNVs partitioned by group tag and create the appropriate SNV group objects
		HashMap<String, ArrayList<SNVEntry>> snvsByTag = db.generateFilteredTAG2SNVsMap(null);
		ArrayList<SNVGroup> groups = new ArrayList<SNVGroup>();
		for(String groupTag : snvsByTag.keySet()) {
			groups.add(new SNVGroup(groupTag, snvsByTag.get(groupTag), db.getNumRobustSNVs(groupTag), db.isRobust(groupTag)));
		}
		
		// 4. cluster SNVs in each group
		AAFClusterer clusterer = new AAFClusterer();
		for(SNVGroup group : groups) {
			Cluster[] clusters = clusterer.clusterSubPopulations(group, ClusteringAlgorithms.EM, 1);
			logger.fine("Clustering results for group: " + group.getTag());
			for(Cluster c : clusters) {
				logger.fine(c.toString());
				logger.log(Level.FINE, c.toString());
			}
			group.setSubPopulations(clusters);
		}
		
		// 5. construct the constraint network
		PHYNetwork constrNetwork = new PHYNetwork(groups, db.getNumofSamples());
		logger.log(Level.FINE, constrNetwork.toString());
		logger.log(Level.FINE, constrNetwork.getNodesAsString());
		
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
		logger.log(Level.FINE, "Found " + spanningTrees.size() + " valid trees");
		
		if(spanningTrees.size() == 0) {
			logger.log(Level.INFO, "No valid trees found. Adjusting the constraints...");	
			// if no valid trees were found, fix the network (e.g. remove group nodes that are not robust)
			constrNetwork = constrNetwork.fixNetwork();
			spanningTrees = constrNetwork.getLineageTrees();  
			logger.log(Level.INFO, "Found " + spanningTrees.size() + " valid trees after constraint adjustment");	
		}
		
		// 7. evaluate/rank the trees
		constrNetwork.evaluateLineageTrees();
		if(spanningTrees.size() > 0) {
			logger.log(Level.FINE, "Best tree error score: " + spanningTrees.get(0).getErrorScore());
		}
		
		// 8. result visualization
		String[] sampleNames = new String[db.getNumofSamples()];
		for(int i = 0; i < db.getNumofSamples(); i++) {
			sampleNames[i] = db.getName(i);
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
			//for(int i = 0; i < db.getNumofSamples(); i++) {
				//System.out.println(spanningTrees.get(0).getLineage(i, sampleNames[i]));
			//}
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
		
		SNVDatabase db = new SNVDatabase(args.inputFileName, args.normalSampleId);
		HashMap<String, ArrayList<SNVEntry>> snvsByTag = db.generateFilteredTAG2SNVsMap(null);
		
		PHYNetwork net = readNetworkFromFile(netFileName);
		String[] sampleNames = new String[net.numSamples];
		ArrayList<PHYTree> trees = readTreesFromFile(treeFileName, sampleNames);
		if(args.numShow == 1) {
			net.displayTree(trees.get(0), sampleNames, snvsByTag, null);
		} else if(args.numShow > 1) {
			for(int i = 0; i < args.numShow; i++) {
				if(trees.size() < i) {
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
			e.printStackTrace();
			System.err.println("Failed to read to the file: " + netFileName);
			System.exit(1);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			System.err.println("Failed to load network from file: " + netFileName);
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
		options.addOption("o", true, "Output file path (default: input file path");
		options.addOption("vcf", false, "Input data type is vcf (default: validation)");
		options.addOption("nid", "normal", true, "Normal sample id (default: 0)");
		options.addOption("s", "save", false, "Save the output to file");
		options.addOption("d", "showBestTree", false, "Display the best lineage tree");
		options.addOption("n", "showNetwork", false, "Display the constraint network");
		options.addOption("v", "verbose", false, "Verbose mode");
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmdLine = null;
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		// Set-up input args
		Args params = new Args();	
		
		// required options
		if(!cmdLine.hasOption("i")) {
			new HelpFormatter().printHelp("lineage", options);
			System.exit(-1);
		}
		
		if(cmdLine.hasOption("nid")) {
			params.normalSampleId = Integer.parseInt(cmdLine.getOptionValue("nid"));
		}
		if(cmdLine.hasOption("s")) {
			params.persist = true;
		}
		if(cmdLine.hasOption("d")) {
			params.showBest = true;
		}
		if(cmdLine.hasOption("n")) {
			params.showNetwork = true;
		}
		if(cmdLine.hasOption("vcf")) {
			params.isVCF = true;
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
		
		if(cmdLine.hasOption("build")) {
			// input file
			params.inputFileName = cmdLine.getOptionValue("i");
			// output file	
			if(cmdLine.hasOption("o")) {
				params.outputFilePrefix = cmdLine.getOptionValue("o");	
			} else {
				params.outputFilePrefix = params.inputFileName;	
			}
			buildLineage(params);
		} else if (cmdLine.hasOption("show")){
			params.showFileNamePrefix = cmdLine.getOptionValue("i");
			showTrees(params);
		} else {
			new HelpFormatter().printHelp("lineage", options);
			System.exit(-1);
		}
		//buildLineage("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_2/Patient_2.validation.txt", 0);
	}
	
	protected static class Args {
		// --- 'build' command ---
		String inputFileName;
		String outputFilePrefix;	
		int normalSampleId = 0;
		
		// --- 'show' command ---
		String showFileNamePrefix;
		int treeId = 0;
		int numShow = 1;
		
		// flags
		boolean isVCF = false;
		boolean showBest = true;
		boolean showNetwork = false;
		boolean verbose = false;
		boolean persist = true;
	}

	protected static class LogFormatter extends Formatter {
		public String format(LogRecord rec) {
			return rec.getMessage() + "\r\n";
		}
	}
}
