package lineage;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Comparator;
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
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import lineage.AAFClusterer.Cluster;
import lineage.AAFClusterer.ClusteringAlgorithms;
import lineage.PHYTree;

/**
 * Main cell lineage construction pipeline
 * 
 * @autor viq
 */
public class LineageEngine {

	private static final String TREES_TXT_FILE_EXTENSION = ".trees.txt";
	private static final String NET_OBJ_FILE_EXTENSION = ".net.obj";
	private static final String TREES_OBJ_FILE_EXTENSION = ".trees.obj";
	private static final String TREE_JPG_FILE_EXTENSION = ".tree.jpg";
	protected static final Logger logger = Logger.getLogger("lineage.engine");
	
	/**
	 * The main pipeline for constructing the cell lineage
	 */
	public static void buildLineage(Args args) {
				
		// 1. load validation/VCF data
		SNVDataStore db = new SNVDataStore(args.inputFileName, args.normalSampleId);
		db.annotateSNVs(args.cnvFileName, args.annFileName, args.cosmicFileName, args.tcgaFileName);
		
		// [+ any additional filtering, pre-processing]
		
		// 3. get the SNVs partitioned by group tag and create the appropriate SNV group objects
		HashMap<String, ArrayList<SNVEntry>> snvsByTag = db.getSomaticGroups();
		ArrayList<SNVGroup> groups = new ArrayList<SNVGroup>();
		for(String groupTag : snvsByTag.keySet()) {
			groups.add(new SNVGroup(groupTag, snvsByTag.get(groupTag), db.isRobustGroup(groupTag)));
		}
		if(groups.size() == 0) {
			logger.log(Level.WARNING, "All SNV groups have been filtered out.");
			return;
		}
		
		// 4. cluster SNVs in each group
		AAFClusterer clusterer = new AAFClusterer();
		for(SNVGroup group : groups) {
			Cluster[] clusters = clusterer.clusterSubPopulations(group, ClusteringAlgorithms.EM, 3);
			logger.fine("Clustering results for group: " + group.getTag());
			for(Cluster c : clusters) {
				logger.log(Level.FINE, c.toString());
			}
			group.setSubPopulations(clusters);
		}
		
		// 5. construct the constraint network
		PHYNetwork constrNetwork = new PHYNetwork(groups, db.getNumSamples());
		String[] sampleNames = new String[db.getNumSamples()];
		logger.log(Level.INFO, "Samples: ");
		for(int i = 0; i < db.getNumSamples(); i++) {
			sampleNames[i] = db.getSampleName(i);
			logger.log(Level.INFO, i + ": " + sampleNames[i]);
		}
		logger.log(Level.FINE, constrNetwork.toString());
		
		// 6. find all the lineage trees that pass the AAF constraints
		ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees();  
		logger.log(Level.INFO, "Found " + spanningTrees.size() + " valid trees");
		
		if(spanningTrees.size() == 0) {
			logger.log(Level.INFO, "No valid trees found. Adjusting the network...");	
			// if no valid trees were found, fix the network (e.g. remove group nodes that are not robust)
			int delta = 0;
			do {
				int numNodes = constrNetwork.numNodes;
				constrNetwork = constrNetwork.fixNetwork();
				spanningTrees = constrNetwork.getLineageTrees();  
				logger.log(Level.INFO, "Found " + spanningTrees.size() + " valid trees after network adjustment");	
				delta = numNodes - constrNetwork.numNodes; 
			} while((delta != 0) && (spanningTrees.size() <= 0));
		}
		
		// 7. evaluate/rank the trees
		constrNetwork.evaluateLineageTrees();
		if(spanningTrees.size() > 0) {
			logger.log(Level.INFO, "Best tree error score: " + spanningTrees.get(0).getErrorScore());	
			logger.log(Level.INFO, spanningTrees.get(0).toString());
		} 
		
		// 8. result visualization
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
			writeNetworkToObjFile(constrNetwork, args);
			writeTreesToObjFile(spanningTrees, sampleNames, args);
			writeTreesToTxtFile(constrNetwork, spanningTrees, sampleNames, args);
		}	
	}
	
	/** 
	 * The main pipeline for tree result visualization
	 */
	public static void showTrees(Args args) {
		String netFileName = args.showFileNamePrefix + NET_OBJ_FILE_EXTENSION;
		String treeFileName = args.showFileNamePrefix + TREES_OBJ_FILE_EXTENSION;
		
		SNVDataStore db = new SNVDataStore(args.showFileNamePrefix, args.normalSampleId);
		HashMap<String, ArrayList<SNVEntry>> snvsByTag = db.getSomaticGroups();

		PHYNetwork net = readNetworkFromObjFile(netFileName);
		String[] sampleNames = new String[net.numSamples];
		ArrayList<PHYTree> trees = readTreesFromObjFile(treeFileName, sampleNames);
		for(int i = 0; i < args.numShow; i++) {
			if(trees.size() > i) {
				net.displayTree(trees.get(i), sampleNames, snvsByTag, null);
			} else {
				break;
			}
		}
		
	}
	
	// I/O
	private static void writeTreesToTxtFile(PHYNetwork net, ArrayList<PHYTree> trees, String[] sampleNames, Args args) {
		String treeFileName = args.outputFilePrefix + TREES_TXT_FILE_EXTENSION;
		try {
			FileWriter fw = new FileWriter(treeFileName);
			String nodes = net.getNodesWithMembersAsString();
			if(trees.size() > 0 && args.numSave > 0) {
				fw.write("Tree Nodes: \n");
				fw.write(nodes);
				fw.write("\n");
			}
			for(int i = 0; i < args.numSave; i++) {
				if(trees.size() > i) {
					fw.write("****Tree " + i + "****\n");
					String edges = trees.get(i).toString();
					fw.write(edges);
					fw.write("Error score: " + trees.get(i).getErrorScore()+"\n\n");	
					fw.write("Sample decomposition: \n");
					String lineage = "";
					for(int j = 0; j < sampleNames.length; j++) {
						lineage += trees.get(i).getLineage(j, sampleNames[j]);
						lineage += "\n";
					}
					fw.write(lineage);
					fw.write("\n");
				}
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + treeFileName);
			System.exit(-1);
		}
	}
	
	private static void writeNetworkToObjFile(PHYNetwork net, Args args) {
		String netFileName = args.outputFilePrefix + NET_OBJ_FILE_EXTENSION;
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
	
	private static void writeTreesToObjFile(ArrayList<PHYTree> trees, String[] sampleNames, Args args) {
		String treeFileName = args.outputFilePrefix + TREES_OBJ_FILE_EXTENSION;
		try {
			FileOutputStream fos = new FileOutputStream(treeFileName);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(sampleNames);
			int nSave = trees.size() >= args.numSave ? args.numSave : trees.size(); 
			oos.writeObject(new ArrayList<PHYTree>(trees.subList(0, nSave)));
			oos.close();
			fos.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + treeFileName);
			System.exit(-1);
		}
	}
	
	private static PHYNetwork readNetworkFromObjFile(String netFileName) {
		try {
			FileInputStream fis = new FileInputStream(netFileName);
			ObjectInputStream ois = new ObjectInputStream(fis);
			PHYNetwork net = (PHYNetwork) ois.readObject();
			ois.close();
			fis.close();
			return net;
		} catch (IOException e) {
			System.err.println("Failed to read to the file: " + netFileName + "\nNote: The -build command (with -s option) has to be run before the -show command for the given file.");
			System.exit(1);
		} catch (ClassNotFoundException e) {
			System.err.println("Failed to load network from file: " + netFileName + "\nNote: The -build command (with -s option) has to be run before the -show command for the given file.");
			System.exit(1);
		}
		return null;
	}
	
	private static ArrayList<PHYTree> readTreesFromObjFile(String treeFileName, String[] outSampleNames) {
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
		options.addOption("build", false, "Construct the sample lineage trees");
		options.addOption("show", false, "Display the saved lineage tree(s)");
		options.addOption("i", true, "Input file path");
		options.addOption("o", true, "Output file path (default: input file path)");
		options.addOption("n", "normal", true, "Normal sample id (default: 0)");
		options.addOption("e", true, "VAF error margin (default: 0.1)");
		options.addOption("maxVAFValid", true, "Maximum allowed VAF in a sample (default: 0.6)");
		options.addOption("maxVAFAbsent", true, "Maximum VAF to robustly consider an SSNV as absent from a sample (default: 0.03)");
		options.addOption("minVAFPresent", true, "Minimum VAF to robustly consider an SSNV as present in a sample (default: 0.05)");
		options.addOption("minGroupSize", true, "Minimum SSNV group size (default: 1)");
		options.addOption("minRobustGroup", true, "Minimum number of robust SSNVs in group to be robust (default: 2)");
		options.addOption("minRobustGroupKeep", true, "Minimum number of robust SSNVs per group to keep the group in the network initially (default: 0 - keeps all the groups)");
		options.addOption("minTargetDistRatio", true, "Minimum non-robust SSNV to target group's SSNV ratio per sample (default: 0.5)");
		options.addOption("minClusterSize", true, "Minimum size a cluster must have to be a considered a node in the network (default: 2)");
		options.addOption("minPrivateClusterSize", true, "Minimum size a private mutation cluster must have to be a considered a node in the network (default: 1)");
		options.addOption("maxClusterDist", true, "Maximum mean VAF difference up to which two clusters can be collapsed (default: 0.2)");
		options.addOption("c", "completeNetwork", false, "Add all possible edges to the constraint network (default: private nodes are connected only to closest level parents; only nodes with no other parents are descendants of root)");
		options.addOption("cp", false, "Input data represents cell prevalaence (CP) values");
		options.addOption("s", "save", true, "Maximum number of output trees to save (default: 1)");
		options.addOption("tree", "showTree", false, "Display the top ranking lineage tree");
		options.addOption("net", "showNetwork", false, "Display the constraint network");
		options.addOption("top", true, "Number of top ranking lineage trees to display");
		options.addOption("v", "verbose", false, "Verbose mode");
		options.addOption("h", "help", false, "Print usage");
		//options.addOption("vcf", false, "Input data type is vcf (default: validation)");
		//options.addOption("cnv", true, "File path to CNV regions used to mark SNVs");
		//options.addOption("ann", true, "File path to ANNOVAR SNV annotations");
		//options.addOption("cosmic", true, "File path to COSMIC SNV annotations");
		//options.addOption("tcga", true, "File path to TCGA SNV annotations");
		
		// order
		ArrayList<Option> optionsList = new ArrayList<Option>();
		optionsList.add(options.getOption("build"));
		optionsList.add(options.getOption("show"));
		optionsList.add(options.getOption("i"));
		optionsList.add(options.getOption("o"));
		optionsList.add(options.getOption("n"));
		optionsList.add(options.getOption("e"));
		optionsList.add(options.getOption("maxVAFValid"));
		optionsList.add(options.getOption("maxVAFAbsent"));
		optionsList.add(options.getOption("minVAFPresent"));
		optionsList.add(options.getOption("minGroupSize"));
		optionsList.add(options.getOption("minRobustGroup"));
		optionsList.add(options.getOption("minRobustGroupKeep"));
		optionsList.add(options.getOption("minTargetDistRatio"));
		optionsList.add(options.getOption("minClusterSize"));
		optionsList.add(options.getOption("minPrivateClusterSize"));
		optionsList.add(options.getOption("maxClusterDist"));
		optionsList.add(options.getOption("c"));
		optionsList.add(options.getOption("cp"));
		optionsList.add(options.getOption("s"));
		optionsList.add(options.getOption("net"));
		optionsList.add(options.getOption("tree"));
		optionsList.add(options.getOption("top"));
		optionsList.add(options.getOption("v"));
		optionsList.add(options.getOption("h"));
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmdLine = null;
		HelpFormatter hf = new HelpFormatter();
		hf.setOptionComparator(new OptionComarator<Option>(optionsList));
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			hf.printHelp("lichee", options);
			System.exit(-1);
		}
		
		// Set-up input args
		Args params = new Args();	
		
		// required options
		if(!cmdLine.hasOption("i")) {
			System.out.println("Required parameter: input file path [-i]");
			hf.printHelp("lichee", options);
			System.exit(-1);
		}
		if(!cmdLine.hasOption("n")) {
			System.out.println("Required parameter: normal sample id [-n]");
			hf.printHelp("lichee", options);
			System.exit(-1);
		}
		
		if(cmdLine.hasOption("n")) {
			params.normalSampleId = Integer.parseInt(cmdLine.getOptionValue("n"));
		}
		if(cmdLine.hasOption("e")) {
			Parameters.AAF_ERROR_MARGIN = Double.parseDouble(cmdLine.getOptionValue("e"));
		}
		if(cmdLine.hasOption("maxVAFValid")) {
			Parameters.MAX_ALLOWED_VAF = Integer.parseInt(cmdLine.getOptionValue("maxVAFValid"));
		}
		if(cmdLine.hasOption("maxVAFAbsent")) {
			Parameters.VALIDATION_SOFT_THR = Double.parseDouble(cmdLine.getOptionValue("maxVAFAbsent"));
		}
		if(cmdLine.hasOption("minVAFPresent")) {
			Parameters.VALIDATION_THR = Double.parseDouble(cmdLine.getOptionValue("minVAFPresent"));
		}
		if(cmdLine.hasOption("minGroupSize")) {
			Parameters.GROUP_SIZE_THR = Integer.parseInt(cmdLine.getOptionValue("minGroupSize"));
		}
		if(cmdLine.hasOption("minRobustGroup")) {
			Parameters.ROBUSTGROUP_SIZE_THR = Integer.parseInt(cmdLine.getOptionValue("minRobustGroup"));
		}
		if(cmdLine.hasOption("minRobustGroupKeep")) {
			Parameters.GROUP_ROBUST_NUM_THR = Integer.parseInt(cmdLine.getOptionValue("minRobustGroupKeep"));
		}
		if(cmdLine.hasOption("minTargetDistRatio")) {
			Parameters.MIN_VAF_TARGET_RATIO_PER_SAMPLE = Double.parseDouble(cmdLine.getOptionValue("minTargetDistRatio"));
		}
		if(cmdLine.hasOption("minClusterSize")) {
			Parameters.MIN_CLUSTER_SIZE = Integer.parseInt(cmdLine.getOptionValue("minClusterSize"));
		}
		if(cmdLine.hasOption("minPrivateClusterSize")) {
			Parameters.MIN_PRIVATE_CLUSTER_SIZE = Integer.parseInt(cmdLine.getOptionValue("minPrivateClusterSize"));
		}
		if(cmdLine.hasOption("maxClusterDist")) {
			Parameters.MAX_COLLAPSE_CLUSTER_DIFF = Double.parseDouble(cmdLine.getOptionValue("maxClusterDist"));
		}
		if(cmdLine.hasOption("c")) {
			Parameters.ALL_EDGES = true;
		}
		if(cmdLine.hasOption("cp")) {
			Parameters.CP = true;
			Parameters.AAF_MAX = 1.0;
			Parameters.MAX_ALLOWED_VAF = 1.0;
		}
		if(cmdLine.hasOption("vcf")) {
			Parameters.INFORMAT = Parameters.Format.VCF;
		}
		if(cmdLine.hasOption("tree")) {
			params.showBest = true;
		}
		if(cmdLine.hasOption("net")) {
			params.showNetwork = true;
		}
		if(cmdLine.hasOption("top")) {
			params.numShow = Integer.parseInt(cmdLine.getOptionValue("top"));
			if(params.numShow <= 0) {
				params.showBest = false;
			}
		}
		if(cmdLine.hasOption("s")) {
			params.numSave = Integer.parseInt(cmdLine.getOptionValue("s")); 
			if(params.numSave > 0) {
				params.persist = true;
			}
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
			new HelpFormatter().printHelp("lichee", options);
			System.exit(-1);
		}
	}
	
	protected static class Args {
		// --- 'build' command ---
		String inputFileName;
		String outputFilePrefix;	
		int normalSampleId = 0;
		String cnvFileName;
		String annFileName;
		String cosmicFileName;
		String tcgaFileName;
		
		// --- 'show' command ---
		String showFileNamePrefix;
		int numShow = 1;
		int numSave = 1;
		
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
	
	protected static class OptionComarator<T extends Option> implements Comparator<T> {
	    protected ArrayList<Option> orderedOptions;
	    public OptionComarator(ArrayList<Option> options) {
	    	orderedOptions = options;
	    }
	    public int compare(T o1, T o2) {
	        return orderedOptions.indexOf(o1) - orderedOptions.indexOf(o2);
	    }
	}
}
