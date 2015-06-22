/*
 * Program LICHeE for multi-sample cancer phylogeny reconstruction
 * by Victoria Popic (viq@stanford.edu) 2014
 *
 * MIT License
 *
 * Copyright (c) 2014 Victoria Popic.
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/


package lineage;

import java.io.FileWriter;
import java.io.IOException;
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
import lineage.Parameters.Format;
import lineage.PHYTree;

/**
 * Main cell lineage reconstruction pipeline
 * @autor viq
 */
public class LineageEngine {
	protected static final Logger logger = Logger.getLogger("lineage.engine");
	
	/**
	 * The main pipeline for reconstructing the cell lineage trees
	 */
	public static void buildLineage(Args args) {
				
		// 1. load SNV data
		SNVDataStore db = new SNVDataStore(args.inputFileName, args.clustersFileName, args.normalSampleId);
		
		// 2. get the SNVs partitioned by group tag and create the appropriate SNV group objects
		HashMap<String, ArrayList<SNVEntry>> snvsByTag = db.getSomaticGroups();
		ArrayList<SNVGroup> groups = new ArrayList<SNVGroup>();
		for(String groupTag : snvsByTag.keySet()) {
			groups.add(new SNVGroup(groupTag, snvsByTag.get(groupTag), db.isRobustGroup(groupTag)));
		}
		if(groups.size() == 0) {
			logger.warning("All SNV groups have been filtered out.");
			return;
		}
		
		// 3. cluster SNVs in each group
		AAFClusterer clusterer = new AAFClusterer();
		for(SNVGroup group : groups) {
			if(args.clustersFileName == null) {
				Cluster[] clusters = clusterer.clusterSubPopulations(group, ClusteringAlgorithms.EM, 1);
				logger.fine("Clustering results for group: " + group.getTag());
				for(Cluster c : clusters) {
					logger.fine(c.toString());
				}
				group.setSubPopulations(clusters);
			} else {
				ArrayList<Cluster> groupClusters = db.getClusters().get(group.getTag());
				group.subPopulations = new Cluster[groupClusters.size()];
				group.subPopulations = groupClusters.toArray(group.subPopulations);
			}
		}
		
		// 4. construct the constraint network
		PHYNetwork constrNetwork = new PHYNetwork(groups, db.getNumSamples());
		logger.fine(constrNetwork.toString());
		
		// 5. find all the lineage trees that pass the VAF constraints
		ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees();  
		logger.info("Found " + spanningTrees.size() + " valid tree(s)");
		
		if(spanningTrees.size() == 0) {
			logger.info("Adjusting the network...");	
			// if no valid trees were found, fix the network (e.g. remove group nodes that are not robust)
			int delta = 0;
			do {
				int numNodes = constrNetwork.numNodes;
				constrNetwork = constrNetwork.fixNetwork();
				spanningTrees = constrNetwork.getLineageTrees();  
				delta = numNodes - constrNetwork.numNodes; 
			} while((delta != 0) && (spanningTrees.size() <= 0));
			logger.info("Found " + spanningTrees.size() + " valid trees after network adjustments");	
		}
		
		// 6. evaluate/rank the trees
		if(spanningTrees.size() > 0) {
			constrNetwork.evaluateLineageTrees();
			logger.fine("Top tree\nError score: " + spanningTrees.get(0).getErrorScore());	
			logger.fine(spanningTrees.get(0).toString());
		} 
		
		// 7. result visualization
		if(args.showNetwork) {
			constrNetwork.displayNetwork();
		}
		if(spanningTrees.size() > 0) {
			for(int i = 0; i < args.numShow; i++) {
				if(spanningTrees.size() > i) {
					constrNetwork.displayTree(spanningTrees.get(i), db.getSampleNames(), null, null);
				} else {
					break;
				}
			}
			// 8. persistent storage
			if(args.numSave > 0) {
				writeTreesToTxtFile(constrNetwork, spanningTrees, db.getSampleNames(), args);
			}	
		} 
	}
	
	///// I/O /////
	
	private static void writeTreesToTxtFile(PHYNetwork net, ArrayList<PHYTree> trees, ArrayList<String> sampleNames, Args args) {
		String treeFileName = args.outputFileName;
		try {
			FileWriter fw = new FileWriter(treeFileName);
			fw.write("Nodes:\n" + net.getNodesWithMembersAsString() + "\n");
			for(int i = 0; i < args.numSave; i++) {
				if(trees.size() > i) {
					fw.write("****Tree " + i + "****\n");
					String edges = trees.get(i).toString();
					fw.write(edges);
					fw.write("Error score: " + trees.get(i).getErrorScore()+"\n\n");	
					fw.write("Sample decomposition: \n");
					String lineage = "";
					for(int j = 0; j < sampleNames.size(); j++) {
						lineage += trees.get(i).getLineage(j, sampleNames.get(j));
						lineage += "\n";
					}
					fw.write(lineage);
					fw.write("\n");
				}
			}
			fw.write("SNV info:\n" + net.getNodeMembersOnlyAsString() + "\n");
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + treeFileName);
			System.exit(-1);
		}
	}
	
	// ---- LAUNCH ----
	
	private static final String TREES_TXT_FILE_EXTENSION = ".trees.txt";
	public static void main(String[] args) {
		Options options = new Options(); 
		// Commands
		options.addOption("build", false, "Construct the sample lineage trees");
		
		// Input/Output/Display
		options.addOption("i", true, "Input file path [required]");
		options.addOption("o", true, "Output file path (default: input file with suffix .trees.txt)");
		options.addOption("cp", false, "Input data represents cell prevalaence (CP) values");
		options.addOption("sampleProfile", false, "Input file contains the SSNV sample presence-absence profile (this will disable the default SSNV calling step)");
		options.addOption("n", "normal", true, "Normal sample column id in the list of samples, 0-based (e.g 0 is the first column) [required without -sampleProfile]");
		options.addOption("clustersFile", true, "SSNV clusters file path");
		options.addOption("s", "save", true, "Maximum number of output trees to save (default: 1)");
		options.addOption("showNetwork", "net", false, "Display the constraint network");
		options.addOption("showTree", "tree", true, "Number of top-ranking trees to display (default: 0)");
	
		// SSNV filtering / calling
		options.addOption("maxVAFAbsent", "absent", true, "Maximum VAF to robustly consider an SSNV as absent from a sample [required without -sampleProfile]");
		options.addOption("minVAFPresent", "present", true, "Minimum VAF to robustly consider an SSNV as present in a sample [required without -sampleProfile]");
		options.addOption("maxVAFValid", true, "Maximum allowed VAF in a sample (default: 0.6)");
		options.addOption("minProfileSupport", true, "Minimum number of robust SSNVs required for a group presence-absence profile to be labeled robust (default: 2)");
		
		// Network Construction / Tree Search
		options.addOption("minClusterSize", true, "Minimum size a cluster must have to be a considered a node in the network (default: 2)");
		options.addOption("minPrivateClusterSize", true, "Minimum size a private mutation cluster must have to be a considered a node in the network (default: 1)");
		options.addOption("minRobustNodeSupport", true, "Minimum number of robust SSNVs required for a node to be labeled robust during tree search: non-robust nodes can be removed from the network when no valid lineage trees are found (default: 2)");
		options.addOption("maxClusterDist", true, "Maximum mean VAF difference up to which two clusters can be collapsed (default: 0.2)");
		options.addOption("c", "completeNetwork", false, "Add all possible edges to the constraint network (default: private nodes are connected only to closest level parents; only nodes with no other parents are descendants of root)");
		options.addOption("e", true, "VAF error margin (default: 0.1)");
		options.addOption("nTreeQPCheck", true, "Number of top-ranking trees on which the QP consistency check is run, we have not seen this check fail in practice (default: 0, for best performance)");
		
		options.addOption("v", "verbose", false, "Verbose mode");
		options.addOption("h", "help", false, "Print usage");
		
		// display order
		ArrayList<Option> optionsList = new ArrayList<Option>();
		optionsList.add(options.getOption("build"));

		optionsList.add(options.getOption("i"));
		optionsList.add(options.getOption("o"));
		optionsList.add(options.getOption("cp"));
		optionsList.add(options.getOption("sampleProfile"));
		optionsList.add(options.getOption("n"));
		optionsList.add(options.getOption("clustersFile"));
		optionsList.add(options.getOption("s"));
		optionsList.add(options.getOption("net"));
		optionsList.add(options.getOption("tree"));
		optionsList.add(options.getOption("maxVAFAbsent"));
		optionsList.add(options.getOption("minVAFPresent"));
		optionsList.add(options.getOption("maxVAFValid"));
		optionsList.add(options.getOption("minProfileSupport"));
		optionsList.add(options.getOption("minClusterSize"));
		optionsList.add(options.getOption("minPrivateClusterSize"));
		optionsList.add(options.getOption("minRobustNodeSupport"));
		optionsList.add(options.getOption("maxClusterDist"));
		optionsList.add(options.getOption("c"));
		optionsList.add(options.getOption("e"));
		optionsList.add(options.getOption("nTreeQPCheck"));
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
		if(cmdLine.hasOption("i")) {
			params.inputFileName = cmdLine.getOptionValue("i");
		} else {
			System.out.println("Required parameter: input file path [-i]");
			hf.printHelp("lichee", options);
			System.exit(-1);
		}
		if(cmdLine.hasOption("o")) {
			params.outputFileName = cmdLine.getOptionValue("o");	
		} else {
			params.outputFileName = params.inputFileName + TREES_TXT_FILE_EXTENSION;
		}
		if(cmdLine.hasOption("clustersFile")) {
			params.clustersFileName = cmdLine.getOptionValue("clustersFile");
		}
		if(cmdLine.hasOption("sampleProfile")) {
			Parameters.INPUT_FORMAT = Format.SNV_WITH_PROFILE;
		}
		
		if(cmdLine.hasOption("n")) {
			params.normalSampleId = Integer.parseInt(cmdLine.getOptionValue("n"));
		} else if(!cmdLine.hasOption("sampleProfile")) {
			System.out.println("Required parameter: normal sample id [-n]");
			hf.printHelp("lichee", options);
			System.exit(-1);
		}	
		if(cmdLine.hasOption("showTree")) {
			params.numShow = Integer.parseInt(cmdLine.getOptionValue("showTree"));
		}
		if(cmdLine.hasOption("showNetwork")) {
			params.showNetwork = true;
		}	
		if(cmdLine.hasOption("s")) {
			params.numSave = Integer.parseInt(cmdLine.getOptionValue("s")); 
		}
		
		if(cmdLine.hasOption("maxVAFAbsent")) {
			Parameters.MAX_VAF_ABSENT = Double.parseDouble(cmdLine.getOptionValue("maxVAFAbsent"));
		} else if(!cmdLine.hasOption("sampleProfile")) {
			System.out.println("Required parameter: -maxVAFAbsent");
			hf.printHelp("lichee", options);
			System.exit(-1);
		}
		if(cmdLine.hasOption("minVAFPresent")) {
			Parameters.MIN_VAF_PRESENT = Double.parseDouble(cmdLine.getOptionValue("minVAFPresent"));
		} else if(!cmdLine.hasOption("sampleProfile")) {
			System.out.println("Required parameter: -minVAFPresent");
			hf.printHelp("lichee", options);
			System.exit(-1);
		}
		if(cmdLine.hasOption("maxVAFValid")) {
			Parameters.MAX_ALLOWED_VAF = Integer.parseInt(cmdLine.getOptionValue("maxVAFValid"));
		}
		if(cmdLine.hasOption("minProfileSupport")) {
			Parameters.MIN_GROUP_PROFILE_SUPPORT = Integer.parseInt(cmdLine.getOptionValue("minProfileSupport"));
		}
		if(cmdLine.hasOption("minClusterSize")) {
			Parameters.MIN_CLUSTER_SIZE = Integer.parseInt(cmdLine.getOptionValue("minClusterSize"));
		}
		if(cmdLine.hasOption("minPrivateClusterSize")) {
			Parameters.MIN_PRIVATE_CLUSTER_SIZE = Integer.parseInt(cmdLine.getOptionValue("minPrivateClusterSize"));
		}
		if(cmdLine.hasOption("minRobustNodeSupport")) {
			Parameters.MIN_ROBUST_CLUSTER_SUPPORT = Integer.parseInt(cmdLine.getOptionValue("minRobustNodeSupport"));
		}
		if(cmdLine.hasOption("maxClusterDist")) {
			Parameters.MAX_COLLAPSE_CLUSTER_DIFF = Double.parseDouble(cmdLine.getOptionValue("maxClusterDist"));
		}
		if(cmdLine.hasOption("c")) {
			Parameters.ALL_EDGES = true;
		}
		if(cmdLine.hasOption("cp")) {
			Parameters.CP = true;
			Parameters.VAF_MAX = 1.0;
			Parameters.MAX_ALLOWED_VAF = 1.0;
		}
		if(cmdLine.hasOption("e")) {
			Parameters.VAF_ERROR_MARGIN = Double.parseDouble(cmdLine.getOptionValue("e"));
		}
		if(cmdLine.hasOption("nTreeQPCheck")) {
			Parameters.NUM_TREES_FOR_CONSISTENCY_CHECK = Integer.parseInt(cmdLine.getOptionValue("nTreeQPCheck"));
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
			buildLineage(params);
			
		} else {
			new HelpFormatter().printHelp("lichee", options);
			System.exit(-1);
		}
	}
	
	protected static class Args {
		// --- 'build' command ---
		String inputFileName;
		String outputFileName;	
		String clustersFileName;	
		int normalSampleId = 0;
		String cnvFileName;
		String annFileName;
		String cosmicFileName;
		String tcgaFileName;
		
		// --- 'show' command ---
		String showFileNamePrefix;
		int numShow = 0;
		int numSave = 1;
		
		// flags
		boolean showNetwork = false;
		boolean verbose = false;
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
