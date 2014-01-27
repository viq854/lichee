package ppt;

import util.*;

import java.text.DecimalFormat;
import java.util.*;

import org.apache.commons.cli.*;

import edu.uci.ics.jung.graph.*;
import edu.uci.ics.jung.graph.util.*;

/**
 * Class: TreeBuilder
 * Constructor: (None)
 * Last Edited: August 13, 2013
 * ----
 * This class is the main program. It calls other classes to build
 * the tree from the initial input. This input should be a
 * VCF file.
 */
public class PerfectPhylogenyEngine {
	
	private static enum Samples {
		LPJ040 (3),
		LPJ041 (2),
		LPJ128 (2),
		LPM001 (2),	
		LPM005 (2),	
		LPM011 (2),	
		LPM012 (2),	
		LPM013 (2),	
		LPM016 (2),	
		LPM017 (2),	
		LPM018 (2),	
		LPM019 (2),	
		LPM020 (2),	
		LPM021 (2),	
		LPM022 (2),	
		LPM023 (2),	
		LPM024 (2),	
		LPM025 (2),	
		LPM026 (3),	
		LPM027 (2),	
		LPM028 (2),	
		LPM029 (3);
		
		public final int normal;
		Samples(int normal) {
	        this.normal = normal;
	    }
		
	}
			
	/**
	 * Function: main(String[] args)
	 * Usage: (Main Method)
	 * ----
	 * Main method
	 * 
	 * Essentially, a wrapper function which first checks if building a tree is possible,
	 * then generates the tree.
	 * 
	 * @param args	Main's normal args parameter
	 */
	public static void main(String[] args) {
		CommandLine cmdLineArgs = setOptions(args);
		args = cmdLineArgs.getArgs();

		/*if (args.length == 0) {
			System.out.println("Must pass input matrix!");
			System.exit(1);
		}*/
		

		//for (Samples s: Samples.values())
		{
		Configs.testName = "Patient_2";//s.toString();
		int normalSample = 0;
		//Configs.path =  "/Users/rahelehs/Work/ash/AllMutations/";
		Configs.path =  "/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/"+Configs.testName+"/";
		//path =  "/Users/rahelehs/Work/cancerTree/simulation_vcfs/RECOMB2013/";
		String inputFile = Configs.path+Configs.testName+".validation.txt";
		SNVDatabase snvDB = new SNVDatabase(inputFile, normalSample);
		snvDB.resolveNonRobustconflicts();

		TreeChecker tc = new TreeChecker(snvDB);
		if (tc.isPhyTree()) {
			System.out.println("This is a Perfect Phylogeny Tree.");
			HashMap<Integer, Integer> LColFuncMap = getLColFuncMap(tc.getMatrixPrime());
			buildTree(tc.getMatrixPrime(), LColFuncMap,snvDB);
		} else {
			System.out.println("This is not a Perfect Phylogeny Tree!");
			
			/*
			 * Since there is no Perfect Phylogeny there are some conflict groups, 
			 * conflict resolution 
			 * 1. find minimum number of conflicting mutations
			 * 2. edit them to non-conflicting groups
			 * 3. find sub-population ### this has been removed currently!!!!
			 * 4. build the conflict free tree
			 */
			
		
			
			Set<String> conflicts = tc.getConflicts();
			
			snvDB.resolveConflicts(conflicts);
			System.out.println("Mutation Map after editing to perfect phytree groups!");
			snvDB.reportGroups();
			
			//Step 3 ####
			//
            //SPBuilder spb = new SPBuilder(mutMap, conflictFreeMatrixPrime, conflicts, vcfDB);

			//Step 4
			ArrayList<ArrayList<Integer>> conflictFreeMatrixPrime = tc.getCFMatrixPrime();	
			System.out.println("Conflict Free Matrix Prime:"); 

			TreeChecker.printMatrix(conflictFreeMatrixPrime);
			HashMap<Integer, Integer> LColFuncMap = getLColFuncMap(conflictFreeMatrixPrime);
			buildTree(conflictFreeMatrixPrime, LColFuncMap,snvDB);

		}
		
		snvDB.printSNVs(Configs.path+Configs.testName+".validSNVs.txt");
		
		
		}
		
	}


	/**
	 * Function: setOptions(String[] args)
	 * Usage: cmdLineArgs = setOptions(String[] args)
	 * ----
	 * @param args	Arguments passed into the command line
	 * @return	a CommandLine object specifying 
	 */
	private static CommandLine setOptions(String[] args) {
		Options options = new Options();
		options.addOption("a", "altInput", false, 
				"Take in an alternate transposed matrix with rows being objects/mutations");
		CommandLineParser parser = new PosixParser();
		CommandLine cmd = null;
		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return cmd;
	}


	/**
	 * Builds the PhyTree when given M' and the LColFuncMap
	 * 
	 * Uses Gusfield's 1991 phylogenetic tree algorithm to construct
	 * a PhyTree. The graph data structures are taken from the JUNG
	 * library, a graph package for Java.
	 * 
	 * @param matrixPrime	M'
	 * @param lColFuncMap	The LColFuncMap made from getLColFuncMap
	 */
	private static void buildTree(ArrayList<ArrayList<Integer>> matrixPrime,
			HashMap<Integer, Integer> lColFuncMap, SNVDatabase db) {
		DirectedGraph<Integer, Integer> g = new DirectedSparseGraph<Integer, Integer>();
		
		int numRows = matrixPrime.size();
		int numCols = matrixPrime.get(0).size();
		HashMap<Integer, String> nodeLabels = new HashMap<Integer, String>();
		
		 Map<String, Integer> mutMap = db.getTAG2SNVNum();
		//Add root
		g.addVertex((Integer) 0);
		for (int i = 0; i < numCols; i++){
			Integer nodeNum = i + 1;
			g.addVertex(nodeNum);
			if (lColFuncMap.get(i + 1) > 0){
				g.addEdge(i + 1, lColFuncMap.get(nodeNum), nodeNum, EdgeType.DIRECTED);
			} else {
				g.addEdge(i + 1, 0, nodeNum, EdgeType.DIRECTED);
			}
		}
		for (int i = 0; i < numRows; i++){
			int maxIndex = 0;
				
			for (int j = numCols - 1; j >= 0; j--){
				if (matrixPrime.get(i).get(j) == 1){
					maxIndex = j + 1;
					break;
				}
				
				
				
			}
			Integer currParent = g.getDest(maxIndex);
			if (maxIndex == 0) currParent = 0;
			g.addVertex(-1 * (i + 1));
			g.addEdge(-1 * (i + 1), currParent, -1 * (i + 1), EdgeType.DIRECTED);
			
		}
		DecimalFormat df = new DecimalFormat("#.##");
		
		for (int i = 0; i < numCols; i++){
			String edgeStr = "";
			for (int j = 0; j < numRows ; j++)
				edgeStr += matrixPrime.get(j).get(i);
			if (mutMap.get(edgeStr).intValue() > 0){
				/*
				double [] m =db.getMeanAAF(edgeStr);
				String avgAAF = ""+df.format(m[0]);
				for (int j = 1; j < numRows ; j++)
					avgAAF += ","+df.format(m[j]);*/
					
				nodeLabels.put(new Integer(i+1), mutMap.get(edgeStr).toString());
			}
		}
		nodeLabels.put(0, "germline");
		for(int i=0; i < db.getNumofSamples(); i++){
			nodeLabels.put(-i-1, db.getName(i));
		}
		
		Visualizer.TreeVisualizer(g, nodeLabels);
	}

	
	/**
	 * Creates the LColFuncMap which maps the column to the largest L-function value 
	 * in the column
	 * 
	 * This function is different from TreeChecker's getLColFuncMap because it 
	 * assumes the input matrix is a legitimate tree. Thus, we can just fill in the
	 * L-values and find each column's max.
	 * 
	 * @param matrixPrime	M'
	 * @return				A HashMap of column indices to its respective greatest L-value
	 */
	public static HashMap<Integer, Integer> getLColFuncMap(
			ArrayList<ArrayList<Integer>> matrixPrime) {
		HashMap<Integer, Integer> colFuncMap = new HashMap<Integer, Integer>();
		ArrayList<ArrayList<Integer>> matCopy = new ArrayList<ArrayList<Integer>>();
		int numRows = matrixPrime.size();
		int numCols = matrixPrime.get(0).size();
		for (int i = 0; i < numRows; i++){
			int counter = 0;
			matCopy.add(new ArrayList<Integer>());
			for (int j = 0; j < numCols; j++){
				if (matrixPrime.get(i).get(j) != 0){
					matCopy.get(i).add(counter);
					counter = j + 1;
				} else matCopy.get(i).add(0);
			}
		}
		for (int j = 0; j < numCols; j++){
			int colValue = 0;
			for (int i = 0; i < numRows; i++){
				int cellValue = matCopy.get(i).get(j);
				if (cellValue != 0){
					colValue = cellValue;
					break;
				}
			}
			colFuncMap.put(j + 1, colValue);
		}
		return colFuncMap;
	}
	
}