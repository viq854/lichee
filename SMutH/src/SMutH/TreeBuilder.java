package SMutH;

import io.*;


import java.io.*;
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
public class TreeBuilder {
	
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
		LPM029 (3)
		;
		
		
		
		
		public final int normal;
		Samples(int normal) {
	        this.normal = normal;
	    }
		
	}
	
	public static String path;
	public static String testName;
	public static int normalSample;
			
	private static CommandLine cmdLineArgs;
	
	private DirectedGraph<Integer, Integer> graph;
	
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
		cmdLineArgs = setOptions(args);
		//System.out.println(cmdLineArgs.getArgs()[0]);
		args = cmdLineArgs.getArgs();
		//System.out.println(args[0]);
		/*
		 * For updated input:
		 * 1. Build test cases
		 * 2. Change build configurations
		 * 3. Update line by line processing
		 * 4. Set columns/row sizes early(?)
		 * Looks like it works for now. Can go back later.
		 * 
		 */
		/*if (args.length == 0) {
			System.out.println("Must pass input matrix!");
			System.exit(1);
		}
		*/

		//for (Samples s: Samples.values())
		{
		testName = "tree_5_01";//s.toString();
		normalSample = 0;//s.normal - 1;
		//path =  "/Users/rahelehs/Work/ash/"+testName+"/"
		//path =  "/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/"+testName+"/";
		path =  "/Users/rahelehs/Work/cancerTree/simulation_vcfs/RECOMB2013/";
		String inputFile = path+testName+".raw.vcf";
		SNVDatabase vcfDB = new SNVDatabase(inputFile, normalSample);
		vcfDB.generateGATKFile(path+testName + ".GATK-output.txt");
		ArrayList<ArrayList<Integer>> matrixPrime = TreeChecker.checkIfTree(vcfDB);
		if (matrixPrime != null) {
			System.out.println("This can be a PhyTree!");
			TreeChecker.printMatrix(matrixPrime);
			TreeBuilder tb = new TreeBuilder();
			HashMap<Integer, Integer> LColFuncMap = tb.getLColFuncMap(matrixPrime);
			Map<String, Integer> mutMap = TreeChecker.getMutMap();
			vcfDB.printSNVs(path+testName+".validSNVs.txt");
			//System.out.println(LColFuncMap.toString());
			tb.buildTree(matrixPrime, mutMap, LColFuncMap,vcfDB);
		} else {
			System.out.println("This cannot be a PhyTree!");
			matrixPrime = TreeChecker.getMatrixPrime(vcfDB);
			ArrayList<ArrayList<Integer>> noConflictMatrixPrime = TreeChecker.getCFMatrixPrime(vcfDB);
			Set<ArrayList<Integer>> conflicts = TreeChecker.getConflicts(matrixPrime, noConflictMatrixPrime);
			Map<String, Integer> mutMap = TreeChecker.getMutMap();
			/*
			 * -->mutMap update - Add in all 1's as legit code
			 * -->Edit SNVs
			 * -->Subpopulation Replacement
			 */
			vcfDB.editSNV(conflicts);
			vcfDB.printSNVs(path+testName+".validSNVs.txt");

			SPBuilder spb = new SPBuilder(mutMap, noConflictMatrixPrime, conflicts, testName, vcfDB);
			TreeBuilder tb = new TreeBuilder();	
			tb.graph = spb.getTree();
				//HashMap<Integer, Integer> LColFuncMap = tb.getLColFuncMap(noConflictMatrixPrime);
				//tb.graph = tb.buildNetwork(matrixPrime, mutMap, LColFuncMap,vcfDB);
				
				/**/
			new TreeVisualizer(tb.graph, spb.getNodeLabels(),spb.getEdgeLabels(),vcfDB);
				
//				TreeBuilder tb = new TreeBuilder();
//				HashMap<Integer, Integer> LColFuncMap = tb.getLColFuncMap(noConflictMatrixPrime);
//				System.out.println(LColFuncMap.toString());
//				tb.buildTree(matrixPrime, LColFuncMap);

			
			/*
			 * For conflicts: 
			 * 1. Get conflict-free matrix - DONE
			 * 2. Get conflict columns - Done
			 * 2b. Get codes for conflict columns. - Don't need. Can just look up columns
			 * with matrixPrime
			 * 3. Place in tree 
			 */
//			TreeBuilder tb = new TreeBuilder();
//			HashMap<Integer, Integer> LColFuncMap = tb.getLColFuncMap(noConflictMatrixPrime);
//			tb.buildTree(noConflictMatrixPrime, LColFuncMap);
		}
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
	 * Function: getPossibleCodes(String conflictStr, int k)
	 * Usage: Set<String> codes = getPossibleCodes(conflictStr, k)
	 * ----
	 * A wrapper for recCodes.
	 * 
	 * @param conflictStr The string for which possible codes are being generated.
	 * @param k	Edit distance (CURRENTLY ONLY WORKS IF == 1)
	 * @return	A set of possible strings which conflictStr can turn into
	 */
	/*private static Set<String> getPossibleCodes(String conflictStr, int k) {
		Set<String> codes = new HashSet<String>();
		if (k > conflictStr.length()) k = conflictStr.length();
		recCodes(conflictStr, k, codes, new ArrayList<Integer>());
		return codes;
	}
*/
	/**
	 * Function: recCodes(String conflictStr, int k, Set<String> codes, ArrayList<Integer> usedIndices)
	 * Usage: recCodes(conflictStr, k, codes, usedIndices)
	 * ----
	 * Recursively generates all possible codes conflictStr could become with edit distance == k.
	 * The codes are stored in the set passed into as the third argument. usedIndices keep track
	 * of which indices have been changed previously.
	 * @param conflictStr	The string for which we are getting changes
	 * @param k				The edit distance
	 * @param codes			The set of string codes which we store
	 * @param usedIndices	An ArrayList of integers which store changed indices.
	 */
	/*private static void recCodes(String conflictStr, int k, Set<String> codes, ArrayList<Integer> usedIndices) {
		if (k == 0){
			codes.add(conflictStr);
		} else {
			for (int i = 0; i < conflictStr.length(); i++){
				if (!usedIndices.contains(i)){
					char newChar = (conflictStr.charAt(i) == '0' ? '1' : '0');
					ArrayList<Integer> newIndices = new ArrayList<Integer>(usedIndices);
					newIndices.add(i);
					if (i == 0){
						recCodes(newChar + conflictStr.substring(1), k-1, codes, newIndices);
					} else if (i == conflictStr.length() - 1){
						recCodes(conflictStr.substring(0, i) + newChar, k-1, codes, newIndices);
					} else {
						recCodes(conflictStr.substring(0, i) + newChar + conflictStr.substring(i+1),k-1,
								codes, newIndices);
					}
				}
			}
		}
	}
*/


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
	private void buildTree(ArrayList<ArrayList<Integer>> matrixPrime, Map<String, Integer> mutMap,
			HashMap<Integer, Integer> lColFuncMap, SNVDatabase db) {
		DirectedGraph<Integer, Integer> g = new DirectedSparseGraph<Integer, Integer>();
		int numRows = matrixPrime.size();
		int numCols = matrixPrime.get(0).size();
		HashMap<Integer, String> edgeLabels = new HashMap<Integer, String>();
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
		for (int i = 0; i < numCols; i++){
			String edgeStr = "";
			for (int j = 0; j < numRows ; j++)
				edgeStr += matrixPrime.get(j).get(i);
			if (mutMap.get(edgeStr).intValue() > 0)
				edgeLabels.put(new Integer(i+1), mutMap.get(edgeStr).toString());
		}
		graph = g;
		new TreeVisualizer(g, null,edgeLabels,db);
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
	public HashMap<Integer, Integer> getLColFuncMap(
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