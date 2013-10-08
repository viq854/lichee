package SMutH;

import io.VCFConstants;
import io.VCFDatabase;
import io.VCFEntry;

import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;

import javax.imageio.ImageIO;
import javax.swing.*;

import org.apache.commons.cli.*;
import org.apache.commons.collections15.Transformer;


import edu.uci.ics.jung.algorithms.layout.*;
import edu.uci.ics.jung.graph.*;
import edu.uci.ics.jung.graph.util.*;
import edu.uci.ics.jung.visualization.*;
import edu.uci.ics.jung.visualization.control.*;
import edu.uci.ics.jung.visualization.decorators.*;
import edu.uci.ics.jung.visualization.layout.LayoutTransition;
import edu.uci.ics.jung.visualization.util.Animator;

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
	
	private VisualizationViewer<Integer, Integer> visServer;
	private TreeLayout<Integer, Integer> treeLayout;
	private RadialTreeLayout<Integer, Integer> radialLayout;
	private VisualizationServer.Paintable rings;
	
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
		testName = "tree_4_11";//s.toString();
		normalSample = 0;//s.normal - 1;
		//path =  "/Users/rahelehs/Work/ash/"+testName+"/"
		//path =  "/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/"+testName+"/";
		path =  "/Users/rahelehs/Work/cancerTree/simulation_vcfs/";
		String inputFile = path+testName+".raw.vcf";
		VCFDatabase vcfDB = new VCFDatabase(inputFile, normalSample);
		vcfDB.generateMatrix("output.txt");
		vcfDB.generateGATKFile(path+testName + ".GATK-output.txt");
		ArrayList<ArrayList<Integer>> matrixPrime = TreeChecker.checkIfTree("output.txt");
		if (matrixPrime != null) {
			System.out.println("This can be a PhyTree!");
			TreeChecker.printMatrix(matrixPrime);
			TreeBuilder tb = new TreeBuilder();
			HashMap<Integer, Integer> LColFuncMap = tb.getLColFuncMap(matrixPrime);
			Map<String, Integer> mutMap = TreeChecker.getMutMap();
			printSNVs(mutMap,vcfDB);
			//System.out.println(LColFuncMap.toString());
			tb.buildTree(matrixPrime, mutMap, LColFuncMap,vcfDB);
		} else {
			System.out.println("This cannot be a PhyTree!");
			matrixPrime = TreeChecker.getMatrixPrime("output.txt");
			ArrayList<ArrayList<Integer>> noConflictMatrixPrime = TreeChecker.getCFMatrixPrime("output.txt");
			Set<ArrayList<Integer>> conflicts = TreeChecker.getConflicts(matrixPrime, noConflictMatrixPrime);
			Map<String, Integer> mutMap = TreeChecker.getMutMap();
			updateMutMap(mutMap, matrixPrime.size());
			/*
			 * while tree not satisfiable
			 * -->mutMap update - Add in all 1's as legit code
			 * -->Edit SNVs
			 * -->Subpopulation Replacement
			 */
			for (int i = 0; i < 1; i++){
				editSNV(conflicts, mutMap, vcfDB, i);
				SPBuilder spb = new SPBuilder(mutMap, noConflictMatrixPrime, conflicts, testName, vcfDB);
				TreeBuilder tb = new TreeBuilder();	
				tb.graph = spb.getTree();
				tb.VisualizeTree(tb.graph, spb.getNodeLabels(),spb.getEdgeLabels(),vcfDB);
				
//				TreeBuilder tb = new TreeBuilder();
//				HashMap<Integer, Integer> LColFuncMap = tb.getLColFuncMap(noConflictMatrixPrime);
//				System.out.println(LColFuncMap.toString());
//				tb.buildTree(matrixPrime, LColFuncMap);
			}
			
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

	private static void updateMutMap(Map<String, Integer> mutMap, int size) {
		String all1s = "";
		for (int i = 0; i < size; i++) all1s += "1";
		mutMap.put(all1s, 0);
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
	
	private static void printSNVs(Map<String, Integer> mutMap, VCFDatabase vcfDB) {
		ArrayList<String> codes = new ArrayList<String>(mutMap.keySet());
		Collections.sort(codes);
		Collections.reverse(codes);
		PrintWriter pwv = null, pwh=null;
		try{
			/*pwh = new PrintWriter(new FileWriter(path+testName+".LOH.txt"));
			vcfDB.printLOH(pwh);
			pwh.close();
			*/
			pwv = new PrintWriter(new FileWriter(path+testName+".validSNVs.txt"));
			pwv.write(vcfDB.getHeader()+"\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
		for (int i = 0; i < codes.size(); i++){
			if (mutMap.get(codes.get(i)) == 0)
				continue;
			System.out.println(codes.get(i) + ": " + mutMap.get(codes.get(i)));
			try {
				vcfDB.printEntriesByGATK(pwv, codes.get(i));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}	
		pwv.close();
	}
	
	/**
	 * Function: editSNV(Set<ArrayList<Integer>> conflicts, Map<String, Double> mutMap, VCFDatabase vcfDB)
	 * Usage: editSNV(conflicts, mutMap, vcfDB)
	 * ----
	 * @param conflicts	The set of conflicting binary codes in the original music
	 * @param mutMap	The map of codes to number of mutations with that GATK code
	 * @param vcfDB		The VCFDatabse corresponding to the VCF file for this run of matrix building
	 * @param iterCounter Which iteration of editSNV this is
	 * @param testName 
	 */
	private static void editSNV(Set<ArrayList<Integer>> conflicts, Map<String, Integer> mutMap, VCFDatabase vcfDB, int iterCounter) {
		System.out.println("Original Mutation Map");
		printSNVs(mutMap,vcfDB);
		ArrayList<String> codes = new ArrayList<String>(mutMap.keySet());
		Collections.sort(codes);
		Collections.reverse(codes);

		PrintWriter pw = null;
		try{
			pw = new PrintWriter(new FileWriter(path+testName+".edit.txt"));
			pw.write("#CHROM\tPOS\tOriginal Code\tNew Code\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
		//System.out.append("CONFLICTS: "+conflicts.toString());
		for (ArrayList<Integer> conflict: conflicts){
			String conflictStr = "";
			for (int i = 0; i < conflict.size(); i++){
				conflictStr += ("" + conflict.get(i));
			}
			if (conflictStr.equals("111001001")){
				int i = 1;
				i++;
			}
			ArrayList<VCFEntry> failCodes = null;
			//System.out.println("----");
			//System.out.println(conflictStr + ": " + mutMap.get(conflictStr));
			//Set<String> possible = getPossibleCodes(conflictStr, 1);
			//System.out.println("--POSSIBLE--");
			//System.out.println(possible.toString());
			Set<String> allCodes = new HashSet<String>(mutMap.keySet());
			//Remove all conflicts
			for (ArrayList<Integer> currConflict: conflicts){
				String currConflictStr = "";
				for (int i = 0; i < currConflict.size(); i++){
					currConflictStr += ("" + currConflict.get(i));
				}
				allCodes.remove(currConflictStr);
			}
			//--
			////allCodes.remove(conflictStr);
			////System.out.println(allCodes.toString());
			//System.out.println("--ALL CODES--");
			//System.out.println(allCodes.toString());
			////allCodes.retainAll(possible);
			Map<String, Integer> editDistMap = getEditDistMap(conflictStr, allCodes);
			//System.out.println(editDistMap.toString());
			allCodes = filterEditDistance(allCodes, editDistMap, VCFConstants.EDIT_DISTANCE);
			//get edit distances between conflict str and all codes
			//get all possible by iterating through all int from 1 to desired dist
			Set<String> conflictMatches = allCodes;
			//System.out.println("--POSSIBLE ALL CODES--");
			//System.out.println(allCodes.toString());
			
			//If conflictMatches is empty, no possible things can be moved,
			//should write the relevant codes to EditSNV, and continue
			if (conflictMatches.isEmpty()){
				ArrayList<VCFEntry> entries = vcfDB.getEntriesByGATK(conflictStr);
				for (VCFEntry entry : entries){
					//System.out.println("**"+entry.toString());
					pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + conflictStr + "\t" + conflictStr + "T\n");
			//		System.out.println("Counter: " + counter++);
				}
			}
			
			ArrayList<String> conflictMatchesList = new ArrayList<String>(conflictMatches);
			//int mutConverted = 0;
			Map<String, ArrayList<VCFEntry>> conflictToPossMutMap = new HashMap<String, ArrayList<VCFEntry>>();
			for (String conflictMatch: conflictMatchesList){
				ArrayList<VCFEntry> possMutations = new ArrayList<VCFEntry>();
				ArrayList<VCFEntry> currFailCodes = new ArrayList<VCFEntry>();
				//boolean is0to1 = checkIf0to1(conflictStr, conflictMatch);
				int currMutConverted =
					vcfDB.getValidEntries(conflictStr, conflictMatch, VCFConstants.EDIT_DISTANCE, possMutations, currFailCodes, null);
				//System.out.println("For code " + conflictStr + " and dest " + conflictMatch + " we can convert " +
					//	currMutConverted + " and cannot convert " + currFailCodes.size() + ".");
				conflictToPossMutMap.put(conflictMatch, possMutations);
				if (failCodes == null && !currFailCodes.isEmpty()) failCodes = new ArrayList<VCFEntry>(currFailCodes);
				else if (!currFailCodes.isEmpty()) failCodes.retainAll(currFailCodes);
				//mutMap.put(conflictMatch, mutMap.get(conflictMatch) + currMutConverted);
				//mutConverted += currMutConverted;
				//System.out.println("Intermediate Mutation Map");
				//System.out.println(mutMap.toString());
			}
			ArrayList<VCFEntry> movedEntries = new ArrayList<VCFEntry>();
			for (int i = 0; i < conflictMatchesList.size(); i++){
				String conflictMatch = findLargestUncoveredSet(conflictToPossMutMap, movedEntries);
				if (conflictMatch == null) break;
				ArrayList<VCFEntry> matchEntries = conflictToPossMutMap.get(conflictMatch);
				matchEntries.removeAll(movedEntries);
				int numMoved = matchEntries.size();
				for (int j = 0; j < matchEntries.size(); j++){
					VCFEntry entry = matchEntries.get(j);
					pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + conflictStr + "\t" + conflictMatch + "\n");
					entry.updateGATK(conflictMatch);
				}
				mutMap.put(conflictStr, mutMap.get(conflictStr) - numMoved);
				mutMap.put(conflictMatch, mutMap.get(conflictMatch) + numMoved);
				//System.out.println("We converted " + numMoved + " entries with GATK code [" + conflictStr + "] to GATK code [" + conflictMatch + "].");
				conflictToPossMutMap.remove(conflictMatch);
				movedEntries.addAll(matchEntries);
			}
			//System.out.println("Total Converted: " + mutConverted);
			/*
			 * 1. Find all possible 1-edit distance changes
			 * 2. See which ones are valid; discard rest
			 * 3. 
			 */
			if (failCodes != null){
				failCodes.removeAll(movedEntries);
				for (int j = 0; j < failCodes.size(); j++){
					VCFEntry entry = failCodes.get(j);
					//pw.write("FAILED\n");
					pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + conflictStr + "\t" + conflictStr + "\n");
				}
			}
		}
		pw.close();
		System.out.println("Final Mutation Map");
		printSNVs(mutMap,vcfDB);

}

	private static Set<String> filterEditDistance(Set<String> allCodes,
			Map<String, Integer> editDistMap, double editDistance) {
		// TODO Auto-generated method stub
		Set<String> resultCodes = new HashSet<String>();
		for (String code : allCodes){
			if (editDistMap.get(code) <= editDistance && editDistMap.get(code) != 0)
				resultCodes.add(code);
		}
		return resultCodes;
	}

	private static Map<String, Integer> getEditDistMap(String conflictStr,
			Set<String> allCodes) {
		// TODO Auto-generated method stub
		Map<String, Integer> edMap = new HashMap<String, Integer>();
		for (String code: allCodes){
			edMap.put(code, getEditDist(conflictStr, code));
		}
		return edMap;
	}

	private static Integer getEditDist(String conflictStr, String code) {
		// TODO Auto-generated method stub
		int dist = 0;
		for (int i = 0; i < conflictStr.length(); i++)
			if (conflictStr.charAt(i) != code.charAt(i)) dist++;
		return dist;
	}

	/**
	 * Function: findLargestUncovered(Map<String, ArrayList<VCFEntry>> conflictToPossMutMap, ArrayList<VCFEntry> movedEntries)
	 * Usage: String largestConflict = findLargestUncoveredSet(conflictToPossMutMap, movedEntries)
	 * ----
	 * This algorithm finds the set with the most "uncovered" objects (or entries) in this case.
	 * This is the brunt of the set cover algorithm. 
	 * 
	 * @param conflictToPossMutMap	The map of a conflict to all possible codes it could be moved to
	 * @param movedEntries			All entries that have been moved already
	 * @return	The binary code where the maximal number of entries can be moved
	 */
	private static String findLargestUncoveredSet(
			Map<String, ArrayList<VCFEntry>> conflictToPossMutMap,
			ArrayList<VCFEntry> movedEntries) {
		String bestMatch = null;
		int maxEntries = 0;
		ArrayList<String> conflictMatchesList = new ArrayList<String>(conflictToPossMutMap.keySet());
		for (String conflictMatch: conflictMatchesList){
			ArrayList<VCFEntry> matchEntries = conflictToPossMutMap.get(conflictMatch);
			matchEntries.removeAll(movedEntries);
			if (matchEntries.size() > maxEntries) {
				bestMatch = conflictMatch;
				maxEntries = matchEntries.size();
			}
		}
		return bestMatch;
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
			HashMap<Integer, Integer> lColFuncMap, VCFDatabase db) {
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
		VisualizeTree(g, null,edgeLabels,db);
	}

	/**
	 * Creates the actual graphics and sets up visualizing the tree
	 * 
	 * Uses a combination of JUNG and Java's Swing libraries to display
	 * tree. This might be made into it's own class if necessary.
	 * 
	 * @param g	A specific instance of a DirectedGraph built already
	 * @param hashMap 
	 */
	private void VisualizeTree(DirectedGraph<Integer, Integer> g, HashMap<Integer, String> nodeLabels, HashMap<Integer, String> edgeLabels, VCFDatabase db) {
		final HashMap<Integer, String> nodeLabelsFinal;
		final HashMap<Integer, String> edgeLabelsFinal;
		
		if (nodeLabels == null) nodeLabelsFinal = new HashMap<Integer, String>();
		else nodeLabelsFinal = new HashMap<Integer, String>(nodeLabels);
		nodeLabelsFinal.put(0, "germline");
		for(int i=0; i < db.getNumofSamples(); i++){
			nodeLabelsFinal.put(-i-1, db.getName(i));
		}
		
		
		if (edgeLabels == null) edgeLabelsFinal = new HashMap<Integer, String>();
		else edgeLabelsFinal = new HashMap<Integer, String>(edgeLabels);
		
		JFrame frame = new JFrame(testName+" Tree View");
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree,100,70);
		//treeLayout.setSize(new Dimension(600, 600));
		radialLayout = new RadialTreeLayout<Integer, Integer>(tree);
		radialLayout.setSize(new Dimension(700, 700));
//		BasicVisualizationServer<Integer, Integer> visServer =
//			new BasicVisualizationServer<Integer, Integer>(treeLayout);
		visServer = 
			new VisualizationViewer<Integer, Integer>(treeLayout);
		visServer.setPreferredSize(new Dimension(600, 500));
		rings = new Rings();
		
		//Setting up transformers for JUNG
		
		Transformer<Integer, String> PhyVertexLabelTransformer = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) return nodeLabelsFinal.get(num);
				/*if (num < 0)
					return db.getName(-1*num - 1);
				else if (num == 0)
					return "Root";*/
				//else return "N-" + num.toString();
				else return null;//"N-" + num.toString();
			}
		};
		
		Transformer<Integer, String> PhyEdgeLabelTransformer = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (edgeLabelsFinal != null && edgeLabelsFinal.containsKey(num)) return edgeLabelsFinal.get(num); 
				else
					return null;
			}
		};
		
//		Transformer<Integer, EdgeShape<Integer, Integer>> PhyEdgeShapeTransformer = new Transformer<Integer, EdgeShape<Integer, Integer>>(){
//			public Line<Integer, Integer> transformer(Integer num){
//				if (num >= 0) return (new EdgeShape.Line<Integer, Integer>());
//				else return (new EdgeShape.BentLine<Integer,Integer>());
//			}
//
//			@Override
//			public EdgeShape<Integer, Integer> transform(Integer num) {
//				if (num >= 0) return (new EdgeShape.Line<Integer, Integer>());
//				else return (new EdgeShape.BentLine<Integer,Integer>());
//			}
//		}
		
		float dash[] = {5.0f};
		final Stroke edgeStroke = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,
				BasicStroke.JOIN_MITER, 10.0f, dash, 0.0f);
		Transformer<Integer, Stroke> PhyEdgeStrokeTransformer =
			new Transformer<Integer, Stroke>() {
				public Stroke transform(Integer num) {
					if (num < 0) return edgeStroke;
					else return null;
				}
			};
			
		Transformer<Integer, Paint> PhyVertexPaintTransformer =
			new Transformer<Integer, Paint>() {
				public Paint transform(Integer num){
					if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) return Color.BLUE;
					if (num < 0) return Color.GREEN;
					else return Color.RED;
				}
			};
		
		//visServer.getRenderContext().setVertexLabelTransformer(new ToStringLabeller<Integer>());
		//visServer.setBackground(Color.WHITE);
		visServer.getRenderContext().setVertexLabelTransformer(PhyVertexLabelTransformer);
		visServer.getRenderContext().setVertexFillPaintTransformer(PhyVertexPaintTransformer);
		visServer.getRenderContext().setEdgeShapeTransformer(new EdgeShape.Line<Integer, Integer>());
		visServer.getRenderContext().setEdgeLabelTransformer(PhyEdgeLabelTransformer);
		visServer.getRenderContext().setEdgeStrokeTransformer(PhyEdgeStrokeTransformer);
		
		//Creating graph mouse
		DefaultModalGraphMouse graphMouse = new DefaultModalGraphMouse();
		graphMouse.setMode(ModalGraphMouse.Mode.TRANSFORMING);
		visServer.setGraphMouse(graphMouse);
		
		Container content = frame.getContentPane();
		//final GraphZoomScrollPane panel = new GraphZoomScrollPane(visServer);
		content.add(visServer);
		
		JToggleButton radial = new JToggleButton("Radial");
		radial.addItemListener(new ItemListener() {

			public void itemStateChanged(ItemEvent e) {
				if(e.getStateChange() == ItemEvent.SELECTED){
					LayoutTransition<Integer, Integer> lt = 
						new LayoutTransition<Integer, Integer>(visServer, treeLayout, radialLayout);
					Animator a = new Animator(lt);
					a.start();
					visServer.getRenderContext().getMultiLayerTransformer().setToIdentity();
					visServer.addPreRenderPaintable(rings);
				} else {
					LayoutTransition<Integer, Integer> lt =
						new LayoutTransition<Integer, Integer>(visServer, radialLayout, treeLayout);
					Animator a = new Animator(lt);
					a.start();
					visServer.getRenderContext().getMultiLayerTransformer().setToIdentity();
					visServer.removePreRenderPaintable(rings);
				}
				visServer.repaint();
			}
		});
		
		JPanel controls = new JPanel();
		//controls.setBackground(Color.WHITE);
		controls.add(radial);
		content.add(controls, BorderLayout.SOUTH);
		
		//JFrame frame = new JFrame("Simple Tree View");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		frame.getContentPane().add(content);
		frame.pack();
		frame.setVisible(true);
		
		Dimension size = frame.getSize();
	      //BufferedImage image = new BufferedImage(size.width, size.height, BufferedImage.TYPE_INT_RGB);
	      BufferedImage image = (BufferedImage)frame.createImage(size.width, size.height);
	      Graphics gr = image.getGraphics();
	      frame.paint(gr);
	      gr.dispose();
	      
	      try
	      {
	        ImageIO.write(image, "jpg", new File(path+testName+".tree.jpg"));
	      }
	      catch (IOException e)
	      {
	        e.printStackTrace();
	      }
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
	
	/**
	 * Function: Rings
	 * Constructor: Rings r = new Rings()
	 * ----
	 * A simple set of shapes taken from the JUNG library
	 * to draw the graph in a radial views.
	 *
	 */
	class Rings implements VisualizationServer.Paintable {
		
		Collection<Double> depths;
		
		public Rings() {
			depths = getDepths();
		}
		
		private Collection<Double> getDepths() {
			Set<Double> depths = new HashSet<Double>();
			Map<Integer,PolarPoint> polarLocations = radialLayout.getPolarLocations();
			for(Integer v : graph.getVertices()) {
				PolarPoint pp = polarLocations.get(v);
				depths.add(pp.getRadius());
			}
			return depths;
		}

		public void paint(Graphics g) {
			g.setColor(Color.lightGray);
		
			Graphics2D g2d = (Graphics2D)g;
			Point2D center = radialLayout.getCenter();

			Ellipse2D ellipse = new Ellipse2D.Double();
			for(double d : depths) {
				ellipse.setFrameFromDiagonal(center.getX()-d, center.getY()-d, 
						center.getX()+d, center.getY()+d);
				Shape shape = visServer.getRenderContext().getMultiLayerTransformer().getTransformer(Layer.LAYOUT).transform(ellipse);
				g2d.draw(shape);
			}
		}

		public boolean useTransform() {
			return true;
		}
	}
}