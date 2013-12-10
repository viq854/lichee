package ppt;

import util.*;


import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import edu.uci.ics.jung.graph.*;
import edu.uci.ics.jung.graph.util.*;

public class SPBuilder {
	
	private static final double SUBPOP_PVALUE = Configs.SUBPOP_PVALUE;
	private static final int EDIT_DISTANCE = Configs.EDIT_DISTANCE;
	private DirectedGraph<Integer, Integer> currTree;
	private HashMap<Integer, String> nodeLabels;
	private HashMap<Integer, String> edgeLabels;

	private SNVDatabase vcfDB;

	public SPBuilder( ArrayList<ArrayList<Integer>> noConflictMatrixPrime, SNVDatabase db){
		PerfectPhylogenyEngine tb = new PerfectPhylogenyEngine();
		vcfDB = db;
		HashMap<Integer, Integer> LColFuncMap = tb.getLColFuncMap(noConflictMatrixPrime);
		Map<Integer, ArrayList<Integer>> edgeToCodeMap = new HashMap<Integer, ArrayList<Integer>>();
		nodeLabels = new HashMap<Integer, String>();
		edgeLabels = new HashMap<Integer, String>();
		currTree = assembleGraph(noConflictMatrixPrime, LColFuncMap, db.getTAG2SNVNum(), edgeToCodeMap);
		
		
		//return currTree;
	}
	
	public HashMap<Integer, String> getNodeLabels(){
		return nodeLabels;
	}
	
	public HashMap<Integer, String> getEdgeLabels(){
		return edgeLabels;
	}
	
	public DirectedGraph<Integer, Integer> getTree(){
		return currTree;
	}
	
	private DirectedGraph<Integer, Integer> assembleGraph(ArrayList<ArrayList<Integer>> matrixPrime, 
			HashMap<Integer, Integer> lColFuncMap, Map<String, Integer> mutMap, Map<Integer, ArrayList<Integer>> edgeToCodeMap) {
		DirectedGraph<Integer, Integer> g = new DirectedSparseGraph<Integer, Integer>();
		System.out.println("SUBPOP MatPrime");
		TreeChecker.printMatrix(matrixPrime);
		System.out.println("-----");
		int numRows = matrixPrime.size();
		int numCols = matrixPrime.get(0).size();
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
			//System.out.println(i + "---");
			//getMatrixColumn(matrixPrime, i);
			edgeToCodeMap.put(i + 1, getMatrixColumn(matrixPrime, i));
			
		}
		System.out.println(edgeToCodeMap.toString());
		//At this point successfully have list of edges in tree and what codes they line up to.
		//Each column of matrix corresponds to edge/intermediate. So column 1 == node1 and edge1
		//Each row corresponds to sample (green nodes). So row 1 = sample A.
		//Next step is to get all legitimate red edges. Only will be between red nodes.
		Map<String, Integer> redEdgeMap = new HashMap<String, Integer>(mutMap);
		System.out.println(redEdgeMap.toString());
		int id = 1;
		for(ArrayList<Integer> greenEdge : edgeToCodeMap.values()){
			String greenEdgeStr = "";
			for (int i = 0; i < greenEdge.size(); i++){
				greenEdgeStr += greenEdge.get(i);
			}
			//System.out.println("Removing Code " + greenEdgeStr + " From Red Edge Map");
			if (mutMap.get(greenEdgeStr) > 0)
				edgeLabels.put(new Integer(id), mutMap.get(greenEdgeStr).toString());
			id++;
			redEdgeMap.remove(greenEdgeStr);
		}
		//Remove germlineStr if in red edge map
		String germlineStr = "";
		for (int i = 0; i < numRows; i++) germlineStr += "1";
		if (redEdgeMap.containsKey(germlineStr)) redEdgeMap.remove(germlineStr);
		
		System.out.println("Red Edges Before Checking Size Limits: " + redEdgeMap.toString());
		removeBadRedEdges(redEdgeMap, numCols, totalMutations(mutMap));
		//System.out.println("Final Red Edges: " + redEdgeMap.toString());
		//Go through each VCF until find one with legit red edge.
		//Distribute the weight of each red edge to connecting nodes.

		//If we find no red edges, add the samples and return the graph
		if (redEdgeMap.isEmpty()){
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
			return g;
		}
		//Redo distribution to distribute redEdge count to each sample
		Map<Integer, Integer> subPopCounts = distributeRedEdges(redEdgeMap, edgeToCodeMap, matrixPrime, g);
		System.out.println("subPopCounts:" + subPopCounts.toString());
		//Integer nodeA = findA(subPopCounts, g);
		Integer nodeA = findABySampleCounts(redEdgeMap, numRows, matrixPrime, g, null);
		System.out.println("A: " + nodeA.toString());
		Set<Integer> ancestors = new HashSet<Integer>();
		Integer nodeB = findParents(nodeA, redEdgeMap, subPopCounts, ancestors, g);
		if (nodeB != null) System.out.println("Node B: " + nodeB.toString());
		if (ancestors != null) System.out.println("Ancestors: " + ancestors.toString());
		Integer nodeC = descendAncestor(nodeA, nodeB, ancestors, redEdgeMap.keySet(), matrixPrime, subPopCounts, g);
		Set<Integer> oldNodeA = new HashSet<Integer>();
		System.out.println("A: " + nodeA.toString() + " B: " + nodeB.toString() + " C: " + nodeC.toString());
		while (nodeC.intValue() == -1){
			oldNodeA.add(nodeA);
			nodeA = findABySampleCounts(redEdgeMap, numRows, matrixPrime, g, oldNodeA);
			System.out.println("A: " + nodeA.toString());
			ancestors = new HashSet<Integer>();
			nodeB = findParents(nodeA, redEdgeMap, subPopCounts, ancestors, g);
			if (nodeB != null) System.out.println("Node B: " + nodeB.toString());
			if (ancestors != null) System.out.println("Ancestors: " + ancestors.toString());
			nodeC = descendAncestor(nodeA, nodeB, ancestors, redEdgeMap.keySet(), matrixPrime, subPopCounts, g);
			//find new NodeA,
			//find new ancestors
			//...
			//find new nodeC
		}
		System.out.println("A: " + nodeA.toString() + " B: " + nodeB.toString() + " C: " + nodeC.toString());
		outputSubPopFile(nodeA, nodeC, matrixPrime, mutMap, numRows);
		outputConflictEdgeInfo(nodeA, nodeC, matrixPrime, mutMap);
		//Need to check whether A or C is leaf of B
		boolean nodeAIsLeaf = g.getSuccessorCount(nodeA) == 0;
		boolean nodeCIsLeaf = g.getSuccessorCount(nodeC) == 0;
		//boolean nodeAParentIsB = g.isPredecessor(nodeB, nodeA);
		boolean nodeAParentIsB = (new ArrayList<Integer>(g.getPredecessors(nodeA))).get(0).intValue() == nodeB.intValue();
		
		//System.out.println(g.getPredecessors(nodeC));
		boolean nodeCParentIsB = (new ArrayList<Integer>(g.getPredecessors(nodeC))).get(0).intValue() == nodeB.intValue();
		//Attach samples
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
		//boolean isLeaf = g.getSuccessorCount(nodeC) == 0;
		//Node Character
		char nodeAChar, nodeCChar;
		Integer subPop = (new ArrayList<Integer>(g.getSuccessors(nodeA))).get(0);
		if (subPop < 0){
			nodeAChar = (char) ('A' + (-1 * subPop.intValue() - 1));
		} else {
			nodeAChar = ("" + subPop.toString()).charAt(0);
		}
		//nodeAChar = (char) ('A' + (-1 * subPop.intValue() - 1));
		System.out.println("Subpop 1 Char: " + nodeAChar);
		subPop = (new ArrayList<Integer>(g.getSuccessors(nodeC))).get(0);
		//System.out.println("Subpop 2: " + subPop.toString());
		if (subPop < 0){
			nodeCChar = (char) ('A' + (-1 * subPop.intValue() - 1));
		} else {
			nodeCChar = ("" + nodeC.toString()).charAt(0);
		}
		System.out.println("Subpop 2 Char: " + nodeCChar);
		

		//Figure out which row corresponds with conflict and find common mutations
		//int row = (int)(nodeChar - 'A') + 1;
		
		//2 Cases:
		//1) A is leaf of B or C is leaf of B
//		if ((nodeAIsLeaf && nodeAParentIsB) || (nodeCIsLeaf && nodeCParentIsB)){
//			if (nodeAIsLeaf && nodeAParentIsB){
		if (nodeAParentIsB || nodeCParentIsB){
			if (nodeAParentIsB){

				Integer newNodeA = ExtendNode(nodeA, nodeC, g, matrixPrime, mutMap, redEdgeMap);
				if (newNodeA != null) {
					nodeA = newNodeA;
				}
				nodeB = (new ArrayList<Integer>(g.getPredecessors(newNodeA))).get(0);
				
				Integer oldEdgeBA = (new ArrayList<Integer>(g.getInEdges(nodeA))).get(0);
				g.removeEdge(oldEdgeBA);
				Integer newEdgeBI = numCols + 1;
				g.addVertex(newEdgeBI);
				g.addEdge(newEdgeBI, nodeB, newEdgeBI, EdgeType.DIRECTED);
				//add back old edge with A
				g.addEdge(oldEdgeBA, newEdgeBI, nodeA, EdgeType.DIRECTED);
				//add new edge with subpop C'
				Integer newEdgeIC = newEdgeBI + 1;
				//Need to add entire subtrees/not single nodes
				g.addVertex(newEdgeIC);
				g.addEdge(newEdgeIC, newEdgeBI, newEdgeIC, EdgeType.DIRECTED);
				//if (nodeCIsLeaf){
					g.addVertex(-1 * newEdgeIC);
					nodeLabels.put(-1 * newEdgeIC, nodeCChar + "-Subpop");
					g.addEdge(-1 * newEdgeIC, newEdgeIC, -1 * newEdgeIC, EdgeType.DIRECTED);
				//} else {
					//recAddSubtree(g, nodeC, newEdgeIC);
					//To Do
				//}
			} else {
				Integer newNodeC = ExtendNode(nodeC, nodeA, g, matrixPrime, mutMap, redEdgeMap);
				if (newNodeC != null) {
					nodeC = newNodeC;
				}
				nodeB = (new ArrayList<Integer>(g.getPredecessors(newNodeC))).get(0);
				
				Integer oldEdgeBC = (new ArrayList<Integer>(g.getInEdges(nodeC))).get(0);
				g.removeEdge(oldEdgeBC);
				Integer newEdgeBI = numCols + 1;
				g.addVertex(newEdgeBI);
				g.addEdge(newEdgeBI, nodeB, newEdgeBI, EdgeType.DIRECTED);
				//add back old edge with C
				g.addEdge(oldEdgeBC, newEdgeBI, nodeC, EdgeType.DIRECTED);
				//add new edge with subpop A'
				Integer newEdgeIA = newEdgeBI + 1;
				g.addVertex(newEdgeIA);
				g.addEdge(newEdgeIA, newEdgeBI, newEdgeIA, EdgeType.DIRECTED);
				g.addVertex(-1 * newEdgeIA);
				nodeLabels.put(-1 * newEdgeIA, nodeAChar + "-Subpop");
				g.addEdge(-1 * newEdgeIA, newEdgeIA, -1 * newEdgeIA, EdgeType.DIRECTED);
			}
		} else {
			Integer newEdgeBI = numCols + 1;
			g.addVertex(newEdgeBI);
			g.addEdge(newEdgeBI, nodeB, newEdgeBI, EdgeType.DIRECTED);
			//add A'
			Integer newEdgeIA = newEdgeBI + 1;
			g.addVertex(newEdgeIA);
			g.addEdge(newEdgeIA, newEdgeBI, newEdgeIA, EdgeType.DIRECTED);
			g.addVertex(newEdgeIA * -1);
			nodeLabels.put(-1* newEdgeIA, nodeAChar + "-Subpop");
			g.addEdge(newEdgeIA * -1, newEdgeIA, newEdgeIA * -1, EdgeType.DIRECTED);
			//add C'
			Integer newEdgeIC = newEdgeIA + 1;
			g.addVertex(newEdgeIC);
			g.addEdge(newEdgeIC, newEdgeBI, newEdgeIC, EdgeType.DIRECTED);
			g.addVertex(newEdgeIC * -1);
			nodeLabels.put(-1* newEdgeIC, nodeCChar + "-Subpop");
			g.addEdge(newEdgeIC * -1, newEdgeIC, newEdgeIC * -1, EdgeType.DIRECTED);
		}
		//2) Otherwise
		
		//3 Cases:
//		if (nodeC.equals(0) || ancestors.contains(nodeC)){
//			//1) B is ancestor (and A is not leaf to common ancestor), no conflicts in subtree of B
//			//Add A, B, and intermediate
//			Integer newEdgeI, newEdgeA, newEdgeB;
//			newEdgeI = numCols + 1;
//			newEdgeA = newEdgeI + 1;
//			newEdgeB = newEdgeA + 1;
//			//Attach newEdgeI and intermediate vertex
//			g.addVertex(newEdgeI);
//			g.addEdge(newEdgeI, nodeC, newEdgeI, EdgeType.DIRECTED);
//			//Attach newEdgeA
//			g.addVertex(newEdgeA);
//			g.addEdge(newEdgeA, newEdgeI, newEdgeA, EdgeType.DIRECTED);
//			//Attach newEdgeB
//			g.addVertex(newEdgeB);
//			g.addEdge(newEdgeB, newEdgeI, newEdgeB, EdgeType.DIRECTED);
//			//NEED TO ADD IN WHERE SAMPLES GO
//			//nodeLabels.put(-1 * newEdgeA, "" + nodeChar + "-Subpop")
//			
//		} else if (isLeaf){
//			//2) C is leaf
//			Integer edge2 = (new ArrayList<Integer>(g.getInEdges(nodeC))).get(0);
//			Integer parent = g.getSource(edge2);
//			//remove old edge
//			g.removeEdge(edge2);
//			Integer newEdge1 = numCols + 1;
//			//add intermediate
//			g.addVertex(newEdge1);
//			g.addEdge(newEdge1, parent, newEdge1, EdgeType.DIRECTED);
//			Integer newEdge2 = newEdge1 + 1;
//			//add A' 
//			g.addVertex(newEdge2);
//			g.addEdge(newEdge2, newEdge1, newEdge2, EdgeType.DIRECTED);
//			g.addVertex(-1 * newEdge2);
//			//nodeLabels.put(-1 * newEdge2, "" + A.toString() + "S");
//			System.out.println(g.getSuccessors(nodeA).toString());
//			nodeLabels.put(-1 * newEdge2, "" + nodeChar + "-Subpop");
//			g.addEdge(-1 * newEdge2, newEdge2, -1 * newEdge2, EdgeType.DIRECTED);
//			//add back B
//			g.addEdge(edge2, newEdge1, nodeC, EdgeType.DIRECTED);
//			
//		} else {
//			//3) B is intermediate between 2 & 3
//			Integer newEdge;
//			newEdge = numCols + 1;
//			g.addVertex(newEdge);
//			g.addEdge(newEdge, nodeC, newEdge, EdgeType.DIRECTED);
//			
//			g.addVertex(-1 * newEdge);
//			g.addEdge(-1 * newEdge, newEdge, -1 * newEdge, EdgeType.DIRECTED);
//			nodeLabels.put(-1 * newEdge, "" + nodeChar + "-Subpop");
//		}
		//pw.close()
		return g;
		//graph = g;
		//VisualizeTree(g);
	}

	private Integer ExtendNode(Integer currNode,
			Integer nodeA, DirectedGraph<Integer, Integer> g, ArrayList<ArrayList<Integer>> matrixPrime,
			Map<String, Integer> mutMap, Map<String, Integer> redEdgeMap) {
		// TODO Auto-generated method stub
		//1) Get all successors of currNode
		//2) Get highest successor conversion prob
		//3) If over 0.5, shift A to successor
		//Get all successors of 
		ArrayList<Integer> currNodeCode = getMatrixColumn(matrixPrime, currNode - 1);
		String currNodeStr = "";
		for (int j = 0; j < currNodeCode.size(); j++) currNodeStr += currNodeCode.get(j);
		ArrayList<Integer> nodeACode = getMatrixColumn(matrixPrime, nodeA - 1);
		String nodeAStr = "";
		for (int j = 0; j < nodeACode.size(); j++) nodeAStr += nodeACode.get(j);
		int valueCurrNode = Integer.valueOf(currNodeStr, 2).intValue();
		int valueA = Integer.valueOf(nodeAStr, 2).intValue();
		int conflictCode = valueA | valueCurrNode;
		//System.out.println("ConflictCode: " + conflictCode);
		String conflictCodeStr = Integer.toBinaryString(conflictCode);
		while (conflictCodeStr.length() < matrixPrime.size()){
			conflictCodeStr = "0" + conflictCodeStr;
		}
		System.out.println("ConflictCode: " + conflictCodeStr);
		System.out.println("Testing extending currNode: " + currNodeCode.toString() + " to: " + conflictCodeStr.toString());
		ArrayList<SNVEntry> possMutations = new ArrayList<SNVEntry>();
		ArrayList<SNVEntry> currFailCodes = new ArrayList<SNVEntry>();
		Integer totalMut = vcfDB.getEntriesByGATK(currNodeStr).size();
		int currMutConverted = vcfDB.getValidEntries(currNodeStr, conflictCodeStr, possMutations, currFailCodes);
		double conversionRate = (currMutConverted + 0.0) / totalMut;
		System.out.println("We were able to convert " + currMutConverted + " for subpop C-Extension for a rate of " + 
				(conversionRate));
		boolean canMove = false;
		if (conversionRate > 0.5) canMove = true;
		PrintWriter pw = null;
		if (canMove){
			try {
				pw = new PrintWriter(new FileWriter(Configs.path + Configs.testName + "_editSNV.txt", true));
				//pw.write("SUBPOPULATION MOVEMENTS\n-----\n");
				for (int i = 0; i < possMutations.size(); i++){
					SNVEntry entry = possMutations.get(i);
					pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + currNodeStr + "\t" + conflictCodeStr + "\n");
				}
				for (int i = 0; i < currFailCodes.size(); i++){
					SNVEntry entry = currFailCodes.get(i);
					pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + currNodeStr + "\t" + currNodeStr + "\n");
				}
				pw.close();
			} catch(IOException e) {
				System.out.println("COULDN'T REOPEN FILE");
			}
		}
		
		//ArrayList<Integer> codeC = getMatrixColumn(matrixPrime, nodeC - 1);
		Integer newNodeB = null;
		while (canMove){
			int count = g.getSuccessorCount(currNode);
			if (g.getSuccessorCount(currNode) == 1) break;
			ArrayList<Integer> currSuccessors = new ArrayList<Integer>(g.getSuccessors(currNode));
			//double maxRate = 0.0;
			int maxCount = 0;
			Integer oldNode = new Integer(currNode);
			//Need to output for each generation
			//1) Codes which were successfully converted
			//2) Codes which failed
			Map<Integer, ArrayList<SNVEntry>> nodesToPossMutMap = new HashMap<Integer, ArrayList<SNVEntry>>();
			Map<Integer, ArrayList<SNVEntry>> nodesToFailCodes = new HashMap<Integer, ArrayList<SNVEntry>>();
			Map<Integer, String> nodesToSuccStr = new HashMap<Integer, String>();
			Map<Integer, String> nodesToConflictStr = new HashMap<Integer, String>();
			for (int i = 0; i < currSuccessors.size(); i++){
				Integer currSuccessor = currSuccessors.get(i);
				if (!isRedEdgeEndpoints(nodeA, currSuccessor, redEdgeMap.keySet(), matrixPrime)) continue;
				ArrayList<Integer> currSuccessorCode = getMatrixColumn(matrixPrime, currSuccessor - 1);
				System.out.println("Testing extending currNode: " + currNodeCode.toString() + " to: " + currSuccessorCode.toString());
				String succCodeStr = "";
				
				for (int j = 0; j < currSuccessorCode.size(); j++) succCodeStr += currSuccessorCode.get(j);

				//int valueA = Integer.valueOf(currNodeStr, 2).intValue();
				int currSuccValue = Integer.valueOf(succCodeStr, 2).intValue();
				conflictCode = valueA | currSuccValue;
				System.out.println("ConflictCode: " + conflictCode);
				conflictCodeStr = Integer.toBinaryString(conflictCode);
				while (conflictCodeStr.length() < matrixPrime.size()){
					conflictCodeStr = "0" + conflictCodeStr;
				}
				System.out.println("ConflictCode: " + conflictCodeStr);
				possMutations = new ArrayList<SNVEntry>();
				currFailCodes = new ArrayList<SNVEntry>();
				totalMut = vcfDB.getEntriesByGATK(succCodeStr).size();
				currMutConverted = vcfDB.getValidEntries(succCodeStr, conflictCodeStr, possMutations, currFailCodes);
				int totalMutConverted = currMutConverted;
				if (mutMap.containsKey(conflictCodeStr)) {
					totalMutConverted = currMutConverted + mutMap.get(conflictCodeStr).intValue();
					System.out.println("The group originally had " + mutMap.get(conflictCodeStr).intValue() + " and now has " + totalMutConverted + " conversions.");
				} else {
					System.out.println("The group was previously empty.");
				}
			//	conversionRate = (currMutConverted + 0.0) / totalMut;
				System.out.println("We were able to convert " + currMutConverted + "."); 
				nodesToPossMutMap.put(currSuccessor, possMutations);
				nodesToFailCodes.put(currSuccessor, currFailCodes);
				nodesToSuccStr.put(currSuccessor, succCodeStr);
				nodesToConflictStr.put(currSuccessor, conflictCodeStr);
				//"+ for subpop C-Extension for a rate of " + (conversionRate));
				//canMove = false;
				//if (conversionRate > 0.5 && conversionRate > maxRate){ //canMove = true;
				//if(currMutConverted > maxCount){
				if(totalMutConverted > maxCount){
					//maxRate = conversionRate;
					//maxCount = currMutConverted;
					maxCount = totalMutConverted;
					//newNodeB = currNode;
					currNode = currSuccessor;
				}
			}
			if (currNode.intValue() == oldNode.intValue()) break;
			else {
				try {
					pw = new PrintWriter(new FileWriter(Configs.path + Configs.testName + "_editSNV.txt", true));
					//pw.write("-----\n");
					possMutations = nodesToPossMutMap.get(currNode);
					currFailCodes = nodesToFailCodes.get(currNode);
					String succCodeStr = nodesToSuccStr.get(currNode);
					conflictCodeStr = nodesToConflictStr.get(currNode);
					for (int i = 0; i < possMutations.size(); i++){
						SNVEntry entry = possMutations.get(i);
						pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + succCodeStr + "\t" + conflictCodeStr + "\n");
					}
					for (int i = 0; i < currFailCodes.size(); i++){
						SNVEntry entry = currFailCodes.get(i);
						pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + succCodeStr + "\t" + succCodeStr + "\n");
					}
					pw.close();
				} catch(IOException e) {
					System.out.println("COULDN'T REOPEN FILE");
				}
			}
		}
		//nodeA = currNode;
		//currNode = new Integer(currNode);
		//newNodeB = (new ArrayList<Integer>(g.getPredecessors(currNode))).get(0);
		return currNode;
		
	}

	private void outputConflictEdgeInfo(Integer nodeA, Integer nodeC,
			ArrayList<ArrayList<Integer>> matrixPrime,
			Map<String, Integer> mutMap) {
		// TODO Auto-generated method stub
	ArrayList<Integer> codeA = getMatrixColumn(matrixPrime, nodeA - 1);
	ArrayList<Integer> codeC = getMatrixColumn(matrixPrime, nodeC - 1);
	System.out.println("A: " + codeA.toString() + "C: " + codeC.toString());
	String ACodeStr = "", CCodeStr = "";
	for (int i = 0; i < codeA.size(); i++) ACodeStr += codeA.get(i);
	for (int i = 0; i < codeC.size(); i++) CCodeStr += codeC.get(i);
	ArrayList<SNVEntry> possMutations = new ArrayList<SNVEntry>();
	ArrayList<SNVEntry> currFailCodes = new ArrayList<SNVEntry>();
	int valueA = Integer.valueOf(ACodeStr, 2).intValue();
	int valueC = Integer.valueOf(CCodeStr, 2).intValue();
	int conflictCode = valueA | valueC;
	System.out.println("ConflictCode: " + conflictCode);
	String conflictCodeStr = Integer.toBinaryString(conflictCode);
	while (conflictCodeStr.length() < matrixPrime.size()){
		conflictCodeStr = "0" + conflictCodeStr;
	}
	System.out.println("ConflictCode: " + conflictCodeStr);
	Integer totalMut = vcfDB.getEntriesByGATK(CCodeStr).size();
	int currMutConverted = vcfDB.getValidEntries(CCodeStr, conflictCodeStr, possMutations, currFailCodes);
	System.out.println("We were able to convert " + currMutConverted + " for subpop C-Extension for a rate of " + 
			((currMutConverted + 0.0) / totalMut ));
	}

	private Integer findABySampleCounts(Map<String, Integer> redEdgeMap,
			int numRows, ArrayList<ArrayList<Integer>> matrixPrime, DirectedGraph<Integer, Integer> g, Set<Integer> oldNodeASet) {
		// TODO Auto-generated method stub
		Map<Integer, Integer> sampleCounts = new HashMap<Integer, Integer>();
		for (String redEdge : redEdgeMap.keySet()){
			for (int i = 0; i < redEdge.length(); i++){
				if (redEdge.charAt(i) == '1'){
					if (!sampleCounts.containsKey(i)){
						sampleCounts.put(i, redEdgeMap.get(redEdge).intValue());
					} else {
						sampleCounts.put(i, sampleCounts.get(i) + redEdgeMap.get(redEdge).intValue());
					}
				}
			}
		}
		int max = -1;
		Integer bestSample = null;
		Set<Integer> sampleCountsKeySet = new HashSet<Integer>(sampleCounts.keySet());
		//remove previously seen nodeAs
		//if (oldNodeASet != null) sampleCountsKeySet.removeAll(oldNodeASet);
		for (Integer sample : sampleCountsKeySet){
			int count = sampleCounts.get(sample);
			if (max == -1 || count > max){
				max = count;
				bestSample = sample;
			}
		}
		System.out.println("a) Sample Counts: " + sampleCounts.toString());
		Integer nodeA = null;
		for (int i = 0; i < matrixPrime.get(0).size(); i++){
			if (matrixPrime.get(bestSample).get(i).intValue() == 1){
				if (isSingleSampleEdge(bestSample, i, matrixPrime)) {
					nodeA = i + 1;
					break;
				}
			}
		}
		/*
		 * while (oldNodeASet != null && oldNodeASet.contains(nodeA)){
			max = -1;
			sampleCountsKeySet.remove(bestSample);
			for (Integer sample : sampleCountsKeySet){
				int count = sampleCounts.get(sample);
				if (max == -1 || count > max){
					max = count;
					bestSample = sample;
				}
			}
			System.out.println("b) Sample Counts: " + sampleCounts.toString());
			nodeA = null;
			for (int i = 0; i < matrixPrime.get(0).size(); i++){
				if (matrixPrime.get(bestSample).get(i).intValue() == 1){
					if (isSingleSampleEdge(bestSample, i, matrixPrime)) {
						nodeA = i + 1;
						break;
					}
				}
			}
		}
		*/
		//Need to extend A, keep going UP as much as possible;
		//1) currNode = A.parent
		//2) while all of currNodes's children do have samples with value max
		//3) 	currNode = currNode.parent
		Integer currNode = (new ArrayList<Integer>(g.getPredecessors(nodeA))).get(0);
		
		ArrayList<Integer> currColumn = null;
		if (currNode.intValue() != 0) {
			currColumn = getMatrixColumn(matrixPrime, currNode - 1);
			System.out.println("CurrColumn: " + currColumn.toString());
		}
		while (currNode.intValue() != 0){
			for (int i = 0; i < currColumn.size(); i++){
				if (currColumn.get(i).intValue() == 1 && sampleCounts.containsKey(i) && sampleCounts.get(i) != max) return nodeA;
				else if (currColumn.get(i).intValue() == 1 && !sampleCounts.containsKey(i)) return nodeA;
			}
			nodeA = currNode;
			currNode = (new ArrayList<Integer>(g.getPredecessors(nodeA))).get(0);
			if (currNode.intValue() != 0) currColumn = getMatrixColumn(matrixPrime, currNode - 1);
		}
		return nodeA;
	}
//	private Integer recAddSubtree(DirectedGraph<Integer, Integer> g,
//			Integer nodeC, Integer subRoot) {
//		// TODO Auto-generated method stub
//		Queue<Queue<Integer>> successorGroupQueue = new LinkedList<Queue<Integer>>();
//		Queue<Integer> initialQueue = new LinkedList<Integer>()
//		initialQueue.add(nodeC);
//		successorGroupQueue.add(initialQueue);
//		int currNodeLabel = subRoot.intValue();
//		while (!successorGroupQueue.isEmpty()){
//			Queue<Integer> currSuccessorQueue = successorGroupQueue.remove();
//			while (!currSuccessorQueue.isEmpty()){
//				Integer currNode = currSuccessorQueue.remove();
//				ArrayList<Integer> currNodeSuccessors = new ArrayList<Integer>(g.getSuccessors(currNode));
//				Queue<Integer> newSuccessorQueue = new LinkedList<Integer>();
//				for (int i = 0; i < currNodeSuccessors.size(); i++){
//					//Integer currSuccessor = new Integer(currSuccessors.get(i));
//					//create a copy of the successor and attach to newSubRoot
//					Integer currSuccessor = currNodeSuccessors.get(i);
//					if (currSuccessor > 0){
//						Integer successorCopy = new Integer(++currNodeLabel);
//						g.addVertex(successorCopy);
//						g.addEdge(successorCopy, subRoot, successorCopy);
//						newSuccessorQueue.add(currNodeSuccessors.get(i));
//					} else {
//						Integer sampleCopy = new Integer(-1 * currNode);
//						g.addVertex(sampleCopy);
//						g.addEdge(sampleCopy, arg1)
//					}
//				}
//				successorGroupQueue.add(newSuccessorQueue);
//			}
//		}
//		for (int i = 0; i < successorGroupQueue.size(); i++){
//			
//		}
//		return null;
//		
//	}

	private boolean isSingleSampleEdge(Integer bestSample, int col, ArrayList<ArrayList<Integer>> matrixPrime) {
	// TODO Auto-generated method stub
		//boolean isSampleEdge = false;
		ArrayList<Integer> currCol = getMatrixColumn(matrixPrime, col);
		for (int j = 0; j < currCol.size(); j++){
			if (j != bestSample && currCol.get(j).intValue() != 0) return false;
						
		}
		return true;
	}

	private void outputSubPopFile(Integer A, Integer B,
			ArrayList<ArrayList<Integer>> matrixPrime,
			Map<String, Integer> mutMap, int numRows) {
		// TODO Auto-generated method stub
		ArrayList<Integer> ACode = getMatrixColumn(matrixPrime, A - 1);
		ArrayList<Integer> BCode = getMatrixColumn(matrixPrime, B - 1);
		String ACodeStr = "", BCodeStr = "";
		for (int i = 0; i < ACode.size(); i++) ACodeStr += ACode.get(i);
		for (int i = 0; i < BCode.size(); i++) BCodeStr += BCode.get(i);
		int Avalue = Integer.valueOf(ACodeStr, 2).intValue();
		int Bvalue = Integer.valueOf(BCodeStr, 2).intValue();
		int subPopValue = (Avalue | Bvalue);
		String subPopStr = Integer.toString(subPopValue, 2);
		System.out.println("SubPop: " + subPopStr + " NumRows: " + numRows);
		if (subPopStr.length() < numRows){
			while (subPopStr.length() < numRows){
				subPopStr = "0" + subPopStr;
			}
		}
		Integer subPopNum = mutMap.get(subPopStr);
		String outputLine = "SubPop of " + subPopStr + " was found and resolved with count " + subPopNum + ".";
		//System.out.println();
		PrintWriter pw = null;
		try{
			pw = new PrintWriter(new FileWriter(Configs.path + Configs.testName + "_SubPopOutput.txt"));
			pw.write(outputLine);
		} catch (IOException e){
			
		}
		pw.close();
	}

	private Integer descendAncestor(Integer nodeA, Integer nodeB, Set<Integer> ancestors, Set<String> redEdges, 
			ArrayList<ArrayList<Integer>> matrixPrime, Map<Integer, Integer> subPopCounts, DirectedGraph<Integer, Integer> g) {
		// TODO Auto-generated method stub
		ArrayList<Integer> successors = new ArrayList<Integer>(g.getSuccessors(nodeB));
		successors.removeAll(ancestors);
		Queue<Integer> queue = new LinkedList<Integer>();
		for (int i = 0; i < successors.size(); i++){
			queue.add(successors.get(i));
		}
		while (!queue.isEmpty()){
			Integer currNode = queue.remove();
			if (isRedEdgeEndpoints(nodeA, currNode, redEdges, matrixPrime)) return currNode;
			ArrayList<Integer> currSuccessors = new ArrayList<Integer>(g.getSuccessors(currNode));
			for (Integer successor: currSuccessors) queue.add(successor);
		}
		return -1;
//		if (successors.size() > 1){
//			
//			System.out.println("NOT SURE");
//			return null;
//		} else {
//			Integer currNode = successors.get(0);
//			Integer prevNode = null;
//			while (subPopCounts.keySet().contains(currNode)) {
//				if (g.getSuccessorCount(currNode) == 0) return currNode;
//				prevNode = currNode;
//				currNode = (new ArrayList<Integer>(g.getSuccessors(prevNode))).get(0);
//			}
//			return prevNode;
//		}
//		
		
	}

	private boolean isRedEdgeEndpoints(Integer nodeA, Integer currNode,
			Set<String> redEdges, ArrayList<ArrayList<Integer>> matrixPrime) {
		// TODO Auto-generated method stub
		ArrayList<Integer> nodeACode = getMatrixColumn(matrixPrime, nodeA - 1);
		ArrayList<Integer> currNodeCode = getMatrixColumn(matrixPrime, currNode - 1);
		for (String redEdge : redEdges){
			if (isEdgeEndpoints(nodeACode, currNodeCode, redEdge)) return true;
		}
		return false;
		//return isEdgeEndPoints(nodeACode, currNodeCode, )
	}

	private int ALtoInt(ArrayList<Integer> code) {
		// TODO Auto-generated method stub
		String codeStr = "";
		for (int i = 0; i < code.size(); i++) codeStr += code.get(i);
		int value = Integer.valueOf(codeStr, 2).intValue();
		return value;
	}

	private Integer findParents(Integer A,
			Map<String, Integer> redEdgeMap, Map<Integer, Integer> subPopCounts, Set<Integer> ancestors, DirectedGraph<Integer, Integer> g) {
		if (A == null) return null;
		//Set<Integer> parents = new HashSet<Integer>();
		Integer newParent = (new ArrayList<Integer>(g.getPredecessors(A))).get(0);
		//Set<Integer> alreadyCovered = new HashSet<Integer>();
		//alreadyCovered.add(A);
		System.out.println(newParent);
		//find count for first parent
		ArrayList<Integer> currParentsChildren = new ArrayList<Integer>(g.getSuccessors(newParent));
		currParentsChildren.remove(A);
		double score = getScore(currParentsChildren, subPopCounts, g);
		Integer bestAncestor = newParent;
		//if (newParent.equals(0)); //A is leaf of root;
		while (!newParent.equals(0)){
			ancestors.add(newParent);
			newParent = (new ArrayList<Integer>(g.getPredecessors(newParent))).get(0);
			ArrayList<Integer> newSuccessors = new ArrayList<Integer>(g.getSuccessors(newParent));
			newSuccessors.removeAll(ancestors);
			double newScore = getScore(newSuccessors, subPopCounts, g);
			if (newScore > score){
				score = newScore;
				bestAncestor = newParent;
			}
			
		}
		//ancestors.add(newParent);
		//ancestors = parents;
		System.out.println(score);
		//newParent = (new ArrayList<Integer>(g.getPredecessors(newParent))).get(0);
		return bestAncestor;
	}

	private double getScore(ArrayList<Integer> scoringRoots, Map<Integer, Integer> subPopCounts, DirectedGraph<Integer, Integer> g) {
		Set<Integer> allScoringNodes = new HashSet<Integer>();
		allScoringNodes.addAll(scoringRoots);
		//while 
		for (int i = 0; i < scoringRoots.size(); i++){
			Integer currRoot = scoringRoots.get(i);
			Stack<Integer> toCheck = new Stack<Integer>();
			toCheck.addAll(g.getSuccessors(currRoot));
//			if (g.getSuccessorCount(currRoot) != 0){
//				Set<Integer> successors = new HashSet<Integer>(g.getSuccessors(currRoot));
//				allScoringNodes.addAll(successors);
//				
//			}
			while (!toCheck.isEmpty()){
				Integer currNode = toCheck.pop();
				allScoringNodes.add(currNode);
				toCheck.addAll(g.getSuccessors(currNode));
			}
			//System.out.println("Curr: " + currRoot.toString());
			//System.out.println("Is leaf: " + (g.getSuccessorCount(currRoot) == 0));
		}
		//Take out any nodes that do not have red edges
		allScoringNodes.retainAll(subPopCounts.keySet());
		double total = 0.0;
		for (Integer scoringNode : allScoringNodes){
			total += subPopCounts.get(scoringNode);
		}
		System.out.println(allScoringNodes.toString());
		return total;
	}

	private Integer findA(Map<Integer, Integer> subPopCounts,
			DirectedGraph<Integer, Integer> g) {
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		//find highest count
		double max = 0.0;
		Integer A = null;
		for (Map.Entry<Integer, Integer> entry: subPopCounts.entrySet()){
			if (entry.getValue() > max) {
				max = entry.getValue();
				A = entry.getKey();
			} else if (entry.getValue().equals(max)){
				A = getDeeper(entry.getKey(), A, tree);
			}
		}
		//System.out.println(tree.getDepth((Integer) 3));
		//System.out.println(tree.getDepth(4));
		// TODO Auto-generated method stub
		return A;
	}

	private Integer getDeeper(Integer key, Integer A,
			DelegateTree<Integer, Integer> tree) {
		Integer root = tree.getRoot();
		Set<Integer> children = new HashSet<Integer>(tree.getChildren(root));
		while (true){
			if (children.contains(key)) return A;
			else if (children.contains(A)) return key;
			else {
				Set<Integer> newChildren = new HashSet<Integer>();
				for (Integer child: children){
					newChildren.addAll(tree.getChildren(child));
				}
				children = newChildren;
			}
		}
		// TODO Auto-generated method stub
		//return null;
	}

	
	//Fixed?
	private Map<Integer, Integer> distributeRedEdges(Map<String, Integer> redEdgeMap,
			Map<Integer, ArrayList<Integer>> edgeToCodeMap,
			ArrayList<ArrayList<Integer>> matrixPrime, DirectedGraph<Integer, Integer> g) {
		//Get all red edge codes
		ArrayList<String> redEdgeCodes = new ArrayList<String>(redEdgeMap.keySet());
		System.out.println("Red Edge Codes: " + redEdgeMap.toString());
		//Set up map of nodes and counts
		Map<Integer, Integer> subPopCounts = new HashMap<Integer, Integer>();
		//Distribute each edge
		for (int i = 0; i < redEdgeCodes.size(); i++){
			String redEdge = redEdgeCodes.get(i);
			for (int x = 0; x < matrixPrime.get(0).size(); x++){
				for (int y = 0; y < x; y++){
					if (x != y && isEdgeEndpoints(getMatrixColumn(matrixPrime, x), getMatrixColumn(matrixPrime, y), redEdge)){
						Integer endPoint1 = x + 1;
						Integer endPoint2 = y + 1;
						if (subPopCounts.containsKey(endPoint1)){
							subPopCounts.put(endPoint1, subPopCounts.get(endPoint1) + redEdgeMap.get(redEdge));
						} else {
							subPopCounts.put(endPoint1, redEdgeMap.get(redEdge));
						}
						if (subPopCounts.containsKey(endPoint2)){
							subPopCounts.put(endPoint2, subPopCounts.get(endPoint2) + redEdgeMap.get(redEdge));
						} else {
							subPopCounts.put(endPoint2, redEdgeMap.get(redEdge));
						}
					}
				}
			}
			//Try every pair of edges to find conflict
			
//			Integer leftEdgeNode = 0, rightEdgeNode = 0;
//			ArrayList<Integer> currChildren = new ArrayList<Integer>(g.getSuccessors(0));
//			System.out.println(currChildren.toString());
			//Right now works only if root has two children
//			if(currChildren.size() != 2){
//				System.out.println("UNDEFINED BEHAVIOR");
//			} else {
//				leftEdgeNode = currChildren.get(0);
//				rightEdgeNode = currChildren.get(1);
//			}
//			System.out.println(leftEdgeNode);
//			System.out.println(rightEdgeNode);
			
			//traverse left and right
		}
		//Every edge is now distributed among nodes
		//System.out.println(subPopCounts.toString());
		return subPopCounts;
	}

	private boolean isEdgeEndpoints(ArrayList<Integer> code1, ArrayList<Integer> code2, String redEdgeCode) {
		// TODO Auto-generated method stub
		String code1Str = "", code2Str = "";
		for (int i = 0; i < code1.size(); i++) code1Str += code1.get(i);
		for (int i = 0; i < code2.size(); i++) code2Str += code2.get(i);
		int code1Int = Integer.valueOf(code1Str, 2).intValue();
		int code2Int = Integer.valueOf(code2Str, 2).intValue();
		int destCode = Integer.valueOf(redEdgeCode, 2).intValue();
		//if one is subset of other, return false
		if ((code1Int | code2Int) == code1Int || (code1Int | code2Int) == code2Int) return false;
		if ((code1Int | code2Int) == destCode) return true;
		return false;
	}

	private void removeBadRedEdges(Map<String, Integer> redEdgeMap, int numEdges,
			int totalMutations) {
		ArrayList<String> badEdges = new ArrayList<String>();
		System.out.println("n (Total Mutations) = " + totalMutations);
		System.out.println("k (Number of Edges) = " + numEdges);
		int threshold = getConflictThreshold(totalMutations, numEdges);
		System.out.println("Threshold: " + threshold);
		for (String key : redEdgeMap.keySet()){
			if (redEdgeMap.get(key) < threshold){
				badEdges.add(key);
			}
		}
		for (String badEdge: badEdges){
			redEdgeMap.remove(badEdge);
		}
	}
	
	private static int getConflictThreshold(int n, int k){
		//k++;
		double pValue = 0.0;
		int x = n / k;
		for (; x > 1; x--){
			pValue = Math.pow((n - x * k)/(double) n, k - 1);
			if (pValue > SUBPOP_PVALUE) break;
		}
		return x;
	}
	

	private int totalMutations(Map<String, Integer> mutMap) {
		// TODO Auto-generated method stub
		ArrayList<Integer> values = new ArrayList<Integer>(mutMap.values());
		int total = 0;
		for (int i = 0; i < values.size(); i++){
			total += values.get(i);
		}
		return total;
	}

	private ArrayList<Integer> getMatrixColumn(
			ArrayList<ArrayList<Integer>> matrixPrime, int col) {
		ArrayList<Integer> colCode = new ArrayList<Integer>();
		for (int i = 0; i < matrixPrime.size(); i++){
			ArrayList<Integer> row = matrixPrime.get(i);
			//for (int j = 0; j < row.size(); j++) System.out.print(row.get(i));
			colCode.add(row.get(col));
			//System.out.println(row.get(col).toString());
		}
		return colCode;
	}
}
