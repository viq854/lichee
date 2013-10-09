package lineage;

import java.util.ArrayList;
import java.util.HashMap;

import lineage.AAFClusterer.Cluster;

/**
 * PHYGraph is a directed constraint graph representing the phylogenetic relationship
 * among sample sub-populations.
 * Each internal node in the graph represents a sub-population. 
 * The leaves of the graph represent input cancer samples.
 * A directed edge between two nodes denotes the 'happened-before' evolutionary 
 * relationship between the two nodes.
 */
public class PHYGraph {

	/** Nodes in the graph divided by levels (number of samples node SNPs occurred in) */
	private HashMap<Integer, ArrayList<PHYNode>> nodes;
	
	/** Adjacency map of nodes to their down-neighbors */
	private HashMap<PHYNode, ArrayList<PHYNode>> edges;
	
	/** Total number of nodes in the graph */
	private int numNodes;
	
	/** Total number of edges in the graph */
	private int numEdges;
	
	/** Total number of cancer samples */
	private int numSamples;
	
	/** Error margin used for comparing AAF centroid data */
	private static final double AAF_ERROR_MARGIN = 0.08;
	
	/** Maximum AAF */
	private static final double AAF_MAX = 0.5;
	
	/** Debugging-only */
	protected static int nodeCounter = 0;
	
	/**
	 * Constructs a PHYGraph from the sub-populations of the SNP groups
	 */
	public PHYGraph(ArrayList<SNPGroup> groups, int totalNumSamples) {
		numSamples = totalNumSamples;
		nodes = new HashMap<Integer, ArrayList<PHYNode>>();
		edges = new HashMap<PHYNode, ArrayList<PHYNode>>(); 
		
		// add root node
		addNode(new PHYNode(), numSamples+1);
				
		// add group sub-population nodes
		for(SNPGroup g : groups) {
			PHYNode[] groupNodes = new PHYNode[g.getSubPopulations().length];
			for(int i = 0; i < groupNodes.length; i++) {
				PHYNode node = new PHYNode(g, i);
				addNode(node, g.getNumSamples());
				groupNodes[i] = node;
			}
			// add edges between each group's sub-population nodes
			for(int i = 0; i < groupNodes.length; i++) {
				for(int j = i+1; j <  groupNodes.length; j++) {
					checkAndAddEdge(groupNodes[i], groupNodes[j]);
				}
			}
		}
		
		// add sample leaf nodes
		for(int i = 0; i < numSamples; i++) {
			addNode(new PHYNode(i), 0);
		}
		
		// add inter-level edges
		for(int i = numSamples + 1; i > 0; i--) {
			for(int j = i-1; j >= 0; j--) {
				ArrayList<PHYNode> fromLevelNodes = nodes.get(i);
				if(fromLevelNodes == null) continue;
				for(PHYNode n1 : fromLevelNodes) {
					ArrayList<PHYNode> toLevelNodes = nodes.get(j);
					if(toLevelNodes == null) continue;
					for(PHYNode n2: toLevelNodes) {
						checkAndAddEdge(n1, n2);
					}
				}
			}
		}
	}
	
	// ---- Graph Construction ----
	
	/**
	 * Checks if an edge should be added between two nodes in the network based on the AAF data.
	 * If yes, it adds the edge in the appropriate direction.
	 * @requires n1 to be at an equal or higher level than n2
	 * @param n1 - node 1
	 * @param n2 - node 2
	 */
	public void checkAndAddEdge(PHYNode n1, PHYNode n2) {
		if(n2.isLeaf) {
			int sampleId = n2.getLeafSampleId();
			if(n1.getAAF(sampleId) > 0) {
				addEdge(n1, n2);
			}
			return;
		}
		
		int comp = 0;
		for(int i = 0; i < numSamples; i++) {
			if((n1.getAAF(i) == 0) && (n2.getAAF(i) != 0)) break;
			comp += (n1.getAAF(i) >= (n2.getAAF(i) - AAF_ERROR_MARGIN)) ? 1 : 0;
		}
		
		if(comp == numSamples) {
			addEdge(n1, n2);
		} else {
			comp = 0;
			for(int i = 0; i < numSamples; i++) {
				if((n2.getAAF(i) == 0) && (n1.getAAF(i) != 0)) break;
				comp += (n2.getAAF(i) >= (n1.getAAF(i) - AAF_ERROR_MARGIN)) ? 1 : 0;
			}
			if(comp == numSamples) {
				addEdge(n2, n1);
			}
		} 
	}
	
	/** Adds a new node to the graph */
	public void addNode(PHYNode node, int level) {
		ArrayList<PHYNode> nodeList = nodes.get(level);
		if(nodeList == null) {
			nodes.put(level, new ArrayList<PHYNode>());
		}
		nodes.get(level).add(node);
		numNodes++;
	}
	
	/** Adds a new edge to the graph */
	public void addEdge(PHYNode from, PHYNode to) {
		ArrayList<PHYNode> nbrs = edges.get(from);
		if(nbrs == null) {
			edges.put(from, new ArrayList<PHYNode>());
		}
		edges.get(from).add(to);
		numEdges++;
		
		
	}
	
	/**
	 * Returns the nodes to which this node has edges to
	 * or null if the node connects to no edges
	 */
	public ArrayList<PHYNode> getDownNeighbors(PHYNode node) {
		return edges.get(node);
	}
	
	
	// ---- Spanning Tree Generation ----
	// based on Gabow, Myers '78
	
	public void grow() {
		
	}
	
	
	
	
	/**
	 * Place-holder operation for lineage
	 */
	public void getLineageTrees() {
		
	}
	
	/**
	 * Returns a string representation of the graph
	 */
	public String toString() {
		String graph = "--- PHYLOGENETIC CONSTRAINT GRAPH --- \n";
		graph += "numNodes = " + numNodes + ", ";
		graph += "numEdges = " + numEdges + "\n";
		
		// print nodes by level
		graph += "NODES: \n";
		for(int i = numSamples + 1; i >= 0; i--) {
			graph += "level = " + i + ": \n";
			ArrayList<PHYNode> levelNodes = nodes.get(i);
			if(levelNodes == null) {
				graph += "EMPTY \n";
			} else {
				for(PHYNode n : levelNodes) {
					graph += n.toString() + "\n";
				}
			}
		}
		graph += "EDGES: \n";
		for(PHYNode n1 : edges.keySet()) {
			ArrayList<PHYNode> nbrs = edges.get(n1);
			for(PHYNode n2 : nbrs) {
				graph += n1.getNodeId() + " -> " + n2.getNodeId() + "\n";
			}
		}
		
		return graph;
	}
	
	/**
	 * Node in the phylogenetic graph
	 * Represents a sub-population or sample (if leaf)
	 * Is associated with a given AAF (alternative allele frequency)
	 * 
	 */
	protected class PHYNode {
		
		/** Sub-population cluster that the node represents */
		private Cluster cluster;		
		
		/** SNP group the node belongs to */
		private SNPGroup snpGroup;
		
		/** Flag indicating if the node is a sample leaf*/
		private boolean isLeaf;
		
		/** Flag indicating if the node is the germline root */
		private boolean isRoot;
		
		/** The sample id if node is a leaf*/
		private int leafSampleId;
		
		/** Debugging-only node id */
		private int nodeId;
		
		
		
		/** 
		 * Internal node constructor
		 * @param g - SNP group the node belongs to
		 * @param nodeClusterId
		 */
		public PHYNode(SNPGroup g, int nodeClusterId) {
			snpGroup = g;
			cluster = snpGroup.getSubPopulations()[nodeClusterId];
			isLeaf = false;
			nodeId = nodeCounter;
			nodeCounter++;
		}
		
		/**
		 * Leaf node constructor - represents each tumor sample
		 * @param sampleId - ID of the represented tumor sample
		 */
		public PHYNode(int sampleId) {
			isLeaf = true;
			leafSampleId = sampleId;
			nodeId = nodeCounter;
			nodeCounter++;
		}
		
		/**
		 * Root node constructor
		 */
		public PHYNode() {
			isRoot = true;
			nodeId = nodeCounter;
			nodeCounter++;
		}
		
		/**
		 * Returns the SNP cluster that this node represents
		 */
		public Cluster getCluster() {
			return cluster;
		}
		
		/**
		 * Returns true if the node is a leaf
		 */
		public boolean isLeaf() {
			return isLeaf;
		}
		
		/**
		 * Returns the ID of the sample 
		 * @requires node is a leaf
		 */
		public int getLeafSampleId() {
			return leafSampleId;
		}
		
		public int getNodeId() {
			return nodeId;
		}
		
		/**
		 * Returns the cluster centroid AAF for the given sample ID
		 * Returns 0 if the sample is not represented
		 */
		public double getAAF(int sampleId) {
			if(isRoot) {
				return AAF_MAX;
			} 
			if(isLeaf) {
				return 0;
			}
			
			int sampleIndex = snpGroup.getSampleIndex(sampleId);
			if(sampleIndex == -1) {
				return 0;
			}
			return cluster.getCentroid()[sampleIndex];
		}
		
		public String toString() {
			String node = "Node " + nodeId + ": ";
			if(!isLeaf && !isRoot) {
				node += "group tag = " + snpGroup.getTag() + ", ";
				node += cluster.toString();
			} else if(isLeaf) {
				node += "leaf sample id = " + leafSampleId;
			} else {
				node += "root";
			}
			return node;
		}
		
	}
	
	
}
