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
	
	/** Total number of cancer samples */
	private int numSamples;
	
	/** Error margin used for comparing AAF centroid data */
	private static final double AAF_ERROR_MARGIN = 0.08;
	
	/**
	 * Constructs a PHYGraph from the sub-populations of the SNP groups
	 */
	public PHYGraph(ArrayList<SNPGroup> groups, int totalNumSamples) {
		numSamples = totalNumSamples;
		nodes = new HashMap<Integer, ArrayList<PHYNode>>();
		edges = new HashMap<PHYNode, ArrayList<PHYNode>>(); 
		
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
		for(int i = numSamples; i > 0; i--) {
			for(int j = i-1; j >= 0; j--) {
				for(PHYNode n1 : nodes.get(i)) {
					for(PHYNode n2: nodes.get(j)) {
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
			comp += (n1.getAAF(i) >= (n2.getAAF(i) - AAF_ERROR_MARGIN)) ? 1 : 0;
		}
		if(comp == numSamples) {
			addEdge(n1, n2);
		} else if (comp == 0) {
			addEdge(n2, n1);
		} 
	}
	
	/** Adds a new node to the graph */
	public void addNode(PHYNode node, int level) {
		ArrayList<PHYNode> nodeList = nodes.get(level);
		if(nodeList == null) {
			nodes.put(level, new ArrayList<PHYNode>());
		}
		nodes.get(level).add(node);
	}
	
	/** Adds a new edge to the graph */
	public void addEdge(PHYNode from, PHYNode to) {
		ArrayList<PHYNode> nbrs = edges.get(from);
		if(nbrs == null) {
			edges.put(from, new ArrayList<PHYNode>());
		}
		edges.get(from).add(to);
	}
	
	/**
	 * Returns the nodes to which this node has edges to
	 * or null if the node connects to no edges
	 */
	public ArrayList<PHYNode> getDownNeighbors(PHYNode node) {
		return edges.get(node);
	}
	
	
	// ---- Network Flow ----
	
	/**
	 * Place-holder operation for lineage
	 */
	public void getLineageTrees() {
		
	}
	
	/**
	 * Node in the phylogenetic graph
	 * Represents a sub-population or sample (if leaf)
	 * Is associated with a given AAF (alternative allele frequency)
	 * 
	 */
	protected class PHYNode {
		private Cluster cluster;
		private SNPGroup snpGroup;
		private boolean isLeaf;
		private int leafSampleId;
		
		public PHYNode(SNPGroup g, int nodeClusterId) {
			snpGroup = g;
			cluster = snpGroup.getSubPopulations()[nodeClusterId];
			isLeaf = false;
		}
		
		// leaf
		public PHYNode(int sampleId) {
			isLeaf = true;
			leafSampleId = sampleId;
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
		
		public int getLeafSampleId() {
			return leafSampleId;
		}
		
		/**
		 * Returns the cluster centroid AAF for the given sample ID
		 */
		public double getAAF(int sampleId) {
			int sampleIndex = snpGroup.getSampleIndex(sampleId);
			if(sampleIndex == -1) {
				return 0;
			}
			return cluster.getCentroid()[sampleIndex];
		}
		
	}
	
	
}
