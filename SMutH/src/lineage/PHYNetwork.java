package lineage;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import util.SNVEntry;
import util.Visualizer;

/**
 * PHYNetwork is a directed constraint graph representing the phylogenetic relationship
 * among sample sub-populations. It is a DAG.
 * Each internal node in the graph represents a sub-population.
 * A directed edge between two nodes denotes the 'happened-before' evolutionary 
 * relationship between the two nodes.
 */
public class PHYNetwork implements Serializable {
	private static final long serialVersionUID = 1L;
	
	/** Nodes in the graph divided by levels (number of samples node SNVs occurred in) */
	private HashMap<Integer, ArrayList<PHYNode>> nodes;
	
	/** Nodes in the graph indexed by their unique ID */
	private transient HashMap<Integer, PHYNode> nodesById;
	
	/** Adjacency map of nodes to the their neighbors/children */
	private transient HashMap<PHYNode, ArrayList<PHYNode>> edges;
	
	/** Total number of nodes in the graph.
	 *  During construction: used as a counter to assign unique IDs to nodes */
	protected int numNodes;
	
	/** Total number of edges in the graph */
	private int numEdges;
	
	/** Total number of tissue samples */
	protected int numSamples;
	
	/** Maximum AAF (used for the root node) */
	protected static final double AAF_MAX = 0.5;
		
	private static Logger logger = LineageEngine.logger;
	
	// ---- Network Construction ----
	
	/**
	 * Constructs a PHYNetwork from the sub-populations of the SNV groups
	 */
	public PHYNetwork(ArrayList<SNVGroup> groups, int totalNumSamples) {
		numSamples = totalNumSamples;
		numNodes = 0;
		nodes = new HashMap<Integer, ArrayList<PHYNode>>();
		nodesById = new HashMap<Integer, PHYNode>();
		edges = new HashMap<PHYNode, ArrayList<PHYNode>>(); 
	
		// add root node
		PHYNode root = new PHYNode(numSamples+1, numNodes);
		addNode(root, numSamples+1);
				
		// add group sub-population nodes
		for(SNVGroup g : groups) {
			PHYNode[] groupNodes = new PHYNode[g.getSubPopulations().length];
			for(int i = 0; i < groupNodes.length; i++) {
				PHYNode node = new PHYNode(g, i, g.getNumSamples(), numNodes);
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
		
		// add inter-level edges
		for(int i = numSamples + 1; i > 0; i--) {
			ArrayList<PHYNode> fromLevelNodes = nodes.get(i);
			if(fromLevelNodes == null) continue;
			// find the next non-empty level
			int j = i-1;
			ArrayList<PHYNode> toLevelNodes = nodes.get(j);
			while((toLevelNodes == null) && (j > 0)) {
				j--;
				toLevelNodes = nodes.get(j);
			}
			if(toLevelNodes == null) continue;
			for(PHYNode n1 : fromLevelNodes) {
				for(PHYNode n2: toLevelNodes) {
					checkAndAddEdge(n1, n2);
				}
			}
		}
		
		// find the nodes that are not connected and connect them to a valid node in the closest
		// higher level
		int[] nodeMask = new int[numNodes];
		for(PHYNode n : edges.keySet()) {
			for(PHYNode m : edges.get(n)) {
				nodeMask[m.getNodeId()] = 1;
			}
		}
		
		// skips the root
		for(int i = 1; i < nodeMask.length; i++) {
			if(nodeMask[i] == 0) {
				PHYNode n = nodesById.get(i);			
				// find a parent in the closest higher level
				boolean found = false;
				for(int j = n.getLevel() + 2; j <= numSamples + 1; j++) {
					ArrayList<PHYNode> fromLevelNodes = nodes.get(j);
					if(fromLevelNodes == null) continue;
					for(PHYNode n2 : fromLevelNodes) {
						if(checkAndAddEdge(n2, n) == 0) {
							// found a parent
							found = true;
							break;
						}
					}
					if(found) break;
				}
				if(!found) {
					addEdge(root, n);
				}
			}
		}
	}
	
	/**
	 * Checks if an edge should be added between two nodes in the network based on the AAF data.
	 * If yes, it adds the edge in the appropriate direction.
	 * The edge is added in the direction that minimizes the error
	 * @requires n1 to be at an equal or higher level than n2
	 * @param n1 - node 1
	 * @param n2 - node 2
	 */
	public int checkAndAddEdge(PHYNode n1, PHYNode n2) {
		if(n2.isLeaf) {
			int sampleId = n2.getLeafSampleId();
			if(n1.getAAF(sampleId) > 0) {
				addEdge(n1, n2);
				return 0;
			}
			return -1;
		}
		
		int comp_12 = 0;
		int comp_21 = 0;
		double err_12 = 0;
		double err_21 = 0;
		
		for(int i = 0; i < numSamples; i++) {
			if((n1.getAAF(i) == 0) && (n2.getAAF(i) != 0)) break;
			comp_12 += (n1.getAAF(i) >= (n2.getAAF(i) - getAAFErrorMargin(n1, n2, i))) ? 1 : 0;
			if(n1.getAAF(i) < n2.getAAF(i)) {
				err_12 += n2.getAAF(i) - n1.getAAF(i);
			}
		}
		for(int i = 0; i < numSamples; i++) {
			if((n2.getAAF(i) == 0) && (n1.getAAF(i) != 0)) break;
			comp_21 += (n2.getAAF(i) >= (n1.getAAF(i) - getAAFErrorMargin(n2, n1, i))) ? 1 : 0;
			if(n2.getAAF(i) < n1.getAAF(i)) {
				err_21 += n1.getAAF(i) - n2.getAAF(i);
			}
		}
		
		if(comp_12 == numSamples) {
			if (comp_21 == numSamples) {
				if(err_12 < err_21) {
					addEdge(n1, n2);
					return 0;
				} else {
					addEdge(n2, n1);
					return 1;
				}
			} else {
				addEdge(n1, n2);
				return 0;
			}
		} else if(comp_21 == numSamples) {
			addEdge(n2, n1);
			return 1;
		}
		
		return -1;
	}
	
	/**
	 * Returns the AAF error margin on the edge 
	 * between the from and to nodes
	 */
	private double getAAFErrorMargin(PHYNode from, PHYNode to, int i) {
		if(Parameters.STATIC_ERROR_MARGIN) {
			return Parameters.AAF_ERROR_MARGIN;
		}
		
		double parentStdError;
		double childStdError;
		int parentSampleSize = 0;
		int childSampleSize = 0;
		if(from.isRoot()) {
			parentStdError = Parameters.AAF_ERROR_MARGIN;
		} else {
			parentSampleSize = from.getCluster().getMembership().size();
			parentStdError = 1.96*from.getStdDev(i)/Math.sqrt((double)parentSampleSize);
		}	
		if(to.isRoot()) {
			childStdError = Parameters.AAF_ERROR_MARGIN;
		} else {
			childSampleSize = to.getCluster().getMembership().size();
			childStdError = 1.96*to.getStdDev(i)/Math.sqrt((double)childSampleSize);
		}
		
		double standardError = parentStdError + childStdError;
		
		if(standardError > Parameters.AAF_ERROR_MARGIN)
			return standardError;
		
		return Parameters.AAF_ERROR_MARGIN;
		
		// max std dev of the centroid vector
		/*double maxStd = Parameters.AAF_ERROR_MARGIN;
		if(from.getStdDev(i) > maxStd) {
			maxStd = from.getStdDev(i);
		}
		if(to.getStdDev(i) > maxStd) {
			maxStd = to.getStdDev(i);
		}
		return maxStd; */
	}
	
	/** Adds a new node to the graph */
	public void addNode(PHYNode node, int level) {
		ArrayList<PHYNode> nodeList = nodes.get(level);
		if(nodeList == null) {
			nodes.put(level, new ArrayList<PHYNode>());
		}
		nodes.get(level).add(node);
		nodesById.put(node.getNodeId(), node);
		numNodes++;
	}
	
	/** Adds a new edge to the graph */
	public void addEdge(PHYNode from, PHYNode to) {
		ArrayList<PHYNode> nbrs = edges.get(from);
		if(nbrs == null) {
			edges.put(from, new ArrayList<PHYNode>());
		}
		if(!edges.get(from).contains(to)) {
			edges.get(from).add(to);
			numEdges++;
		}
	}
	
	/** Removes an edge from the graph */
	public void removeEdge(PHYNode from, PHYNode to) {
		ArrayList<PHYNode> nbrs = edges.get(from);
		if(nbrs != null) {
			for(PHYNode n : nbrs) {
				if(n.equals(to)) {
					nbrs.remove(n);
					break;
				}
			}
		}
	}
	
	/** Adds all the inter-level edges */
	public void addAllHiddenEdges() {
		for(int i = numSamples + 1; i > 0; i--) {
			ArrayList<PHYNode> fromLevelNodes = nodes.get(i);
			if(fromLevelNodes == null) continue;
			for(int j = i-1; j > 0; j--) {
				ArrayList<PHYNode> toLevelNodes = nodes.get(j);
				if(toLevelNodes == null) continue;
				for(PHYNode n1 : fromLevelNodes) {
					for(PHYNode n2: toLevelNodes) {
						checkAndAddEdge(n1, n2);
					}
				}
			}
		}
	}
	
	// ---- Network Adjustments ----
	
	/**
	 * The network needs to be adjusted when no valid spanning PHYTrees are found.
	 * Adjustments include: 
	 * - removing nodes corresponding to groups that are less robust
	 * - increasing the error margin (to do)
	 * - adding hidden edges (to do)
	 */
	public PHYNetwork fixNetwork() {
		// reconstruct the network from robust groups only
		Set<SNVGroup> filteredGroups = new HashSet<SNVGroup>();
		SNVGroup toRemove = null;
		for(PHYNode n : nodesById.values()) {
			SNVGroup group = n.snvGroup;
			if(group != null) {
				filteredGroups.add(n.snvGroup);
				if((!group.isRobust())) {
					if (toRemove == null) {
						toRemove = group;
					} else if(group.getNumSNVs() < toRemove.getNumSNVs()) {
						toRemove = group;
					}
				}
			}
		}
		if(toRemove != null) {
			filteredGroups.remove(toRemove);
			logger.log(Level.INFO, "Removed group " + toRemove.getTag() + " of size " + toRemove.getNumSNVs());
		}
		return new PHYNetwork(new ArrayList<SNVGroup>(filteredGroups), numSamples);
	}
	
	// ---- Spanning PHYTree Generation ----
	
	// based on the algorithm from Gabow & Myers '78
	
	/** List of all generated spanning trees */
	private transient ArrayList<PHYTree> spanningTrees;
	
	/** Stack of edges directed from vertices in tree T to vertices not in T */
	private transient ArrayList<PHYEdge> f;
	
	/** The last spanning tree output so far */
	private transient PHYTree L;
	
	/**
	 * Finds all spanning trees rooted at r
	 */
	public void grow(PHYTree t) {
		// if the tree t contains all the nodes, it is complete
		if(t.treeNodes.size() == numNodes) {
			L = t;
			spanningTrees.add(L.clone());
		} else {
			// list used to reconstruct the original F
			ArrayList<PHYEdge> ff = new ArrayList<PHYEdge>();
			
			boolean b = false;
			while(!b && (f.size() > 0)) {
				// new tree edge
				PHYEdge e = f.remove(f.size() - 1);
				PHYNode v = e.to;
				t.addNode(v);
				t.addEdge(e.from, v);
				
				// update f
				ArrayList<PHYEdge> edgesAdded = new ArrayList<PHYEdge>();
				ArrayList<PHYNode> vNbrs = edges.get(v);
				if(vNbrs != null) {
					for(PHYNode w : vNbrs) {
						if(!t.containsNode(w)) {
							PHYEdge vw = new PHYEdge(v, w);
							f.add(vw);
							edgesAdded.add(vw);
						}
					}
				}
				
				// remove (w,v) w in T from f
				ArrayList<PHYEdge> edgesRemoved = new ArrayList<PHYEdge>();
				for(int i = 0; i < f.size(); i++) {
					PHYEdge wv = f.get(i);
					if(t.containsNode(wv.from) && (wv.to.equals(v))) {
						edgesRemoved.add(wv);
					}
				}
				f.removeAll(edgesRemoved);
	
				// recurse
				grow(t);
				
				// pop
				f.removeAll(edgesAdded);
				
				// restore
				f.addAll(edgesRemoved);
				
				// remove e from T and G
				t.removeEdge(e.from, e.to);
				this.removeEdge(e.from, e.to);
				
				// add e to FF
				ff.add(e);
				
				// bridge test
				for(PHYNode w : this.edges.keySet()) {
					ArrayList<PHYNode> wNbrs = this.edges.get(w);
					if(wNbrs == null) continue;
					for(PHYNode n : wNbrs) {
						if(n.equals(v)) {
							// check if w is a descendant of v in L
							if(!L.isDescendent(v, w)) {
								b = false;
								break;
							}
						}
					}
					if(!b) break;
				}
			}
			
			// pop from ff, push to f, add to G
			for(int i = ff.size()-1; i >=0; i--) {
				PHYEdge e = ff.get(i);
				f.add(e);
				this.addEdge(e.from, e.to);
			}
			ff.clear();
		}
	}
	
	/**
	 * Generates all the spanning trees from the constraint network
	 * that pass the AAF constraints
	 */
	public ArrayList<PHYTree> getLineageTrees() {
		spanningTrees = new ArrayList<PHYTree>();
		
		PHYNode root = nodes.get(numSamples+1).get(0);
		// initialize tree t to contain the root
		PHYTree t = new PHYTree();
		t.addNode(root);
		// initialize f to contain all edges (root, v)
		f = new ArrayList<PHYEdge>();
		for(PHYNode n : edges.get(root)) {
			f.add(new PHYEdge(root, n));
			
		}
		grow(t);
		applyAAFConstraints(spanningTrees);
		return spanningTrees;
	}
	
	/**
	 * Applies the AAF constraints to all the spanning trees
	 * and removes the trees that don't pass the constraints
	 */
	private void applyAAFConstraints(ArrayList<PHYTree> trees) {
		ArrayList<PHYTree> toBeRemoved = new ArrayList<PHYTree>();
		for(PHYTree t : trees) {
			if(!checkAAFConstraints(t)) {
				toBeRemoved.add(t);
			}
		}
		spanningTrees.removeAll(toBeRemoved);
	}
	
	/**
	 * Returns true if the tree passes the AAF constraints
	 * @param t - spanning tree
	 */
	private boolean checkAAFConstraints(PHYTree t) {
		for(PHYNode n : t.treeEdges.keySet()) {
			ArrayList<PHYNode> nbrs = t.treeEdges.get(n);			
			for(int i = 0; i < numSamples; i++) {
				double affSum = 0;
				double errMargin = 0.0;
				for(PHYNode n2 : nbrs) {
					affSum += n2.getAAF(i);
					errMargin += getAAFErrorMargin(n, n2, i);
				}
				if(affSum >= n.getAAF(i) + errMargin) {
					return false;
				}
			}
		}
		return true;
	}
	
	/** 
	 * Evaluates the spanning trees by computing their error score
	 * and ranking them by this score (lowest error first)
	 */
	public void evaluateLineageTrees() {
		Collections.sort(spanningTrees);
	}
	
	/** Debugging only - tests that all the spanning trees found are different */
	protected void testSpanningTrees() {
		for(int i = 0; i < spanningTrees.size(); i++) {
			for(int j = i + 1; j < spanningTrees.size(); j++) {
				PHYTree t1 = spanningTrees.get(i);
				PHYTree t2 = spanningTrees.get(j);
				
				// compare
				boolean sameEdges = true;
				for(PHYNode n1 : t1.treeEdges.keySet()) {
					for(PHYNode n2 : t1.treeEdges.get(n1)) {
						sameEdges &= t2.containsEdge(n1, n2);
					}
				}
				if(sameEdges) {
					System.out.println("Found same tree");
					System.out.println(t1);
					System.out.println(t2);
					return;
				}
			}
		}
		System.out.println("All trees are distinct");
	}
	
	// ---- Visualization ----
	
	/** Displays the constraint network graph */
	public void displayNetwork() {
		DirectedGraph<Integer, Integer> g = new DirectedSparseGraph<Integer, Integer>();
		HashMap<Integer, String> nodeLabels = new HashMap<Integer, String>();
			
		int edgeId = 0;
		for (PHYNode n : edges.keySet()) {
			g.addVertex(n.getNodeId());
			nodeLabels.put(n.getNodeId(), n.getLabel());
			for(PHYNode n2 : edges.get(n)) {
				if(!g.containsVertex(n2.getNodeId())) {
					g.addVertex(n2.getNodeId());
					nodeLabels.put(n2.getNodeId(), n2.getLabel());
				}
				g.addEdge(edgeId, n.getNodeId(), n2.getNodeId(), EdgeType.DIRECTED);
				edgeId++;
			}
		}
		Visualizer.showNetwork(g, nodeLabels);	
	}
	
	/** Displays a spanning tree of the network */
	public void displayTree(PHYTree t, String[] sampleNames, HashMap<String, ArrayList<SNVEntry>> snvsByTag, String fileOutputName) {			
		DirectedGraph<Integer, Integer> g = new DirectedSparseGraph<Integer, Integer>();
		HashMap<Integer, String> nodeLabels = new HashMap<Integer, String>();
		HashMap<Integer, PHYNode> nodeObj = new HashMap<Integer, PHYNode>();
		
		int edgeId = 0;
		for (PHYNode n : t.treeEdges.keySet()) {
			g.addVertex(n.getNodeId());
			nodeLabels.put(n.getNodeId(), n.getLabel());
			nodeObj.put(n.getNodeId(), n);
			for(PHYNode n2 : t.treeEdges.get(n)) {
				if(!g.containsVertex(n2.getNodeId())) {
					g.addVertex(n2.getNodeId());
					nodeLabels.put(n2.getNodeId(), n2.getLabel());
					nodeObj.put(n2.getNodeId(), n2);
				}
				g.addEdge(edgeId, n.getNodeId(), n2.getNodeId(), EdgeType.DIRECTED);
				edgeId++;
			}
		}
		
		// add sample leaves
		for(int i = 0; i < numSamples; i++) {
			PHYNode n = new PHYNode(0, i, numNodes + i);
			g.addVertex(-n.getNodeId());
			nodeLabels.put(-n.getNodeId(), sampleNames[i]);
			nodeObj.put(-n.getNodeId(), n);
			
			// find a parent in the closest higher level		 
			boolean found = false;
			ArrayList<PHYNode> parents = new ArrayList<PHYNode>();
			ArrayList<PHYNode> sameLevelParents = new ArrayList<PHYNode>();
			for(int j = n.getLevel() + 1; j <= numSamples; j++) {
				ArrayList<PHYNode> fromLevelNodes = nodes.get(j);
				if(fromLevelNodes == null) continue;
				for(PHYNode n2 : fromLevelNodes) {
					if(n2.getAAF(i) > 0) {
						boolean addEdge = true;
						for(PHYNode p : parents) {
							if(t.isDescendent(n2, p)) {
								addEdge = false;
								break;
							}
						}
						if(addEdge) {
							sameLevelParents.add(n2);
							parents.add(n2);
							found = true;
						}
					}
				}
				// remove nodes that are in same level that are connected
				ArrayList<PHYNode> toRemove = new ArrayList<PHYNode>();
				for(PHYNode n1 : sameLevelParents) {
					for(PHYNode n2 : sameLevelParents) {
						if(t.isDescendent(n1, n2)) {
							toRemove.add(n1);
						}
					}
				}
				sameLevelParents.removeAll(toRemove);
				
				for(PHYNode n2 : sameLevelParents) {
					g.addEdge(edgeId, n2.getNodeId(), -n.getNodeId());
					edgeId++;
				}
				sameLevelParents.clear();
			}
			if(!found) {
				g.addEdge(edgeId, 0, -n.getNodeId());
				edgeId++;
			}
		}			
		Visualizer.showLineageTree(g, nodeLabels, snvsByTag, fileOutputName, nodeObj, t);	
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
			if(levelNodes != null) {
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
	 * Returns a string representation of the graph
	 */
	public String getNodesAsString() {
		String s = "";
		
		// print nodes by level
		for(int i = numSamples + 1; i >= 0; i--) {
			ArrayList<PHYNode> levelNodes = nodes.get(i);
			if(levelNodes != null) {
				for(PHYNode n : levelNodes) {
					if(n.isRoot()) continue;
					s += n.getSNVGroup().getTag() + "\t";
					s += n.getCluster().getMembership().size() + "\t";
					double[] c = n.getCluster().getCentroid();
					DecimalFormat df = new DecimalFormat("#.##");
					for(int j = 0; j < c.length; j++) {
						s += df.format(c[j]) + "\t";
					}
					s += "\n";
				}
			}
		}
		return s;
	}
}
