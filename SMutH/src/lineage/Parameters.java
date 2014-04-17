package lineage;

public class Parameters {
	
	// Constraint graph and spanning tree generation
	/** Error margin used for comparing AAF centroid values when adding edges in the network */
	protected static double AAF_ERROR_MARGIN = 0.01;
	protected static boolean STATIC_ERROR_MARGIN = false;
	
	// Clusters
	/** Minimum size a sub-population cluster must have to be a considered a node in the network */
	protected static int MIN_CLUSTER_SIZE = 5;
	
	/** Maximum centroid difference up to which two clusters can be collapsed */
	protected static double MAX_COLLAPSE_CLUSTER_DIFF = 0.1;
	
	protected static boolean ALL_EDGES = true;
	
	/** Stop tree search once this many valid trees have been found */
	protected static int MAX_NUM_TREES = 100000;
}
