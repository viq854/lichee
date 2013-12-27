package lineage;

public class Parameters {
	
	// Constraint graph and spanning tree generation
	
	/** Error margin used for comparing AAF centroid values when adding edges in the network */
	protected static final double AAF_ERROR_MARGIN = 0.08;
	
	// Clusters
	
	/** Minimum size a sub-population cluster must have to be a considered a node in the network */
	protected static final int MIN_CLUSTER_SIZE = 2;
	
	/** Maximum centroid difference up to which two clusters can be collapsed */
	protected static final double MAX_COLLAPSE_CLUSTER_DIFF = 0;
	
	/** Maximum number of clusters allowed per SNV group (additional clusters will be collapsed) */
	protected static final int MAX_CLUSTER_NUM = 2;
	

}
