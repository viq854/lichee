package lineage;

public class Parameters {
	
	// Input type
	protected enum Format { VCF, SNV, SNV_WITH_PROFILE}
	protected static Format INPUT_FORMAT = Format.SNV;
	protected static boolean CP = false;
	
	// SNV group partitioning
	protected static double MAX_ALLOWED_VAF = 0.6;
	protected static double MIN_VAF_PRESENT = 0.005; 
	protected static double MAX_VAF_ABSENT = 0.005;
	protected static int MIN_SNVS_PER_GROUP = 1;  
	protected static int MIN_ROBUST_SNVS_PER_GROUP = 0;  
	protected static double MIN_GROUP_PROFILE_SUPPORT = 2;  
	protected static double MIN_VAF_TARGET_RATIO_PER_SAMPLE = 0.5;
	
	// Clusters
	/** Minimum size a sub-population cluster must have to be a considered a node in the network */
	protected static int MIN_CLUSTER_SIZE = 2;
	protected static int MIN_PRIVATE_CLUSTER_SIZE = 1; 
	protected static double MIN_ROBUST_CLUSTER_SUPPORT = 2;
	
	/** Maximum centroid difference up to which two clusters can be collapsed */
	protected static double MAX_COLLAPSE_CLUSTER_DIFF = 0.2; 
	
	// Constraint graph and spanning tree generation
	/** Maximum VAF (used for the root node) */
	protected static double VAF_MAX = 0.5;
	/** Error margin used for comparing VAF centroid values when adding edges in the network */
	protected static double VAF_ERROR_MARGIN = 0.1;	
	protected static boolean ALL_EDGES = false;
	
	/** Stop tree search once this many valid trees have been found */
	protected static int MAX_NUM_TREES = 100000;
	protected static int MAX_NUM_GROW_CALLS = 100000000;
	protected static int NUM_TREES_FOR_CONSISTENCY_CHECK = 0;
}
