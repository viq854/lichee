package lineage;

public class Parameters {
	
	// Input type
	protected enum Format { VCF, MUT, MUTC}
	protected static Format INFORMAT = Format.MUT;
	protected static boolean CP = false;
	
	// SNV group partitioning
	protected static double MAX_ALLOWED_VAF = 0.6;
	protected static double VALIDATION_THR = 0.005; 
	protected static double VALIDATION_SOFT_THR = 0.005;
	protected static int GROUP_SIZE_THR = 1;  
	protected static int GROUP_ROBUST_NUM_THR = 0;  
	protected static double ROBUSTGROUP_SIZE_THR = 2; 
	protected static double MIN_VAF_TARGET_RATIO_PER_SAMPLE = 0.5;
	protected static int MIN_ALLOWED_COVERAGE = 500; 
	
	// Clusters
	/** Minimum size a sub-population cluster must have to be a considered a node in the network */
	protected static int MIN_CLUSTER_SIZE = 2;
	protected static int MIN_PRIVATE_CLUSTER_SIZE = 1; 
	
	/** Maximum centroid difference up to which two clusters can be collapsed */
	protected static double MAX_COLLAPSE_CLUSTER_DIFF = 0.2; 
	
	// Constraint graph and spanning tree generation
	/** Maximum AAF (used for the root node) */
	protected static double AAF_MAX = 0.5;
	/** Error margin used for comparing AAF centroid values when adding edges in the network */
	protected static double AAF_ERROR_MARGIN = 0.1;	
	protected static boolean ALL_EDGES = false;
	
	/** Stop tree search once this many valid trees have been found */
	protected static int MAX_NUM_TREES = 100000;
	protected static int MAX_NUM_GROW_CALLS = 100000000;
	protected static int NUM_TREES_FOR_CONSISTENCY_CHECK = 5;
	protected static boolean CHECK_CONSISTENCY = true;
}
