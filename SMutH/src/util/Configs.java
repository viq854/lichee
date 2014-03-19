package util;


public final class Configs {

	private Configs(){}
	
	// general variables!!!
	public static String path;
	public static String testName;
		
	// Type of input
	public enum format { VCF, MUT, FL, SIM}
	public static format INFORMAT = format.MUT;
	
	// for validation SNVs
	public static double VALIDATION_THR = 0.04;
	public static double VALIDATION_SOFT_THR = 0.015; 
	public static double SYST_ERR_THR = 0.005;
	
	
	// for WGS SNVs
	//public static final double AVG_COVERAGE = 50;
	public static final double WG_HARD_THR = 0.1;
	public static final double WG_SOFT_THR = 0.01;
	public static final double MIN_COVERAGE = 14;
	public static final double MIN_QUAL = 30;
	
	// For group size; 1 == any size is acceptable
	public enum GroupSizeType {FIXED, AVERAGE, SIGNIFICANT}; 
	public static  GroupSizeType gst = GroupSizeType.FIXED;
	public static double GROUP_PVALUE = 0.2;
	
	public static int GROUP_SIZE_THR = 1;  
	public static int GROUP_ROBUST_NUM_THR = 0;  
	public static double ROBUSTGROUP_SIZE_THR = 2; 
	public static int EDIT_DISTANCE = 15;
	public static double MIN_VAF_TARGET_RATIO_PER_SAMPLE = 0.5;
	public static double MAX_ALLOWED_VAF = 0.6;
	public static int MIN_ALLOWED_COVERAGE = 200;
	
	
	public static final double SUBPOP_PVALUE = 0.0001;
	// to edit SNV 
	public static final double BASE_ERROR = 0.02;
	public static final double EDIT_PVALUE= 0.01;
	
	
	// LOH analysis
	public static final double HETEROZYGOUS = 0.2;
	public static final double LOH_ERROR = 0.1;
	public static final int LOH_WINDOW = 20;
	public static final int SMOOTHING_WINDOW = 20;
	

	/**
	 * Returns the minimum size for a group to be considered significant
	 * @param n total number of mutations
	 * @param k number of groups
	 * @param THR p-value threshold 
	 * @return
	 */
	
	public static int getGroupSizeThreshold(int n, int k){
		switch(gst){
		//a. fixed size
		case FIXED:
			return GROUP_SIZE_THR ;
		//b. simply n/(constant*k)
		case AVERAGE:
			return (int)(n/(10*k));
		//c. using p-value
		default :
			double pValue = 0.0;
			int x = n / k;
			for (; x > 1; x--){
				pValue = Math.pow((n - x * k)/(double) n, k - 1);
				if (pValue > GROUP_PVALUE) break;
			}
			//System.out.println("gourp size "+n+" "+k+" "+x);
			return x;
		}
	}
	
}
