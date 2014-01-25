package util;

public final class Configs {

	private Configs(){}
	
	// general variables!!!
	public static String path;
	public static String testName;
	public static int normalSample;
	
	public static boolean IGNORE_GERMLINE = Boolean.FALSE;
	
	// Type of input
	public enum format { VCF, MUT, FL}
	public static format INFORMAT = format.MUT;
	
	// for validation SNVs
	public static final double VALIDATION_THR = 0.04; //0.04 0.06 0.07 0.08
	public static final double VALIDATION_SOFT_THR = 0.015; // 0.015 0.02
	
	// for WGS SNVs
	//public static final double AVG_COVERAGE = 50;
	public static final double WG_HARD_THR = 0.2;
	public static final double WG_SOFT_THR = 0.01;
	public static final double MIN_COVERAGE = 14;
	public static final double MIN_QUAL = 30;
	
	// For group size; 1 == any size is acceptable
	public static final double GROUP_PVALUE = 0.1; 
	public static final double ROBUSTGROUP_PVALUE = 0.01;

	
	// to edit SNV 
	public static final double BASE_ERROR = 0.02;
	public static final double EDIT_PVALUE= 0.01;
	public static final int EDIT_DISTANCE = 2;
	
	// to identify sub-population _
	public static final double SUBPOP_PVALUE = 0.00001;
	
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
	
	public static int getGroupSizeThreshold(int n, int k, double THR){
		double pValue = 0.0;
		int x = n / k;
		for (; x > 1; x--){
			pValue = Math.pow((n - x * k)/(double) n, k - 1);
			if (pValue > THR) break;
		}
		//System.out.println("gourp size "+n+" "+k+" "+x);
		return x;
	}
	
}
