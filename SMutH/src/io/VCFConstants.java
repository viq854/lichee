package io;

public final class VCFConstants {

	private VCFConstants(){}
	
	// to validate SNV
	public static final double VALIDATION_THR = 0.01;
	public static final double VALIDATION_HARD_THR = 0.025;
	
	//public static final double AVG_COVERAGE = 50;
	public static final double MIN_COVERAGE = 14;
	public static final double GROUP_PVALUE = 0.1;
	public static final double MIN_QUAL = 30;
	
	// to edit SNP 
	public static final double BASE_ERROR = 0.02;
	public static final double EDIT_PVALUE= 0.1;
	public static final int EDIT_DISTANCE = 4;
	
	// to identify sub-population _
	public static final double SUBPOP_PVALUE = 0.00001;
	
	// LOH analysis
	public static final double HETEROZYGOUS = 0.2;
	public static final double LOH_ERROR = 0.1;
	public static final int LOH_WINDOW = 20;
	public static final int SMOOTHING_WINDOW = 20;
	
	
}
