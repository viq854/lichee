package io;

public final class VCFConstants {

	private VCFConstants(){}
	
	public static int NormalSample = 1;
	
	// to validate SNP
	public static final double AVG_COVERAGE = 50;
	public static final double MIN_COVERAGE = 10;
	public static final double GROUP_PVALUE = 0.5;
	public static final double MIN_QUAL = 200;
	
	// to edit SNP 
	public static final double BASE_ERROR = 0.02;
	public static final double EDIT_PVALUE= 0.1;
	public static final int EDIT_DISTANCE = 4;
	
	// to identify sub-population _
	public static final double SUBPOP_PVALUE = 0.01;
	
	// LOH analysis
	public static final double HETEROZYGOUS = 0.2;
	public static final double LOH_ERROR = 0.1;
	public static final int LOH_WINDOW = 20;
	public static final int SMOOTHING_WINDOW = 20;
	
	
}
