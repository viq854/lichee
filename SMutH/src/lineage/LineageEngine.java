package lineage;

/**
 * MaincCell lineage builder pipeline
 *
 */
public class LineageEngine {

	/**
	 * The main pipeline for constructing the cell lineage 
	 * @param vcfFileName - path to the input VCF file
	 */
	public void buildLineage(String vcfFileName) {
		// 1. load VCF
		
		// 2. handle normal cell contamination, CNVs, 
		//    determine the minimum number of clusters using LOH
		//    (+ any additional filtering, pre-processing)
		
		// 3. partition SNPs into appropriate groups
		
		// 4. cluster SNPs in each group
		
		// 5. incorporate CNVs and re-cluster
		
		// 6. construct constraint network
		
		// 7. build and output final cell lineage
		
	}
	
}
