package lineage;

import io.VCFEntry;
import java.util.ArrayList;
import lineage.AAFClusterer.Cluster;

/**
 * A SNP group is a set of SNPs occurring in a given subset of samples.
 * All SNPs are partitioned into SNP groups based on their occurrence across samples.
 * Given S samples, there can be at most 2^S different groups.
 * A SNP group is uniquely identified by an S-bit binary tag (each bit corresponding to
 * a given sample), where a bit is set if that sample contains the SNPs in this group.
 *
 */
public class SNPGroup {

	/** Binary tag identifying the group 
	 * (the length of the tag is equal to the number of input samples) */
	private String tag;
	
	/** Number of samples represented by this group */
	private int numSamples;
	
	/** Number of SNPs assigned to this group */
	private int numSNPs;
	
	/** Alternative allele frequency data matrix (numSNPs x numSamples) */
	private double[][] alleleFreqBySample;
	
	/** SubPopulation clusters */
	private Cluster[] subPopulations;
	
	public SNPGroup(String groupTag, ArrayList<VCFEntry> snps) {
		tag = groupTag;
		numSamples = 0;		
		int[] sampleIndex = new int[tag.length()];
		for(int i = 0; i < tag.length(); i++) {
			if(tag.charAt(i) == '1') {
				sampleIndex[numSamples] = i;
				numSamples++;
			}
		}
		numSNPs = snps.size();
		alleleFreqBySample = new double[numSNPs][numSamples];
		for(int i = 0; i < numSNPs; i++) {
			VCFEntry snp = snps.get(i);
			for(int j = 0; j < numSamples; j++) {
				alleleFreqBySample[i][j] = snp.getAAF(sampleIndex[j]);
			}
		}
	}
	
	public double[][] getAlleleFreqBySample() {
		return alleleFreqBySample;
	}
	
	public int getNumSamples() {
		return numSamples;
	}
	
	public int getNumSNPs() {
		return numSNPs;
	}
	
	public void setSubPopulations(Cluster[] clusters) {
		subPopulations = clusters;
	}
	
	public Cluster[] getSubPopulations() {
		return subPopulations;
	}
}
