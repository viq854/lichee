package lineage;

import io.*;
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
public class SNVGroup {

	/** Binary tag identifying the group 
	 * (the length of the tag is equal to the number of input samples) */
	private String tag;
	
	/** Number of samples represented by this group */
	private int numSamples;
	
	/** Indices of the samples represented in this group (from 0 to |tag|-1 MSF order) */
	private int[] sampleIndex;
	
	/** Alternative allele frequency data matrix (numSNPs x numSamples) */
	private double[][] alleleFreqBySample;
	
	/** SubPopulation clusters */
	private Cluster[] subPopulations;
	
	/** SNPs assigned to this group */
	private ArrayList<SNVEntry> snps;
	
	public SNVGroup(String groupTag, ArrayList<SNVEntry> groupSNPs) {
		tag = groupTag;
		numSamples = 0;		
		sampleIndex = new int[tag.length()];
		for(int i = 0; i < tag.length(); i++) {
			if(tag.charAt(i) == '1') {
				sampleIndex[numSamples] = i;
				numSamples++;
			}
		}
		snps = groupSNPs;
		alleleFreqBySample = new double[snps.size()][numSamples];
		for(int i = 0; i < snps.size(); i++) {
			SNVEntry snp = snps.get(i);
			for(int j = 0; j < numSamples; j++) {
				alleleFreqBySample[i][j] = snp.getAAF(sampleIndex[j]);
			}
		}
		
		System.out.println("Created group: " + this.toString());
	}
	
	/**
	 * Insert additional SNPs into the group
	 * This method will add the SNP to the sub-population that is more likely 
	 * to contain the SNP
	 * The sub-population is picked by finding the closest centroid to the scaled
	 * AAF of the SNP 
	 */
	public void addSNPsFromCNVs(ArrayList<VCFEntry> cnvs, int[] scale) {
		snps.addAll(cnvs);		
		for(VCFEntry snp : cnvs) {
			double[] aaf = new double[numSamples];
			for(int j = 0; j < numSamples; j++) {
				aaf[j] = snp.getAAF(sampleIndex[j]);
				
				// scale AAF and find the cluster to which this point has the smallest distance to
			}
		}
	}
	
	// Getters/Setters
	
	public double[][] getAlleleFreqBySample() {
		return alleleFreqBySample;
	}
	
	public int getNumSamples() {
		return numSamples;
	}
	
	public int getNumSNPs() {
		return snps.size();
	}
	
	public Cluster[] getSubPopulations() {
		return subPopulations;
	}
	
	public String getTag() {
		return tag;
	}
	
	/**
	 * Returns the index of this sample in the centroid/AAF data of the group
	 * @return -1 if this sample is not represented in the group
	 */
	public int getSampleIndex(int sampleId) {
		for(int i = 0; i < numSamples; i++) {
			if(sampleIndex[i] == sampleId) {
				return i;
			}
		}
		return -1;
	}
	
	public void setSubPopulations(Cluster[] clusters) {
		subPopulations = clusters;
	}
	
	public String toString() {
		String group = "";
		group += "tag = " + this.tag + ", ";
		group += "numSamples = " + this.numSamples + ", ";
		group += "numSNPs = " + this.snps.size() + ", ";
		if(this.subPopulations != null) group += "numSubPopulations = " + this.subPopulations.length;
		return group;
	}
}
