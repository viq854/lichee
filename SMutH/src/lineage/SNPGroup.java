package lineage;

import io.*;

import java.util.ArrayList;

import lineage.AAFClusterer.Cluster;
import lineage.AAFClusterer.DistanceMetric;

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
	
	/** Indices of the samples represented in this group (from 0 to |tag|-1 MSF order) */
	private int[] sampleIndex;
	
	/** Alternative allele frequency data matrix (numSNPs x numSamples) */
	private double[][] alleleFreqBySample;
	
	/** SubPopulation clusters */
	private Cluster[] subPopulations;
	
	/** SNPs assigned to this group */
	private ArrayList<SNVEntry> snps;
	
	/** Number of solid/robust mutations in the group */
	private int numRobustSNPs;
	
	/** Flag indicating whether this group is robust */
	private boolean isRobust;
	
	public SNPGroup(String groupTag, ArrayList<SNVEntry> groupSNPs, int groupNumRobustSNPs, boolean isGroupRobust) {
		tag = groupTag;
		numRobustSNPs = groupNumRobustSNPs;
		isRobust = isGroupRobust;
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
	 * TODO
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
	
	public int getNumRobustSNPs() {
		return numRobustSNPs;
	}
	
	public boolean isRobust() {
		return isRobust;
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
	
	/**
	 * Returns true if the given sample contains the mutations of this group
	 */
	public boolean containsSample(int sampleId) {
		return (getSampleIndex(sampleId) != -1);
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof SNPGroup)) {
			return false;
		}
		SNPGroup g = (SNPGroup) o;
		if(this.tag == g.tag) {
			return true;
		} else {
			return false;
		}
	}
	
	// --- Sub-population Cluster Filtering / Collapse ---
	
	/** Entry in the cluster centroid distance minimum priority queue */
	protected class ClusterPairDistance {
		/** Cluster pair */
		protected int clusterId1;
		protected int clusterId2;
		/** Distance between cluster centroids */
		protected double distance;
		
		public ClusterPairDistance(int cluster1, int cluster2, double clusterDistance) {
			clusterId1 = cluster1;
			clusterId2 = cluster2;
			distance = clusterDistance;
		}
	}
	
	
	/**
	 * Set the sub-populations of this group based on clustering results
	 * Performs filtering based on cluster size, as well as collapses clusters
	 * with similar centroids
	 * @param clusters - clustering algorithm results
	 */
	public void setSubPopulations(Cluster[] clusters) {
		
		// 1. filter out clusters that are too small
		ArrayList<Cluster> filteredClusters = new ArrayList<Cluster>();
		for(Cluster c : clusters) {
			if(c.getMembership().size() >= Parameters.MIN_CLUSTER_SIZE) {
				filteredClusters.add(c);
			}
		}
		
		if(filteredClusters.size() < 1) {
			System.out.println("Warning: All clusters in group " + tag + " have been filtered out");
			subPopulations = new Cluster[filteredClusters.size()];
			subPopulations = filteredClusters.toArray(subPopulations);
			return;
		}
		
		// 2. collapse clusters that have similar centroids
		
		// compute the distance matrix between clusters
		// as long as there are clusters to collapse (i.e. cluster centroid distance
		// is less than MAX_COLLAPSE_CLUSTER_DIFF), collapse clusters with smallest distance first
		
		ArrayList<ClusterPairDistance> minDistQueue = new ArrayList<ClusterPairDistance>();
		int numClusters = filteredClusters.size();
		for(int i = 0; i < numClusters; i++) {
			for(int j = i+1; j < numClusters; j++) {
				Cluster c1 = filteredClusters.get(i);
				Cluster c2 = filteredClusters.get(j);
				double dist = c1.getDistanceToCluster(c2.getCentroid(), DistanceMetric.EUCLIDEAN);
				ClusterPairDistance pd = new ClusterPairDistance(c1.getId(), c2.getId(), dist);
				
				int k = 0;
				for(k = 0; k < minDistQueue.size(); k++) {
					if(minDistQueue.get(k).distance > dist) {
						break;
					}
				}
				minDistQueue.add(k, pd);
			}
		}
		
		while((minDistQueue.size() > 0) && 
				((minDistQueue.get(0).distance < Parameters.MAX_COLLAPSE_CLUSTER_DIFF) || (numClusters > Parameters.MAX_CLUSTER_NUM))) {
			ClusterPairDistance pd = minDistQueue.remove(0);
			Cluster c1 = clusters[pd.clusterId1];
			Cluster c2 = clusters[pd.clusterId2];
			
			// collapse into c1
			for(Integer obs : c2.getMembership()) {
				c1.addMember(obs);
			}
			numClusters--;
			
			c1.recomputeCentroid(alleleFreqBySample, snps.size(), numSamples);
			filteredClusters.remove(c2);
			System.out.println("Collapse clusters: group = " + tag + " cluster " + pd.clusterId1 + " and " + pd.clusterId2 + 
					" distance = " + pd.distance);
			
			// remove distances from c1 and c2 from the queue
			ArrayList<ClusterPairDistance> toRemove = new ArrayList<ClusterPairDistance>();
			for(ClusterPairDistance cpd : minDistQueue) {
				if(cpd.clusterId1 == c1.getId() || cpd.clusterId2 == c1.getId() 
				   || cpd.clusterId1 == c2.getId() || cpd.clusterId2 == c2.getId()) {
					toRemove.add(cpd);
				}
			}
			minDistQueue.removeAll(toRemove);
			
			// compute the distance from c1 to all the other clusters
			for(Cluster c : filteredClusters) {
				if(c.getId() == c1.getId() || c.getId() == c2.getId()) {
					continue;
				}
				double dist = c1.getDistanceToCluster(c.getCentroid(), DistanceMetric.EUCLIDEAN);
				ClusterPairDistance cpd = new ClusterPairDistance(c1.getId(), c.getId(), dist);
				
				int k = 0;
				for(k = 0; k < minDistQueue.size(); k++) {
					if(minDistQueue.get(k).distance > dist) {
						break;
					}
				}
				minDistQueue.add(k, cpd);
			}
		}
		
		subPopulations = new Cluster[filteredClusters.size()];
		subPopulations = filteredClusters.toArray(subPopulations);
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
