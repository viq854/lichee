package lineage;

import util.*;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import lineage.AAFClusterer.Cluster;
import lineage.AAFClusterer.DistanceMetric;

/**
 * An SNV group is a set of SNVs occurring in a given subset of samples.
 * All SNVs are partitioned into SNV groups based on their occurrence across samples.
 * Given S samples, there can be at most 2^S different groups.
 * An SNV group is uniquely identified by an S-bit binary tag (each bit corresponding to
 * a given sample), where a bit is set if that sample contains the SNVs in this group.
 *
 */
public class SNVGroup implements Serializable {

	private static final long serialVersionUID = 1L;

	/** Binary tag identifying the group 
	 * (the length of the tag is equal to the number of input samples) */
	private String tag;
	
	/** Number of samples represented by this group */
	private int numSamples;
	
	/** Indices of the samples represented in this group (from 0 to |tag|-1 MSF order) */
	private int[] sampleIndex;
	
	/** Alternative allele frequency data matrix (numSNVs x numSamples) */
	private transient double[][] alleleFreqBySample;
	
	/** SubPopulation clusters */
	private Cluster[] subPopulations;
	
	/** SNVs assigned to this group */
	private transient ArrayList<SNVEntry> snvs;
	
	/** Number of solid/robust mutations in the group */
	private transient int numRobustSNVs;
	
	/** Flag indicating whether this group is robust */
	private transient boolean isRobust;
	
	private static Logger logger = LineageEngine.logger;
	
	public SNVGroup(String groupTag, ArrayList<SNVEntry> groupSNVs, int groupNumRobustSNVs, boolean isGroupRobust) {
		tag = groupTag;
		numRobustSNVs = groupNumRobustSNVs;
		isRobust = isGroupRobust;
		numSamples = 0;		
		sampleIndex = new int[tag.length()];
		for(int i = 0; i < tag.length(); i++) {
			if(tag.charAt(i) == '1') {
				sampleIndex[numSamples] = i;
				numSamples++;
			}
		}
		snvs = groupSNVs;
		alleleFreqBySample = new double[snvs.size()][numSamples];
		for(int i = 0; i < snvs.size(); i++) {
			SNVEntry snv = snvs.get(i);
			for(int j = 0; j < numSamples; j++) {
				alleleFreqBySample[i][j] = snv.getAAF(sampleIndex[j]);
			}
		}
		
		//logger.addHandler(LineageEngine.logger.getHandlers()[0]);
		//logger.setLevel(LineageEngine.logger.getLevel());
		//logger.log(Level.FINE, "Created group: " + this.toString());
	}
	
	// Getters/Setters
	
	public ArrayList<SNVEntry> getSNVs() {
		return snvs;
	}
	
	public double[][] getAlleleFreqBySample() {
		return alleleFreqBySample;
	}
	
	public int getNumSamples() {
		return numSamples;
	}
	
	public int getNumSamplesTotal() {
		return tag.length();
	}
	
	public int getNumSNVs() {
		return snvs.size();
	}
	
	public Cluster[] getSubPopulations() {
		return subPopulations;
	}
	
	public String getTag() {
		return tag;
	}
	
	public int getNumRobustSNVs() {
		return numRobustSNVs;
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
		if(!(o instanceof SNVGroup)) {
			return false;
		}
		SNVGroup g = (SNVGroup) o;
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
		
		// 1. filter out clusters that are too small (unless there is only one cluster)
		ArrayList<Cluster> filteredClusters = new ArrayList<Cluster>();
		for(Cluster c : clusters) {
			if((clusters.length == 1) || (c.getMembership().size() >= Parameters.MIN_CLUSTER_SIZE)) {
				filteredClusters.add(c);
			}
		}
		
		if(filteredClusters.size() < 1) {
			logger.log(Level.WARNING, "All subpopulation clusters in group " + tag + " have been filtered out. "
					+ "The group is now empty.");
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
				((minDistQueue.get(0).distance < Parameters.MAX_COLLAPSE_CLUSTER_DIFF))) {// || (numClusters > Parameters.MAX_CLUSTER_NUM))) {
			ClusterPairDistance pd = minDistQueue.remove(0);
			Cluster c1 = clusters[pd.clusterId1];
			Cluster c2 = clusters[pd.clusterId2];
			
			// collapse into c1
			for(Integer obs : c2.getMembership()) {
				c1.addMember(obs);
			}
			numClusters--;
			
			c1.recomputeCentroid(alleleFreqBySample, snvs.size(), numSamples);
			filteredClusters.remove(c2);
			logger.log(Level.FINE, "Collapse clusters: group = " + tag + " cluster " + pd.clusterId1 + " and " + pd.clusterId2 + 
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
		group += "(tag = " + this.tag + ", ";
		group += "numSamples = " + this.numSamples + ", ";
		group += "numSNVs = " + this.snvs.size() + ") ";
		if(this.subPopulations != null) group += "numSubPopulations = " + this.subPopulations.length;
		return group;
	}
}
