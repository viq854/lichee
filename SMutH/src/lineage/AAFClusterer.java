package lineage;

import java.util.ArrayList;
import java.util.Random;

/**
 * Simple implementation of a few clustering techniques
 * specialized for floating point data  
 * representing the allele frequency values of the samples in the application.
 *
 */
public class AAFClusterer {
	
	/** Convergence threshold [0,1] (will converge if the change from one iteration to the next
	 *  is smaller than this value) */
	private static final double CONVERGENCE_THRESHOLD = 0.00001;
	
	/** Default fuzzifier (determines the level of cluster fuzziness; 
	 * m = 2 in absence of knowledge) */
	protected static final int DEFAULT_FUZZIFIER = 2;
	
	/** Clustering algorithms */
	public enum ClusteringAlgorithms {
		KMEANS,
		FUZZYCMEANS
	}
	/** Distance measures */
	public enum DistanceMetric {
		EUCLIDEAN
	}
	
	/**
	 * Clustering dispatcher
	 * @param group - SNP group to cluster based on AAF data
	 * @param alg - algorithm to use for clustering
	 */
	public void clusterSubPopulations(SNPGroup group, ClusteringAlgorithms alg) {
		AAFClusterer clusterer = new AAFClusterer();
		int numClusters = 2;
		
		switch(alg) {
		case FUZZYCMEANS:
			clusterer.fuzzyCMeans(group.getAlleleFreqBySample(), group.getNumSNPs(), 
					group.getNumSamples(), numClusters, 
					AAFClusterer.DEFAULT_FUZZIFIER, DistanceMetric.EUCLIDEAN);
		case KMEANS:
			System.err.println("Method not implemented");
			System.exit(-1);
		default:
				
		}
		
	}
	
	// ---- Clustering Algorithms ----
	
	// ---- K-Means ----
	
	/**
	 * K-Means Clustering
	 * @param data - matrix of observations (numObs x numFeatures)
	 * @param k - number of clusters
	 */
	public void kmeans(double[][] data, int numObs, int numFeatures, int k) {}
	
	// ---- Fuzzy C-Means ----
	
	/**
	 * http://en.wikipedia.org/wiki/Fuzzy_clustering#Fuzzy_c-means_clustering
	 * Every point has a degree of belonging to each cluster
	 * The centroid is a means of all points weighted by the degree of belonging
	 * to cluster
	 * [Could improve clustering under noise]
	 * @param data - matrix of observations (numObs x numFeatures)
	 * @param c - number of clusters
	 * @param m - fuzzifier (determines the level of cluster fuzziness 
	 * (large m => fuzzier clusters, m = 1 => crisp partitions; m = 2 in absence of knowledge)
	 */
	public Cluster[] fuzzyCMeans(double[][] data, int numObs, int numFeatures, int c, int m, DistanceMetric d) {
	
		// 1. initialization
	
		// centers of each cluster
		double[][] centroids = new double[c][numFeatures];
				
		// will randomly pick the centers and compute the resulting coefficients
		for(int i = 0; i < c; i++) {
			centroids[i] = data[new Random().nextInt(numObs)];
		}
		
		// initial coefficients giving the degree of belonging to each cluster 
		double[][] coeff = computeFCMCoefficients(data, centroids, numObs, c, m, d);
		
		// 2. repeat until convergence
		double delta = Double.MAX_VALUE;
		
		while(delta > CONVERGENCE_THRESHOLD) {
			// 3. compute the centroid for each cluster
			for(int i = 0; i < c; i++) {
				double[] sum = new double[numFeatures];
				double sum_coeff = 0;
				for(int j = 0; j < numObs; j++) {
					double cf = coeff[j][i];
					sum_coeff += cf;
					for(int k = 0; k < numFeatures; k++) {
						sum[k] += cf*data[j][k];
					}
				}
				for(int k = 0; k < numFeatures; k++) {
					centroids[i][k] = sum[k]/sum_coeff;
				}
			}

			// 4. re-compute the coefficients
			double[][] coeff_new = computeFCMCoefficients(data, centroids, numObs, c, m, d);

			// 5. find the max change in coefficients
			delta = Double.MIN_VALUE;
			for(int i = 0; i < numObs; i++) {
				for(int j = 0; j < c; j++) {
					double change = Math.abs(coeff[i][j] - coeff_new[i][j]);
					if(change > delta) {
						delta = change;
					}
				}
			}
		}
		
		return getFCMHardClusters(coeff, centroids, numObs, c);
	}
	
	private Cluster[] getFCMHardClusters(double[][] coeff, double[][] centroids, int numObs, int c) {
		
		Cluster[] clusters = new Cluster[c];
		for(int i = 0; i < c; i++) {
			clusters[i] = new Cluster(centroids[i]);
		}
		for(int i = 0; i < numObs; i++) {
			// find the cluster to which this observation has the highest probability of belonging
			double prob = coeff[i][0];
			int clusterId = 0;
			for(int j = 0; j < c; j++) {
				if(coeff[i][j] > prob) {
					prob = coeff[i][j];
					clusterId = j;
				}
			}
			clusters[clusterId].addMember(i);
		}
		return clusters;
	}
	
	private double[][] computeFCMCoefficients(double[][] data, double[][] centroids, int numObs, 
			int c, int m, DistanceMetric d) {
		double[][] coeff_new = new double[numObs][c];
		for(int i = 0; i < numObs; i++) {
			double[] distToCenters = new double[c];
			for(int j = 0; j < c; j++) {
				// compute the distance to each centroid
				distToCenters[j] = getDistance(data[i], centroids[j], d);
			}
			for(int j = 0; j < c; j++) {
				double distRatioSum = 0;
				for(int k = 0; k < c; k++) {
					distRatioSum += Math.pow(distToCenters[c]/distToCenters[k], 2/(m-1));
				}
				coeff_new[i][c] = 1/distRatioSum;
			}
		}
		return coeff_new;
	}
	
	// ---- Distance Metrics ---- 
	
	private double getDistance(double[] x, double[] y, DistanceMetric d) {
		switch(d) {
		case EUCLIDEAN:
			return getEuclideanDistance(x, y);
		default:
			return 0;	
		}
	}
	
	/**
	 * Computes the Euclidean distance (2-norm distance) between 2 vectors
	 * @requires arrays x and y to have the same length
	 * @return Euclidean distance between x and y
	 */
	private double getEuclideanDistance(double[] x, double[] y) {
		double diffSum = 0;
		for(int i = 0; i < x.length; i++) {
			diffSum += Math.pow(Math.abs(x[i] - y[i]), 2);
		}
		return Math.sqrt(diffSum);
	}

	// ---- Utilities ----
	
	protected class Cluster {
		
		/** Cluster centroid */
		private double[] centroid;
		
		/** List of observations assigned to this cluster */
		private ArrayList<Integer> membership;
		
		public Cluster(double[] clusterCentroid) {
			centroid = clusterCentroid;
			membership = new ArrayList<Integer>();
		}
		
		public Cluster(double[] clusterCentroid, ArrayList<Integer> assignments) {
			centroid = clusterCentroid;
			membership = assignments;
		}
		
		public double[] getCentroid() {
			return centroid;
		}
		
		public ArrayList<Integer> getMembership() {
			return membership;
		}
		
		public void addMember(int obsId) {
			membership.add(new Integer(obsId));
		}
		
	}
	
}
