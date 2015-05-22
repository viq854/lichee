/*
 * Program LICHeE for multi-sample cancer phylogeny reconstruction
 * by Victoria Popic (viq@stanford.edu) 2014
 *
 * MIT License
 *
 * Copyright (c) 2014 Victoria Popic.
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/


package lineage;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

import weka.clusterers.ClusterEvaluation;
import weka.clusterers.EM;
import weka.clusterers.SimpleKMeans;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Simple implementation of a few clustering techniques
 * for allele frequency values of sample SNVs.
 * 
 * @autor viq
 */
public class AAFClusterer implements Serializable {
	
	private static final long serialVersionUID = 1L;

	/** Convergence threshold [0,1] (will converge if the change from one iteration to the next
	 *  is smaller than this value) */
	private static final double CONVERGENCE_THRESHOLD = 0.00001;
	
	/** Default fuzzifier (determines the level of cluster fuzziness; 
	 * m = 2 in absence of knowledge) */
	protected static final int DEFAULT_FUZZIFIER = 2;
	
	/** Clustering algorithms */
	public enum ClusteringAlgorithms {
		KMEANS,
		FUZZYCMEANS,
		EM
	}
	/** Distance measures */
	public enum DistanceMetric {
		EUCLIDEAN,
		AVG_PER_SAMPLE
	}
	
	/**
	 * Clustering dispatcher
	 * @requires the number of SNVs in a group to be bigger than 1
	 * @param group - SNV group to cluster based on AAF data
	 * @param alg - algorithm to use for clustering
	 */
	public Cluster[] clusterSubPopulations(SNVGroup group, ClusteringAlgorithms alg, int minNumClusters) {		
		switch(alg) {
		case FUZZYCMEANS:
			return fuzzyCMeans(group.getAlleleFreqBySample(), group.getNumSNVs(), 
					group.getNumSamples(), minNumClusters, 
					AAFClusterer.DEFAULT_FUZZIFIER, DistanceMetric.EUCLIDEAN);
		case KMEANS:
			return kmeans(group.getAlleleFreqBySample(), group.getNumSNVs(), group.getNumSamples(), minNumClusters);
		case EM:
			return em(group.getAlleleFreqBySample(), group.getNumSNVs(), group.getNumSamples());
		default:
			return null;	
		}
	}
	
	// ---- Clustering Algorithms ----
	
	// ---- K-Means ----
	
	/**
	 * K-Means Clustering
	 * @param data - matrix of observations (numObs x numFeatures)
	 * @param k - number of clusters
	 */
	public Cluster[] kmeans(double[][] data, int numObs, int numFeatures, int k) {
		Instances ds = convertMatrixToWeka(data, numObs, numFeatures);
		
		// uses Euclidean distance by default
		SimpleKMeans clusterer = new SimpleKMeans();
		try {
			clusterer.setPreserveInstancesOrder(true);
			clusterer.setNumClusters(k);
			clusterer.buildClusterer(ds);
			
			// cluster centers
			Instances centers = clusterer.getClusterCentroids();
			Cluster[] clusters = new Cluster[centers.numInstances()];
			for(int i = 0; i < centers.numInstances(); i++) {
				Instance inst = centers.instance(i);
				double[] mean = new double[inst.numAttributes()];
				for(int j = 0; j < mean.length; j++) {
					mean[j] = inst.value(j);
				}
				clusters[i] = new Cluster(mean, i);
			}
			
			// cluster members
			int[] assignments = clusterer.getAssignments();
			for(int i = 0; i < assignments.length; i++) {
				clusters[assignments[i]].addMember(i);
			}
			return clusters;
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
			return null;
		}
		
	}
	
	/**
	 * Expectation Maximization clustering
	 * @param data - matrix of observations (numObs x numFeatures)
	 * @param k - number of clusters
	 */
	public Cluster[] em(double[][] data, int numObs, int numFeatures) {
		Instances ds = convertMatrixToWeka(data, numObs, numFeatures);
		EM clusterer = new EM();
		try {
			clusterer.buildClusterer(ds);
			ClusterEvaluation eval = new ClusterEvaluation();                                        
			eval.setClusterer(clusterer);                                  
			eval.evaluateClusterer(new Instances(ds));                               
			int numClusters = eval.getNumClusters();
			
			Cluster[] clusters = new Cluster[numClusters];
			double[][] clusterCentroids = new double[numClusters][numFeatures];
			int[] clusterCount = new int[numClusters];
			
			double[] assignments = eval.getClusterAssignments();
			for(int i = 0; i < ds.numInstances(); i++) {
				Instance inst = ds.instance(i);
				int clusterId = (int) assignments[i];
				for(int j = 0; j < numFeatures; j++) {
					clusterCentroids[clusterId][j] += inst.value(j);
				}
				clusterCount[clusterId]++;
			}
			
			for(int i = 0; i < numClusters; i++) {
				double[] mean = new double[numFeatures];
				for(int j = 0; j < numFeatures; j++) {
					mean[j] = clusterCentroids[i][j]/clusterCount[i];
				}
				clusters[i] = new Cluster(mean, i);
			}
			
			// cluster members & std dev
			double[][] clusterStdDev = new double[numClusters][numFeatures];
			for(int i = 0; i < ds.numInstances(); i++) {
				int clusterId = (int)assignments[i];
				clusters[clusterId].addMember(i);
				for(int j = 0; j < numFeatures; j++) {
					clusterStdDev[clusterId][j] += Math.pow(ds.instance(i).value(j) - clusters[clusterId].getCentroid()[j], 2);
				}
			}
			
			for(int i = 0; i < numClusters; i++) {
				double[] dev = new double[numFeatures];
				for(int j = 0; j < numFeatures; j++) {
					dev[j] = Math.sqrt(clusterStdDev[i][j]/clusterCount[i]);
				}
				clusters[i].setStdDev(dev);
			}
			
			return clusters;
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
			return null;
		}
	}
	
	// ---- Fuzzy C-Means ----
	
	/**
	 * http://en.wikipedia.org/wiki/Fuzzy_clustering#Fuzzy_c-means_clustering
	 * Every point has a degree of belonging to each cluster
	 * The centroid is a means of all points weighted by the degree of belonging to cluster
	 * @param data - matrix of observations (numObs x numFeatures)
	 * @param c - number of clusters
	 * @param m - fuzzifier (determines the level of cluster fuzziness 
	 * (large m => fuzzier clusters, m = 1 => crisp partitions; m = 2 in absence of knowledge)
	 */
	public Cluster[] fuzzyCMeans(double[][] data, int numObs, int numFeatures, int c, int m, DistanceMetric d) {
	
		// 1. initialization
	
		// centers of each cluster
		double[][] centroids = new double[c][numFeatures];
				
		// pick random coefficients
		double[][] coeff = new double[numObs][c];
		for(int i = 0; i < numObs; i++) {
			double totalProb = 0;
			for(int j = 0; j < c-1; j++) {
				coeff[i][j] = 0.1*new Random().nextInt(10 - (int)totalProb*10);
				totalProb += coeff[i][j];
			}
			coeff[i][c-1] = 1 - totalProb;
		}
		
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
			coeff = coeff_new;
		}
		
		return getFCMHardClusters(coeff, centroids, numObs, c);
	}
	
	private Cluster[] getFCMHardClusters(double[][] coeff, double[][] centroids, int numObs, int c) {
		
		Cluster[] clusters = new Cluster[c];
		for(int i = 0; i < c; i++) {
			clusters[i] = new Cluster(centroids[i], i);
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
					distRatioSum += Math.pow(distToCenters[j]/distToCenters[k], 2/(m-1));
				}

				coeff_new[i][j] = 1/distRatioSum;
			}
		}
		return coeff_new;
	}
	
	// ---- Distance Metrics ---- 
	
	private double getDistance(double[] x, double[] y, DistanceMetric d) {
		switch(d) {
		case EUCLIDEAN:
			return getEuclideanDistance(x, y);
		case AVG_PER_SAMPLE:
			return getAvgSampleDistance(x, y);
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
	
	private double getAvgSampleDistance(double[] x, double[] y) {
		double diffSum = 0;
		for(int i = 0; i < x.length; i++) {
			diffSum += Math.abs(x[i] - y[i]);
		}
		return diffSum/x.length;
	}
	
	// ---- WEKA-related Utilities ----
	public Instances convertMatrixToWeka(double[][] data, int numObs, int numFeatures) {
		// convert the data to WEKA format
		FastVector atts = new FastVector();
		for(int i = 0; i < numFeatures; i++) {
			atts.addElement(new Attribute("Feature" + i, i));
		}
		
		Instances ds = new Instances("AAF Data", atts, numObs);
		for(int i = 0; i < numObs; i++) {
			ds.add(new Instance(numFeatures));
		}
		for(int i = 0; i < numFeatures; i++) {
			for(int j = 0; j < numObs; j++) {
				ds.instance(j).setValue(i, data[j][i]);
			}
		}
		return ds;
	}
	

	/**
	 * Cluster of observation points
	 * Each cluster has an associated centroid point and a list of members
	 */
	protected class Cluster implements Serializable {
		
		private static final long serialVersionUID = 1L;

		/** Cluster id, unique per group */
		private int id;
		
		/** Cluster centroid */
		private double[] centroid;
		
		/** Cluster standard deviation */
		private double[] stdDev = null;
		
		/** List of observations assigned to this cluster */
		private ArrayList<Integer> members;
		
		private boolean robust = false;
		
		public Cluster(double[] clusterCentroid, int clusterId) {
			centroid = clusterCentroid;
			members = new ArrayList<Integer>();
			id = clusterId;
		}
		
		public Cluster(double[] clusterCentroid, ArrayList<Integer> assignments, int clusterId) {
			centroid = clusterCentroid;
			members = assignments;
			id = clusterId;
			
		}
		
		/**
		 * Compute the distance of a given observation to this cluster
		 * Currently the method computes the distance to the centroid
		 * @param x - observation point
		 * @param d - distance function to be used (e.g. Euclidean)
		 * @return distance of observation to cluster
		 */
		public double getDistanceToCluster(double[] x, DistanceMetric d) {
			return getDistance(x, centroid, d);
		}
		
		/**
		 * Add a new observation to the cluster
		 * @param obsId - Id of the observation (index in the data matrix)
		 */
		public void addMember(int obsId) {
			members.add(new Integer(obsId));
		}
		
		/** 
		 * Computes the mean of all the members of the cluster
		 */
		public void recomputeCentroidAndStdDev(double[][] data, int numObs, int numFeatures) {
			double[] newCentroid = new double[numFeatures];
			for(int i = 0; i < members.size(); i++) {
				for(int j = 0; j < numFeatures; j++) {
					newCentroid[j] += data[members.get(i)][j];
				}
			}
			for(int j = 0; j < numFeatures; j++) {
				newCentroid[j] = newCentroid[j]/members.size();
			}
			centroid = newCentroid;
			
			double[] clusterStdDev = new double[numFeatures];
			for(int i = 0; i < members.size(); i++) {
				for(int j = 0; j < numFeatures; j++) {
					clusterStdDev[j] += Math.pow(data[members.get(i)][j] - centroid[j], 2);
				}
			}
			for(int j = 0; j < numFeatures; j++) {
				clusterStdDev[j] =  Math.sqrt(clusterStdDev[j]/members.size());
			}
			stdDev = clusterStdDev;
		}
		
		/**
		 * Returns the cluster centroid (mean) per sample
		 */
		public double[] getCentroid() {
			return centroid;
		}
		
		/**
		 * Returns the standard deviation per sample
		 * @requires setStdDev() method to have been called (currently implemented for EM only),
		 * will return null otherwise
		 */
		public double[] getStdDev() {
			return stdDev;
		}
		
		public void setStdDev(double[] dev) {
			stdDev = dev;
		}
		
		public ArrayList<Integer> getMembership() {
			return members;
		}
		
		public int getId() {
			return id;
		}
		
		public boolean isRobust() {
			return robust;
		}
		
		public void setRobust() {
			robust = true;
		}
		
		public String toString() {
			String c = "";
			c += "Size: " + members.size() + "\n";
			DecimalFormat df = new DecimalFormat("#.##");
			c += "VAF Mean: [";
			for(int i = 0; i < centroid.length; i++) {
				c += " " + df.format(centroid[i]) + " ";
			}
			c += "] \n";
			c += "       Stdev:";
			if(stdDev != null) {
				c += " [";
				for(int i = 0; i < stdDev.length; i++) {
					c += " " + df.format(stdDev[i]) + " ";
				}
				c += "]";
			}
			return c;
		}
	}
	
}
