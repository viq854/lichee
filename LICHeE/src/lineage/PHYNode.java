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
import java.util.ArrayList;

import lineage.AAFClusterer.Cluster;

/**
 * Node in the phylogenetic graph
 * Represents a sub-population or sample (if leaf)
 * Is associated with a given AAF (alternative allele frequency)
 * 
 * @autor viq
 */
public class PHYNode implements Serializable, Comparable<PHYNode> {
	private static final long serialVersionUID = 1L;

	/** Sub-population cluster that the node represents */
	private Cluster cluster;		
	
	/** SNV group the node belongs to */
	protected SNVGroup snvGroup;
	
	/** Flag indicating if the node is a sample leaf*/
	protected boolean isLeaf;
	
	/** Flag indicating if the node is the germline root */
	private boolean isRoot;
	
	/** The sample id if node is a leaf*/
	private int leafSampleId;
	
	/** Debugging-only node id */
	private int nodeId;
	
	/** Level in the constraint network */
	private int level;
	
	/** 
	 * Internal node constructor
	 * @param g - SNV group the node belongs to
	 * @param nodeClusterId
	 */
	public PHYNode(SNVGroup g, int nodeClusterId, int networkLevel, int uniqueId) {
		snvGroup = g;
		cluster = snvGroup.getSubPopulations()[nodeClusterId];
		isLeaf = false;
		nodeId = uniqueId;
		level = networkLevel;
	}
	
	/**
	 * Leaf node constructor - represents each tumor sample
	 * @param sampleId - ID of the represented tumor sample
	 */
	public PHYNode(int networkLevel, int sampleId, int uniqueId) {
		isLeaf = true;
		leafSampleId = sampleId;
		nodeId = uniqueId;
		level = 0;
	}
	
	/**
	 * Root node constructor
	 */
	public PHYNode(int networkLevel, int uniqueId) {
		isRoot = true;
		nodeId = uniqueId;
		level = networkLevel;
	}
	
	/**
	 * Returns the SNV cluster that this node represents
	 */
	public Cluster getCluster() {
		return cluster;
	}
	
	/**
	 * Returns true if the node is a leaf
	 */
	public boolean isLeaf() {
		return isLeaf;
	}
	
	/**
	 * Returns true if the node is a root
	 */
	public boolean isRoot() {
		return isRoot;
	}
	
	/**
	 * Returns the ID of the sample 
	 * @requires node is a leaf
	 */
	public int getLeafSampleId() {
		return leafSampleId;
	}
	
	public int getNodeId() {
		return nodeId;
	}
	
	public int getLevel() {
		return level;
	}
	
	public SNVGroup getSNVGroup() {
		return snvGroup;
	}
	
	/**
	 * Returns the total number of samples in the network
	 */
	public int getNumSamples() {
		if(isRoot) {
			return level-1;
		}
		return snvGroup.getNumSamplesTotal();
	}
	
	/**
	 * Returns the cluster centroid AAF for the given sample id
	 * Returns 0 if the sample is not represented
	 */
	public double getAAF(int sampleId) {
		if(isRoot) {
			return Parameters.AAF_MAX;
		} 
		if(isLeaf) {
			return 0;
		}	
		int sampleIndex = snvGroup.getSampleIndex(sampleId);
		if(sampleIndex == -1) {
			return 0;
		}
		return cluster.getCentroid()[sampleIndex];
	}
	
	/**
	 * Returns the cluster standard deviation for the given sample id
	 * Returns 0 if the sample is not represented
	 */
	public double getStdDev(int sampleId) {
		if(isRoot) {
			return 0;
		} 
		if(isLeaf) {
			return 0;
		}
		
		int sampleIndex = snvGroup.getSampleIndex(sampleId);
		if(sampleIndex == -1) {
			return 0;
		}
		return cluster.getStdDev()[sampleIndex];
	}
	
	/**
	 * Returns the SNV entries in the cluster
	 * corresponding to this node
	 */
	public ArrayList<SNVEntry> getSNVs(ArrayList<SNVEntry> groupSNVs) {
		ArrayList<Integer> ids = cluster.getMembership();
		ArrayList<SNVEntry> nodeSNVs = new ArrayList<SNVEntry>();
		for(int i : ids) {
			nodeSNVs.add(groupSNVs.get(i));
		}
		return nodeSNVs;
	}
	
	public String toString() {
		String node = "Node " + nodeId + ": ";
		if(!isLeaf && !isRoot) {
			node += "group tag = " + snvGroup.getTag() + ", ";
			node += cluster.toString();
		} else if(isLeaf) {
			node += "leaf sample id = " + leafSampleId;
		} 
		return node;
	}
	
	public String getLabel() {
		String node = "";
		if(!isLeaf && !isRoot) {
			node += nodeId + ": \n";
			node += snvGroup.getTag() + "\n";
			node += "("+cluster.getMembership().size()+")";
			//node += cluster.toString();
		} else if(isLeaf) {
			node += "sample " + leafSampleId;
		} else {
			node += "germline";
		}
		return node;
	}
	
	public String getLongLabel() {
		String node = "";
		if(!isLeaf && !isRoot) {
			node += "Group: " + snvGroup.getTag() + "\n";
			node += cluster.toString();
		} else if(isLeaf) {
			node += "sample " + leafSampleId;
		} else {
			node += "germline";
		}
		return node;
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof PHYNode)) {
			return false;
		}
		PHYNode n = (PHYNode) o;
		if(this.nodeId == n.nodeId) {
			return true;
		} else {
			return false;
		}
	}
	
	public int hashCode() {
		return nodeId;
	}

	@Override
	public int compareTo(PHYNode arg0) {
		if(arg0.getNodeId() < this.nodeId) {
			return -1;
		} else if(arg0.getNodeId() > this.nodeId) {
			return 1;
		}
		return 0;
	}
	
}
