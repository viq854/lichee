package lineage;

import java.util.ArrayList;
import java.util.HashMap;

import SMutH.TreeVisualizer;
import unmixing.Unmixing;
import lineage.AAFClusterer.Cluster;
import lineage.AAFClusterer.ClusteringAlgorithms;
import lineage.PHYGraph.Tree;
import io.*;

/**
 * Main cell lineage builder pipeline
 *
 */
public class LineageEngine {

	/**
	 * The main pipeline for constructing the cell lineage 
	 * @param vcfFileName - path to the input VCF file
	 */
	public static void buildLineage(String path, String sampleName, int normalSample) {
		// 1. load VCF data
		SNVDatabase db = new SNVDatabase(path+sampleName+".validation.txt", normalSample);
		
		
		// 2. handle normal cell contamination, CNVs, 
		//    determine the minimum number of clusters using LOH
		//    (+ any additional filtering, pre-processing)
		
		//Unmixing um = new Unmixing(path+sampleName+".BP.txt",normalSample);

		
		// 3. get the SNPs partitioned by group tag
		HashMap<String, ArrayList<SNVEntry>> snpsByTag = db.generateFilteredTAG2SNVsMap(null);
		// create the appropriate SNP group objects
		ArrayList<SNPGroup> groups = new ArrayList<SNPGroup>();
		for(String groupTag : snpsByTag.keySet()) {
			SNPGroup group = new SNPGroup(groupTag, snpsByTag.get(groupTag));
			groups.add(group);
		}
		
		
		// 4. cluster SNPs in each group
		AAFClusterer clusterer = new AAFClusterer();
		for(SNPGroup group : groups) {
			/*clusterer.clusterSubPopulations(group, ClusteringAlgorithms.FUZZYCMEANS, 2);
			System.out.println("FUZZYCMEANS Clustering results group: " + group.getTag());
			for(Cluster c : group.getSubPopulations()) {
				System.out.println(c.toString());
			}
			
			clusterer.clusterSubPopulations(group, ClusteringAlgorithms.KMEANS, 2);
			System.out.println("KMEANS Clustering results group: " + group.getTag());
			for(Cluster c : group.getSubPopulations()) {
				System.out.println(c.toString());
			}*/
			
			System.out.println("EM Clustering results group: " + group.getTag());
			clusterer.clusterSubPopulations(group, ClusteringAlgorithms.EM, 2);
			for(Cluster c : group.getSubPopulations()) {
				System.out.println(c.toString());
			}
		}
		
		// 5. incorporate CNVs
		
		// 6. construct constraint network
		int[] sampleMutationMask = new int[db.getNumofSamples()];
		sampleMutationMask[0] = 1; // sample 0 has no mutations
		PHYGraph constrNetwork = new PHYGraph(groups, db.getNumofSamples(), sampleMutationMask);
		System.out.println(constrNetwork.toString());
		
		// 7. find all the spanning trees
		ArrayList<Tree> spanningTrees = constrNetwork.getLineageTrees();  
		System.out.println("Found " + spanningTrees.size() + " trees");
		if(spanningTrees.size() > 0) {
			System.out.println(spanningTrees.get(0));
		}
		
		// 8. apply AAF constraints and other filters to prune out trees
		constrNetwork.testSpanningTrees();
		constrNetwork.filterSpanningTrees();
		System.out.println("Filtered " + spanningTrees.size() + " trees");	
		if(spanningTrees.size() > 0) {
			System.out.println(spanningTrees.get(0));
		}
		
		// 9. result visualization
		constrNetwork.displayNetwork();
		if(spanningTrees.size() > 0) {
			spanningTrees.get(0).displayTree();;
		}
		
		
	}
	
	public static void main(String[] args) {
		//buildLineage("/Users/viq/smuth/SMutH/data/","patient1", 0);
		buildLineage("/Users/viq/smuth/SMutH/data/","Patient_3", 0);
		//buildLineage("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_2/","Patient_2", 0);

	}
	
}
