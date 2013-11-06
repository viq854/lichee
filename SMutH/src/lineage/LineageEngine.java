package lineage;

import java.util.ArrayList;
import java.util.HashMap;

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
		SNVDatabase db = new SNVDatabase(path+sampleName+".vcf", normalSample);
		//SNVDatabase db = new SNVDatabase(path+sampleName+".validation.txt", normalSample);
		
		
		// 2. handle normal cell contamination, CNVs, 
		//    determine the minimum number of clusters using LOH
		//    (+ any additional filtering, pre-processing)
		
		//Unmixing um = new Unmixing(path+sampleName+".BP.txt",normalSample);

		
		// 3. get the SNPs partitioned by group tag
		HashMap<String, ArrayList<SNVEntry>> snpsByTag = db.generateFilteredTAG2SNVsMap(null);
		// create the appropriate SNP group objects
		ArrayList<SNPGroup> groups = new ArrayList<SNPGroup>();
		for(String groupTag : snpsByTag.keySet()) {
			SNPGroup group = new SNPGroup(groupTag, snpsByTag.get(groupTag), db.getNumRobustSNVs(groupTag), db.isRobust(groupTag));
			groups.add(group);
		}
		
		
		// 4. cluster SNPs in each group
		AAFClusterer clusterer = new AAFClusterer();
		for(SNPGroup group : groups) {
			System.out.println("Clustering results group: " + group.getTag());
			Cluster[] clusters = clusterer.clusterSubPopulations(group, ClusteringAlgorithms.EM, 1);
			for(Cluster c : clusters) {
				System.out.println(c.toString());
			}
			group.setSubPopulations(clusters);
		}
		
		// 5. construct constraint network
		int[] sampleMutationMask = new int[db.getNumofSamples()];
		sampleMutationMask[0] = 1; // sample 0 has no mutations
		PHYGraph constrNetwork = new PHYGraph(groups, db.getNumofSamples(), sampleMutationMask);
		System.out.println(constrNetwork.toString());
		
		// 6. find all the spanning trees
		ArrayList<Tree> spanningTrees = constrNetwork.getLineageTrees();  
		System.out.println("Spanning trees: " + spanningTrees.size() + " trees total");
		
		// 7. apply AAF constraints and other filters to prune out trees
		//constrNetwork.testSpanningTrees();
		constrNetwork.filterSpanningTrees();
		System.out.println("Spanning trees: " + spanningTrees.size() + " trees after AAF constraint *filtering*");	

		if(spanningTrees.size() < Parameters.SPANNING_TREE_MIN_NUMBER) {
			System.out.println("Fixing the network");	
			// if the number of filtered spanning trees is less than acceptable,
			// fix the network (e.g. remove group nodes that are not robust)
			constrNetwork = constrNetwork.fixNetwork();
			spanningTrees = constrNetwork.getLineageTrees();  
			System.out.println("Spanning trees: " + spanningTrees.size() + " trees total");
			constrNetwork.filterSpanningTrees();
			System.out.println("Spanning trees: " + spanningTrees.size() + " trees after AAF constraint *filtering*");	
		}
		if(spanningTrees.size() > 0) {
			System.out.println(spanningTrees.get(0));
		}
		
		// 8. evaluate/rank the trees
		
		// 9. result visualization
		constrNetwork.displayNetwork();
		if(spanningTrees.size() > 0) {
			spanningTrees.get(0).displayTree();
			for(int i = 0; i < db.getNumofSamples(); i++) {
				System.out.println(spanningTrees.get(0).getLineage(i));
			}
		}
		
	}
	
	public static void main(String[] args) {
		buildLineage("/Users/viq/smuth/SMutH/data/","patient1", 0);
		//buildLineage("/Users/viq/smuth/SMutH/data/","Patient_1", 0);
		//buildLineage("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_2/","Patient_2", 0);

	}
	
}
