package lineage;

import java.util.ArrayList;
import java.util.HashMap;

import unmixing.Unmixing;

import lineage.AAFClusterer.Cluster;
import lineage.AAFClusterer.ClusteringAlgorithms;
import io.VCFDatabase;
import io.VCFEntry;

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
		VCFDatabase db = new VCFDatabase(path+sampleName+".raw.vcf", normalSample);
		
		
		// 2. handle normal cell contamination, CNVs, 
		//    determine the minimum number of clusters using LOH
		//    (+ any additional filtering, pre-processing)
		
		//Unmixing um = new Unmixing(path+sampleName+".BP.txt",normalSample);

		
		// 3. get the SNPs partitioned by group tag
		HashMap<String, ArrayList<VCFEntry>> snpsByTag = db.generateFilteredTAG2SNVsMap(null);
		// create the appropriate SNP group objects
		ArrayList<SNPGroup> groups = new ArrayList<SNPGroup>();
		for(String groupTag : snpsByTag.keySet()) {
			SNPGroup group = new SNPGroup(groupTag, snpsByTag.get(groupTag));
			groups.add(group);
		}
		
		
		// 4. cluster SNPs in each group
		AAFClusterer clusterer = new AAFClusterer();
		for(SNPGroup group : groups) {
			clusterer.clusterSubPopulations(group, ClusteringAlgorithms.FUZZYCMEANS, 2);
			System.out.println("FUZZYCMEANS Clustering results group: " + group.getTag());
			for(Cluster c : group.getSubPopulations()) {
				System.out.println(c.toString());
			}
			
			clusterer.clusterSubPopulations(group, ClusteringAlgorithms.KMEANS, 2);
			System.out.println("KMEANS Clustering results group: " + group.getTag());
			for(Cluster c : group.getSubPopulations()) {
				System.out.println(c.toString());
			}
		}
		
		// 5. incorporate CNVs
		
		// 6. construct constraint network
		//PHYGraph constrNetwork = new PHYGraph(groups, db.getNumofSamples());
		
		// 7. build and output final cell lineage
		//constrNetwork.getLineageTrees();  
		
		// 8. result visualization
		
	}
	
	public static void main(String[] args) {
		buildLineage("/Users/rahelehs/Work/cancerTree/simulation_vcfs/","tree_3_03", 0);
		//buildLineage("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_2/","Patient_2", 0);

	}
	
}
