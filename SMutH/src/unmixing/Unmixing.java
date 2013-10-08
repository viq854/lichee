package unmixing;

import io.VCFConstants;
import io.VCFDatabase;
import io.VCFEntry;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Class: Unmixing
 * Constructor: (None)
 * Last Edited: February 22, 2013
 * ----
 * This class is the main program. It calls other classes to build
 * the tree from the initial input. This input should be a
 * VCF file.
 * Input: 
 *  3 samples with 4 LOHs
 *   |S1|S2|S3|
 * L1|  |  |  |
 * L2|  |  |  |
 * L3|  |  |  |
 * L4|  |  |  |
 * 
 */

public class Unmixing {
	
	private static final double MAXNLOHs = 10;
	private static final double LAF_ERROR = .05;
	public static final double Depth_ERROR = .05;

	private static final double LAF_THR = .4;
	private static final double EPSILON = 1e-10;

	private int numofLOHs;
	private int numofSamples;
	private ArrayList<CNVregion> CNVs;
	private double[][] LOHs;
	private double[][] fractions;
	private ArrayList<Integer> components;
	private boolean[] NORMALsamples; 
	
	public static void main(String[] args) {
		/**
		 * Create HG files
		 */
		
		//1,5 - 2,0 - 3,2 -4,3 -5,5 -6,5
		/*String testName = "Patient_6";//s.toString();
		VCFConstants.NormalSample = 5;//s.normal - 1;
		String path =  "/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/"+testName+"/";
		String inputFile = path+testName+".recalibrated.hardfiltered.vcf";
		String hgFile = path+testName+".HG.txt";
		VCFDatabase vcfDB = new VCFDatabase(inputFile,hgFile);
		vcfDB.generateMatrix("output.txt");
		*/
		
		//double[][] listofLohs = new double[][] {{0.6,0.5,0.2,0.43},{1,0.5,0.2,0.43},{1,0.8,1,1},{1,0.8,0.2,1}};
		//double[][] listofLohs = new double[][] {{0.6,0.5},{1,0.5},{1,0.8}};
		new Unmixing("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_2/Patient_2.BP.txt",0);
		

	}
	public Unmixing(String BPFILE, int normalSample){
		
	    // 1. build the matrix
		/*double[][] LAFlist = new double [db.getNumofSamples()][db.getHGEntries().size()];
		int[][] Depthlist = new int [db.getNumofSamples()][db.getHGEntries().size()];

		
		int counter = 0;
		for (VCFEntry entry: db.getHGEntries()){
			//entry.getChromosome()+"\t"+entry.getPosition();
			for (int i = 0; i < numofSamples; i++){
				Depthlist[i][counter] = entry.getReadDepth(i);
				LAFlist[i][counter] = (entry.getRefCount(i) < entry.getAltCount(i)? entry.getRefCount(i):entry.getAltCount(i))/(double)(Depthlist[i][counter]) ;
			}
			counter++;
		}
		*/
		// 2. smoothing the LAFlist...
		// 3. segmentations ...	-> create CNVs
		CNVs = new ArrayList<CNVregion>();
		BufferedReader rd;
		String currLine;
		try {
			rd = new BufferedReader(new FileReader(BPFILE));
			currLine = rd.readLine();
		
		int lastChr =1, lastPoint = 1;
		
		while (currLine != null){
	
			String[] parts = currLine.split("\t");
			
			int currChr= new Integer(parts[0]).intValue();
			int currPoint= new Integer(parts[1]).intValue(); //currPoint-lastPoint
			double[] avgLAF = new double [parts.length/2-1];
			double[] avgDepth = new double [parts.length/2-1];
			System.out.println("Point:"+ currChr+" "+currPoint);
			
			double normalDepth = new Double(parts[normalSample*2+3]).doubleValue();
				
				
			for (int i = 1; i <= parts.length/2-1; i++){
				avgLAF[i-1] = new Double(parts[i*2]).doubleValue();
				avgDepth[i-1] =  new Double(parts[i*2+1]).doubleValue()/normalDepth;
				System.out.printf("%.3f %.3f \t", avgLAF[i-1], avgDepth[i-1]);
			}
			System.out.print("\n");
			
			if (validCNV(avgLAF)){
				if (CNVs.size() >0 &&
					//CNVregion lastCNV = CNVs.get(CNVs.size()-1);
					 CNVs.get(CNVs.size()-1).getEndPos() == lastPoint && 
					 similarLAFs(CNVs.get(CNVs.size()-1).getAvgLAF(),avgLAF) &&
					 similarDepths(CNVs.get(CNVs.size()-1).getAvgDepth(),avgDepth)){
						CNVs.get(CNVs.size()-1).append(currChr,currPoint,avgLAF,avgDepth);
				}else{
					CNVs.add(new CNVregion(lastChr,lastPoint,currChr,currPoint,avgLAF,avgDepth));
				}
			}
			
			
			lastChr = currChr;
			lastPoint = currPoint;
			currLine = rd.readLine();
		}
		rd.close();
		} catch (IOException e) {
			System.out.println("LOH input file Reading Error!");
		}
		
		if (CNVs.size() == 0){
			System.out.println("NO CNV!");
		}
		// 4. Identify unique LOHs
		//preprocessing
		/*
		 * cluster similar LOHs
		 */		
		//error margin???
		numofSamples = CNVs.get(CNVs.size()-1).getAvgLAF().length;
		LOHs = new double [20][numofSamples];
		
		for (CNVregion cnv: CNVs){
			//System.out.println(cnv);

			if (cnv.isLoss()){
				System.out.println("LOSS" + cnv);
				double[] LAFs = cnv.getAvgLAF();
				for (int j=0; j< LAFs.length;j++)
					if (LAFs[j] >= LAF_THR) LOHs[numofLOHs][j] = 1;
					else LOHs[numofLOHs][j] = LAFs[j]/(1-LAFs[j]);
				
				int i=0;
				for (; i<numofLOHs;i++){
					//System.out.println("check similarity "+ LOHs[i].toString()+ " "+ cnv.getAvgLAF().toString());

					if (similarLAFs(LOHs[i], LOHs[numofLOHs])){
							
						//updateLOH(i,cnv);
						break;
					}
				}
				if (i == numofLOHs){
				numofLOHs++;
				}
					
			}
		}
		
		/*if (numofLOHs == 1){
			System.out.println("Normal contamination in each sample:");
			for (int s=0; s< numofSamples; s++){
				if (LOHs[0][s] == 1) System.out.print("NA ");
				else System.out.printf("%.2f ", LOHs[0][s]);

			}
			return;
		}*/
	
		System.out.print("List of LOHs:\n");
		for (int i=0; i< numofLOHs;i++){
			for (int j=0; j< LOHs[0].length;j++){
				System.out.print(LOHs[i][j]+" ");
			}
			System.out.print("\n");

		}
		
		///Fast Crazy code
		fractions = new double [numofSamples][numofLOHs+1];
		for (int j=0; j< numofSamples;j++){
			System.out.println("\n The fractions of components in sample "+j);

			ArrayList<Double> sl = new ArrayList<Double>();
			sl.add(new Double(0));
			sl.add(new Double(1));
			for (int i=0; i< numofLOHs;i++)
				if (LOHs[i][j] != 1) sl.add(new Double(LOHs[i][j]));
			Collections.sort(sl);
			for (int i=0; i< sl.size()-1;i++ ){
				if (sl.get(i+1) - sl.get(i) > LAF_ERROR ){
					System.out.printf("%.2f ", sl.get(i+1) - sl.get(i));
					fractions[j][i] = sl.get(i+1) - sl.get(i);
				}
				//if(sl.get(i+1) == 1) break;
			}
			
		
		}
		
		
		
	}
	
	
	private void decomposition(){
// compress columns?
		
		//report LOHs
		//LOHs[2][3]=LOHs[0][3]; //?????
		
		NORMALsamples = new boolean[numofSamples];
		for (int j=0; j< LOHs[0].length;j++)
			NORMALsamples[j] = true;
		

	
		// 5. decomposition 
		findComponenets();
		System.out.println("List of components,");
		for (int i=1; i < components.size(); i++){
			System.out.println(components.get(i));
		}
		
		fractions = new double [numofSamples][components.size()];
		
		System.out.println("The fractions of components each sample,");
		for (int s=0; s< numofSamples; s++){
			if (NORMALsamples[s]) continue;
			System.out.println("Sample "+s+":");
			solve(s);
		}
		
		
		/*for (int s=0; s< numofSamples; s++){
			System.out.println("The fractions of components in Sample "+s+" are:");
	
		}*/
	}
	
	
	private void findComponenets(){
		components = new ArrayList<Integer>();
		int numofpoints = (int)Math.pow(2,numofLOHs);
		components.add(numofpoints-1);
		
		for (int i=0; i < numofpoints-1; i++){
			if (isNeeded(int2binary(i))){
				System.out.println("needed: "+i);
				components.add(new Integer(i));
			}
		}
		int i = 0;
		while (components.size() < numofLOHs+1){
			if (!components.contains(new Integer(i)))
					components.add(new Integer(i));
			i++;
		}
	}
	
	private boolean[] int2binary(int input){
		 boolean[] bits = new boolean[numofLOHs];
		 for (int i = numofLOHs-1; i >= 0; i--) {
		     bits[i] = (input & (1 << i)) != 0;
		 }
		 return bits;
	}
	
	private boolean isNeeded(boolean[] vertex){
		for (int s=0; s< numofSamples; s++){
			if (NORMALsamples[s]) continue;
			boolean needed = true;
			checkpoint:
			for (int u=0; u < numofLOHs-1; u++)
			for (int v=u+1; v < numofLOHs; v++){
				if  (!vertex[u]) {
						if  (!vertex[v]){
							if (LOHs[u][s] + LOHs[v][s] > 1){
								needed =false; break checkpoint;
							}
						}else{
							if (LOHs[u][s] >= LOHs[v][s]){
								needed =false; break checkpoint;
							}
						}
				}else{
					if  (!vertex[v]){
						if (LOHs[u][s] <= LOHs[v][s] ){
							needed =false; break checkpoint;
						}
					}else{
						if (LOHs[u][s] + LOHs[v][s] < 1){
							needed =false; break checkpoint;
						}
					}
				}
			}
			if (needed){
				System.out.print("for sample "+s);
				return true;
			}
		}
		return false;
	}
	
	public  void solve(int s) 
    {
		double [][] A =  new double [numofLOHs+1][components.size()];
		double [] B = new double [numofLOHs+1];
		
		for (int i=0; i < components.size(); i++){
			boolean b[] = int2binary(components.get(i).intValue());
			for (int j=0; j < numofLOHs; j++){			
				A[j][i] = b[j]? 1:0;
				B[j] = LOHs[j][s];
			}
			A[numofLOHs][i] = 1;
			B[numofLOHs] = 1;
		}
		
		int n = components.size();
		int N = numofLOHs+1;
    

		for (int p = 0; p < n; p++) 
        {
          int max = p;
          for (int i = p + 1; i < N; i++) 
          {
              if (Math.abs(A[i][p]) > Math.abs(A[max][p]))
                  max = i;
          }
          double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
          double   t    = B[p]; B[p] = B[max]; B[max] = t;


          if (Math.abs(A[p][p]) <= EPSILON) {
              throw new RuntimeException("Matrix is singular or nearly singular");
          }

          for (int i = p + 1; i < N; i++) {
              double alpha = A[i][p] / A[p][p];
              B[i] -= alpha * B[p];
              for (int j = p; j < n; j++) {
                  A[i][j] -= alpha * A[p][j];
              }
          }
        }
   
      //double[] x = new double[n];
      for (int i = n - 1; i >= 0; i--) 
        {
          double sum = 0.0;
          for (int j = i + 1; j < n; j++) 
              sum += A[i][j] * fractions[s][j];
          fractions[s][i] = (B[i] - sum) / A[i][i];
      }
   
      for (int i = 1; i < n; i++) 
          System.out.println("C" + i + "  " +fractions[s][i]);

    }
	
	
	public double[] getFractions(int sample){
		return fractions[sample]; 
	}
	
	public double getNormalContamination(int sample){
		return fractions[sample][0]; 
	}
	
	public int minNumComponents(int sample){
		int comp = 0;
		while (fractions[sample][comp] > 0) comp++; 
		return comp;
	}
	
	public boolean validCNV(double[] a){
		for (int i =0; i< a.length; i++){
			if (a[i]  < LAF_THR)
				return true;
		}
		return false;
	}
	
	public boolean similarLAFs(double[] a, double[] b){
		for (int i =0; i< a.length; i++){
			if (Math.abs(a[i] - b[i]) > LAF_ERROR)
				return false;
		}
		return true;
	}
	
	
	public boolean similarDepths(double[] a, double[] b){
		for (int i =0; i< a.length; i++){
			if (Math.abs(a[i] - b[i]) > Depth_ERROR)
				return false;
		}
		return true;
	}
	
	public ArrayList<CNVregion> getCNVs(){
		return CNVs;
	}
}
