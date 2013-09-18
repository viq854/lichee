package unmixing;

import io.VCFDatabase;

import java.util.ArrayList;

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
	private static final double EPSILON = 1e-10;
	
	private int numofLOHs;
	private int numofSamples;
	//private ArrayList<LOHEntry> LOHs;
	private double[][] LOHs;
	private double[][] fractions;
	ArrayList<boolean[]> components;
	
	public static void main(String[] args) {
		
		double[][] listofLohs = new double[][] {{0.6,0.5,0.2,0.43},{1,0.5,0.2,0.43},{1,0.8,1,1},{1,0.8,0.2,1}};
		//double[][] listofLohs = new double[][] {{0.6,0.5},{1,0.5},{1,0.8}};

		new Unmixing(listofLohs);
		
	}
	public Unmixing(VCFDatabase db){
	//TODO	
	}
	
	public Unmixing(double [][] listofLOHs){
		LOHs = (double[][])listofLOHs.clone();
		numofLOHs = listofLOHs.length;
		numofSamples = listofLOHs[0].length;
		
		//preprocessing
		/*
		 * cluster similar LOHs
		 */
		
		//error
		
		findComponenets();
		
		fractions = new double [numofSamples][components.size()];
		
		for (int s=0; s< numofSamples; s++){
			System.out.println("The fractions of components in Sample "+s+" are:");
			solve(s);
		}
	}
	
	private void findComponenets(){
		components = new ArrayList<boolean[]>();
		int numofpoints = (int)Math.pow(2,numofLOHs);
		components.add(int2binary(numofpoints-1));
		
		for (int i=0; i < numofpoints-1; i++){
			if (isNeeded(int2binary(i))){
				System.out.println("needed: "+i);
				components.add(int2binary(i));
			}
		}
		int i = 0;
		while (components.size() < numofLOHs+1){
			if (!components.contains(int2binary(i)))
					components.add(int2binary(i));
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
							if (LOHs[u][s] > LOHs[v][s]){
								needed =false; break checkpoint;
							}
						}
				}else{
					if  (!vertex[v]){
						if (LOHs[u][s] <= LOHs[v][s] ){
							needed =false; break checkpoint;
						}
					}else{
						if (LOHs[u][s] + LOHs[v][s] <= 1){
							needed =false; break checkpoint;
						}
					}
				}
			}
			if (needed) return true;
		}
		return false;
	}
	
	public  void solve(int s) 
    {
		double [][] A =  new double [numofLOHs+1][components.size()];
		double [] B = new double [numofLOHs+1];
		
		for (int i=0; i < components.size(); i++){
			boolean b[] = components.get(i);
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
   
      for (int i = 0; i < n; i++) 
          System.out.println("C" + i + "  " +fractions[s][i]);

    }
	
	public double getNormalContamination(int sample){
		return fractions[sample][0]; //components.get(sample);
	}
	
	public int minNumComponents(int sample){
		int comp = 0;
		for (int i = 1; i < components.size(); i++) 
		if (fractions[sample][i] > 0) comp++; //components.get(sample);
		return comp;
	}
}
