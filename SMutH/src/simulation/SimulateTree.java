package simulation;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

public class SimulateTree {
	
	private int numSamples;
	private PrintWriter pw;
	private int numGroups;
	private HashMap<Integer, ArrayList<int[]>> groups;
	
	private final static int MinSize = 5;
	private final static int MaxSize = 200;
	
	private final static int MIN_DIST_AAF = 5;
	private final static double PROB_CLUSTER = .1;

	
	public SimulateTree(int n, int m){
		numSamples = n;
		int rootTag = (int)(Math.pow(2, n)-1);
		int [] root = new int[numSamples];
		Arrays.fill(root, 50);
		for (int i=0; i < m; i++)
		{
			String fileName = "/Users/viq/smuth/SMutH/data/trees/tree"+n+"_"+i;
			numGroups = 0;
			groups = new HashMap<Integer, ArrayList<int[]>>();
			try {
				pw = new PrintWriter(new FileWriter(fileName));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("tree"+n+"_"+i);
			Arrays.fill(root, 50);
			expand(rootTag,root);
			completePrivateGroups();
			pw.close();
			System.out.println("num of groups: "+numGroups+" "+groups.size());
			if (numGroups < 3) i--; 
		}
	}
	
	private void expand(int nodeTag, int node[]){
		/*System.out.print("EXPAND:");
		for (int i=0; i<numSamples; i++){
			System.out.print(node[i]+"\t");
		}
		System.out.println();*/
		
		HashSet<Integer> children = new HashSet<Integer>();
		//children.add(0);
	
		
		Random generator = new Random();
		int randomIndex = generator.nextInt(100);
		/*with 5% probability same cluster!!!!
		if (randomIndex < 10){
			
		}*/
		
		
		
		
		boolean next;
		
		do{ 
			int child[] = new int[numSamples];
			boolean tag[] = new boolean[numSamples];
			child[0] = 0; 
			tag[0] = Boolean.FALSE; //normal sample
			
			int samp = 0;
			next = Boolean.FALSE;
			for (int i=1; i<numSamples; i++){
				randomIndex = 0;
				if (node[i] >= MIN_DIST_AAF)
					randomIndex = generator.nextInt(node[i]);
				if (randomIndex >= MIN_DIST_AAF){
					child[i] =  randomIndex;
					tag[i] = Boolean.TRUE;
					samp++;
					
				}else{
					child[i] = 0;
					tag[i] = Boolean.FALSE;
				}
				if (node[i] - child[i] >= MIN_DIST_AAF)
					next = Boolean.TRUE;
			}
				
			int c = booleansToInt(tag);
			if (c==0) return;
			
			if (groups.containsKey(c)){
			  
			
				/*boolean similar = false;
				for (int[] group : groups.get(c)){
					similar = true; 
					for (int i=1; i<numSamples; i++){					
						if ( Math.abs(group[i]-child[i]) > MIN_DIST_AAF){
							similar = false; break;
						}
					}
					System.out.println(booleansToString(tag)+" "+similar);	
					if (similar) break; 
				}
				if (similar) continue;*/
				continue;
			} else {
				groups.put(new Integer(c), new ArrayList<int[]>());
			}
			groups.get(c).add(child);
			
			pw.write(booleansToString(tag)+"\t"+(generator.nextInt(MaxSize-MinSize)+MinSize)+"\t");
				for (int i=0; i<numSamples; i++){
					node[i] = node[i] - child[i];
					pw.write((double)child[i]/100 +"\t");
					//System.out.print(((double)child[i]/100 +"\t"));
				}
			pw.write("\n");
			//System.out.println();
			if (samp > 1) {numGroups++; expand(c, child);}
					
					
			next = next&&(generator.nextBoolean());
			 
		}while(next);
		
		//System.out.println("done");
		
	}
	
	private void completePrivateGroups(){
		for (int i=0; i<numSamples-1; i++){
			int c = (int)Math.pow(2, i);
			if (!groups.containsKey(c)){
				int[] child = new int[numSamples];
				boolean[] tag = new boolean[numSamples];
				child[numSamples-1-i] =  10;
				tag[numSamples-1-i] = true;
				
				pw.write(booleansToString(tag)+"\t"+(new Random().nextInt(90)+10)+"\t");
				for (int j=0; j<numSamples; j++){
					pw.write((double)child[j]/100 +"\t");
					//System.out.print(((double)child[j]/100 +"\t"));
				}
			pw.write("\n");
			
			
			}
		}
	}
	
	int booleansToInt(boolean[] arr){
	    int n = 0;
	    for (boolean b : arr)
	        n = (n << 1) | (b ? 1 : 0);
	    return n;
	}
	
	
	String booleansToString(boolean[] arr){
	    String n = "";
	    for (boolean b : arr)
	        n = n+ (b ? 1 : 0);
	    return n;
	}
	
	public static void main(String[] args) {
		for (int i=5; i<=15; i++ )
			new SimulateTree(i, 100);
		
		//new SimulateTree(4, 5);
		
	}

}
