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
	
	
	public SimulateTree(int n, int m){
		numSamples = n;
		int rootTag = (int)(Math.pow(2, n)-1);
		int [] root = new int[numSamples];
		Arrays.fill(root, 50);
		for (int i=0; i < m; i++)
		//int i = 4;
		{
			String fileName = "/Users/rahelehs/Work/cancertree/LineageTree/simulation/trees/tree"+n+"_"+i;
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
			if (numGroups < 3) i--; 
		}
	}
	
	private void expand(int nodeTag, int node[]){
		System.out.print("EXPAND:");
		for (int i=0; i<numSamples; i++){
			System.out.print(node[i]+"\t");
		}
		System.out.println();
		HashSet<Integer> children = new HashSet<Integer>();
		//children.add(0);
	
		
		Random generator = new Random();
		boolean next;
		
		do{ 
			int child[] = new int[numSamples];
			boolean tag[] = new boolean[numSamples];
			child[0] = 0; 
			tag[0] = Boolean.FALSE; //normal sample
			
			int samp = 0;
			next = Boolean.FALSE;
			for (int i=1; i<numSamples; i++){
				int randomIndex = 0;
				if (node[i] >= 10)
					randomIndex = generator.nextInt(node[i]);
				if (randomIndex >= 10){
					child[i] =  randomIndex;
					tag[i] = Boolean.TRUE;
					samp++;
					if (node[i] - child[i] >= 10)
						next = Boolean.TRUE;
				}else{
					child[i] = 0;
					tag[i] = Boolean.FALSE;
				}
			}
			
			int c = booleansToInt(tag);
			if (c==0) return;
			
			if (groups.containsKey(c)){
				boolean similar = false;
				for (int[] group : groups.get(c)){
					similar = true; 
					for (int i=1; i<numSamples; i++){					
						if ( Math.abs(group[i]-child[i]) > 10){
							similar = false; break;
						}
					}
					if (similar) break; 
				}
				if (similar) continue;
			}else{
				groups.put(new Integer(c), new ArrayList<int[]>());
			}
			groups.get(c).add(child);
			
			pw.write(booleansToString(tag)+"\t"+(generator.nextInt(180)+20)+"\t");
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
				
				pw.write(booleansToString(tag)+"\t"+(new Random().nextInt(180)+20)+"\t");
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
		for (int i=3; i<=15; i++ )
			new SimulateTree(i, 100);
		
		//new SimulateTree(4, 5);
		
	}

}
