package simulation;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import util.SNVEntry;

public class Evaluate {
	public static String pathS;
	public static String pathT;
	
	public int TP_group;

	public int TP_mut;

	public double msError;
	

	private HashMap<String,ArrayList<GroupInfo>> giSim, giTest;
	private GroupSummary sim, test;
	
	
	public class GroupSummary{
		public int mutations;
		public int groups;
		
		public GroupSummary(){
			mutations=0;
			groups=0;
		}
	}
	public static class GroupInfoComparator implements Comparator<GroupInfo> {

        public int compare(GroupInfo o1, GroupInfo o2) {
            return o1.getTag().compareTo(o2.getTag());
        }
    }	
	
	public Evaluate(int n, int id){
		
		giSim = new HashMap<String,ArrayList<GroupInfo>>();
		sim = new GroupSummary();
		load(pathS+"/tree"+n+"_"+id, giSim, sim);
		giTest = new HashMap<String,ArrayList<GroupInfo>>();
		test = new GroupSummary();
		load(pathT+"/tree"+n+"_"+id+".txt.nodes",giTest,test);
		
		distance();
		//System.out.println("tree"+n+"_"+id+"\t"+((double)TP_group/sim_group)+"\t"+((double)TP_group/test_group)+"\t"+dist);
	
		
	}
	
	public void load(String inputFile, HashMap<String,ArrayList<GroupInfo>> gis, GroupSummary gs){
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String line = rd.readLine();
			while(line != null){
				GroupInfo gi = new GroupInfo(line);
				gs.groups++;
				gs.mutations += gi.getSize();
				if (!gis.containsKey(gi.getTag())){
					gis.put(gi.getTag(), new ArrayList<GroupInfo>());				
				}
				gis.get(gi.getTag()).add(gi); 
			
				line = rd.readLine();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void distance(){
		ArrayList<String> codes1 = new ArrayList<String>(giSim.keySet());
		ArrayList<String> codes2 = new ArrayList<String>(giTest.keySet());

		msError = 0;
		TP_group = 0;
		
		for (int i = 0; i < codes1.size(); i++){
			ArrayList<GroupInfo> a1 = giSim.get(codes1.get(i));
			if (codes2.contains(codes1.get(i))){
			
				ArrayList<GroupInfo> a2 = (ArrayList<GroupInfo>) giTest.get(codes1.get(i)).clone();
				TP_group += Math.min(a1.size(), a2.size());
				
				//System.out.println(codes1.get(i)+" "+a1.size()+" "+ a2.size());
				for (GroupInfo gi : a1){
					GroupInfo closestgi = a2.remove(gi.L2disatnceClosest(a2));
					
					msError += Math.pow(gi.L2disatnce(closestgi),0.5)/GroupInfo.numofSamples;
					TP_mut += Math.min(gi.getSize(), closestgi.getSize());
					//System.out.println("dist = "+dist);
					if (a2.isEmpty()) break;
				}
				
			}else{
				//FN_mut += gis1.get(codes1.get(i)).size();
				/*for (GroupInfo gi : a1){
					//dist += gi.L2disatnce();
				}*/
				
			}
			
		}	
		
		msError = Math.pow(msError, 0.5)/TP_group;
	}
	
	
	public int getSimGroups(){
		return sim.groups;
	}
	
	public int getSimMutations(){
		return sim.mutations;
	}
	
	public int getTestGroups(){
		return test.groups;
	}
	public int getTestMutations(){
		return test.mutations;
	}
	
	
	public static void main(String[] args){
		pathS = "/Users/rahelehs/Work/cancertree/LineageTree/simulation/trees";
		pathT = "/Users/rahelehs/Work/cancertree/LineageTree/simulation/trees/var.01";
		//System.out.println("\tPrecision\tSensitivity\tsq-error");
		
		for (int j =3; j<10;j++){
			GroupInfo.numofSamples = j;
			
			double error = 0;
			double sensitivity = 0, sensitivityM = 0;
			double precision = 0, precisionM = 0; 
			int sim_group = 0, test_group = 0;
			int sim_mut = 0, test_mut =0;
			
			//System.out.println("\tPrecision\tSensitivity\tsq-error");
			for (int i =0; i<100;i++){
				Evaluate e  = new Evaluate(j,i);
				
				precision += ((double)e.TP_group/e.getSimGroups());
				sensitivity += ((double)e.TP_group/e.getTestGroups());
				sim_group += e.getSimGroups();
				test_group += e.getTestGroups();
								
				precisionM += ((double)e.TP_mut/e.getSimMutations());
				sensitivityM += ((double)e.TP_mut/e.getTestMutations());
				sim_mut += e.getSimMutations();
				test_mut += e.getTestMutations();
				
				error += e.msError/e.TP_group;

			}
		
			
		
			System.out.println(precision+"\t"+sensitivity+"\t"+sim_group/100.00+"\t"+test_group/100.00+"\t"+
					precisionM+"\t"+sensitivityM+"\t"+sim_mut/100.00+"\t"+test_mut/100.00+"\t"+error/100);
		}
	}

}
