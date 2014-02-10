package simulation;

import java.util.ArrayList;
import java.util.Comparator;

public class GroupInfo {
	public static int numofSamples;
	private String tag;
	private int size;
	private double[] centroid;
	
	public GroupInfo(String entry){
		String[] entryParts = entry.split("\t");
	
		centroid = new double[numofSamples];
		
		tag = entryParts[0];
		size = new Integer(entryParts[1]).intValue();
		
		if( entryParts.length-2 == numofSamples){
			for (int i = 0; i < numofSamples; i++)
				centroid[i] = Double.parseDouble(entryParts[i+2]);
		}else{
			int j = 2;
			for (int i = 0; i < numofSamples; i++){
				if (tag.charAt(i) == '0') centroid[i] = 0;
				else centroid[i] = Double.parseDouble(entryParts[j++]);
			}
		}
		
	}
		
	public String getTag() {
		return tag;
	}
	
	
	 public int getSize(){
		 return size;
	 }
	public double[] getCentroid(){
		return centroid;
	}
	
	public double L1disatnce(GroupInfo gi){
		double dist =0;
		for (int i=0; i<numofSamples; i++)
			dist += Math.abs(centroid[i] - gi.getCentroid()[i]);
		return dist;
	}
	
	public double L1disatnce(){
		double dist =0;
		for (int i=0; i<numofSamples; i++)
			dist += centroid[i];
		return dist;
	}
	
	public int L2disatnceClosest(ArrayList<GroupInfo> gil){
		double Mindist = L2disatnce(gil.get(0));
		int  Mini = 0;
		for (int i=1; i<gil.size(); i++)
			if (L2disatnce(gil.get(i)) < Mindist){
				Mini = i;
				Mindist = L2disatnce(gil.get(i));
			}
		return Mini;
	}
	
	public double L2disatnce(GroupInfo gi){
		double dist =0;
		for (int i=0; i<numofSamples; i++){
			
			dist += Math.pow(centroid[i] - gi.getCentroid()[i],2.0);
			//System.out.println(dist +"\t"+centroid[i] +" "+ gi.getCentroid()[i]);
		}
		return dist;
	}
	
	public double L2disatnce(){
		double dist =0;
		for (int i=0; i<numofSamples; i++)
			dist += Math.pow(centroid[i] ,2.0);
		return dist;
	}
}
