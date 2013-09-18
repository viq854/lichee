package unmixing;

import java.util.ArrayList;

public class LOHEntry {

	private String chrom;
	private String start;
	private String end;
	private ArrayList<Double> LAFs;
	
	
	public LOHEntry(String entry){
		int numofSamples = 0;
		LAFs =  new ArrayList<Double>(numofSamples);
	}
	
	public LOHEntry(double[] inputLAFs){
		int numofSamples = 0;
		LAFs =  new ArrayList<Double>(numofSamples);
	}
	
	public void setChrom(String chrom) {
		this.chrom = chrom;
	}
	public String getChrom() {
		return chrom;
	}
	public void setStart(String start) {
		this.start = start;
	}
	public String getStart() {
		return start;
	}
	public void setEnd(String end) {
		this.end = end;
	}
	public String getEnd() {
		return end;
	}
}
