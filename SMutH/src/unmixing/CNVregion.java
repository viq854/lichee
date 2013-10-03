package unmixing;

import io.*;

public class CNVregion {

	private int startChr;
	private int startPos;

	private int endChr;
	private int endPos;

	private double[] avgLAF;
	private double[] avgDepth;
	
	
	public CNVregion( int schr, int s, int echr, int e,double[] LAF, double[] Dep){
		startChr = schr;
		startPos = s;
		endChr = echr;
		endPos = e;
		avgLAF = LAF.clone();
		avgDepth = Dep.clone();
		System.out.print("NEW "+this+"\n");
	}
	
	public void append(int echr,int epos, double[] LAF, double[] Dep){
		
		for (int i=0; i< LAF.length;i++){
			avgLAF[i] = (avgLAF[i]+LAF[i])/2;
			avgDepth[i] = (avgDepth[i]+Dep[i])/2;
		}
		endChr = echr;
		endPos = epos;
		System.out.print("APPEND "+this+"\n");
	}


	public void setAvgLAFs(double[] avgLAFs) {
		this.avgLAF = avgLAFs;
	}

	public double[] getAvgLAF() {
		return avgLAF;
	}


	public double[] getAvgDepth() {
		return avgDepth;
	}



	public double getAvgDepth(int sample) {
		return avgDepth[sample];
	}
	

	public boolean isLoss() {
		for (int i=0; i < avgDepth.length; i++)
			if (avgDepth[i] < 1-Unmixing.Depth_ERROR)
				return true;
		return false;
	}
	
	public String toString(){
		String str = startChr +" "+ startPos+" "+ endChr+" "+endPos;
		for (int i=0; i < avgDepth.length; i++)
			str += " "+avgLAF[i];
		
		return str;
	}

	public int getEndPos() {
		// TODO Auto-generated method stub
		return endPos;
	}
	
	/**
	 * -1 before, 0 inside, 1 after
	 * @param entry
	 * @return
	 */
	public int compareLocation(VCFEntry entry){
		if (endChr < entry.getChromNum())
			return 1;
				
		if (startChr > entry.getChromNum())
			return -1;
		
		if (startChr == entry.getChromNum()){
			if (entry.getPosition() < startPos)
				return -1;
			
			if (endChr != entry.getChromNum() || entry.getPosition() < endPos)
				return 0;
			return 1;
		}
		
		if (endChr == entry.getChromNum() && entry.getPosition() > endPos)
				return 1;

		return 0;

	}

}
