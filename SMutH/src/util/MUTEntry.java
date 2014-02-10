/**
 * Class: VCFEntry
 * Constructor: VCFEntry(String entry)
 * ----
 * This class represents a particular VCF entry from a 
 * VCF file. Note that all samples are 0-indexed.
 */
package util;

import util.Configs;
import util.SNVEntry;


public class MUTEntry extends SNVEntry{
	
    private String EOAtag;
	private double[] AAF;
	
	public MUTEntry(String entry, int numofSamples){
		row = entry;
		
		String[] entryParts = entry.split("\t");
		AAF = new double[numofSamples];
		genotype = new String[numofSamples];
		robust = true;
		
		int i =0;
		chrom = entryParts[i++];
		pos = new Integer(entryParts[i++]).intValue();
		ref = entryParts[i++].charAt(0);
		alt = entryParts[i++].charAt(0);
		EOAtag = entryParts[i++];
		/*String all0s = "";
		for (int x = 0; x<numofSamples - EOAtag.length(); x++) all0s +="0";
		EOAtag = all0s + EOAtag;*/
		int AAF_suffix = i;
		for (; i < entryParts.length; i++){
				AAF[i-AAF_suffix] = Double.parseDouble(entryParts[i]);
				if (AAF[i-AAF_suffix] <  Configs.VALIDATION_THR){
					genotype[i-AAF_suffix] = "0/0";
					if (AAF[i-AAF_suffix] >=  Configs.VALIDATION_SOFT_THR)
						robust = false;
						
				}else 
					genotype[i-AAF_suffix] = "1/1";
		}

	}
	
	@Override
	public double getAAF(int i) {
		return AAF[i];
	}
	
	public String getEOATag() {
		return EOAtag;
	}
	
	public boolean EvidenceOfPresence(int sample){
		return (AAF[sample] > Configs.VALIDATION_SOFT_THR );
	}
	
	public boolean EvidenceOfAbsence(int sample){
		return false;//(AAF[sample] < Configs.VALIDATION_HARD_THR );
	}
	
	public void updateGroup(String code){
		for (int i = 0; i < genotype.length; i++){
			
			if (genotype[i].equals("0/0") && code.charAt(i) =='1') {
				genotype[i] = "0/1";
			}else if (!genotype[i].equals("0/0") && code.charAt(i) =='0') {
				genotype[i] = "0/0";
			}
		}
	}
	
}
