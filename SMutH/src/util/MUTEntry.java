/**
 * Class: VCFEntry
 * Constructor: VCFEntry(String entry)
 * ----
 * This class represents a particular VCF entry from a 
 * VCF file. Note that all samples are 0-indexed.
 */
package util;


public class MUTEntry extends SNVEntry{
	
    private String EOAtag;
	private double[] AAF;
	
	public MUTEntry(String entry, int numofSamples){
		row = entry;
		
		
		String[] entryParts = entry.split("\t");
		AAF = new double[numofSamples];
		genotype = new String[numofSamples];
		robust = true;
		
		for (int i = 0; i < entryParts.length; i++){
			switch(i){
			case 0:
				chrom = entryParts[i];
				break;
			case 1:
				pos = new Integer(entryParts[i]).intValue();
				break;
			case 2:
				ref = entryParts[i].charAt(0);
				break;
			case 3:
				alt = entryParts[i].charAt(0);
				break;
			case 4:
				EOAtag = entryParts[i];
				/*String all0s = "";
				for (int x = 0; x<numofSamples - EOAtag.length(); x++) all0s +="0";
				EOAtag = all0s + EOAtag;*/
				break;
			default:
				AAF[i-5] = Double.parseDouble(entryParts[i]);
				if (AAF[i-5] <  Configs.VALIDATION_THR){
					genotype[i-5] = "0/0";
					if (AAF[i-5] >=  Configs.VALIDATION_SOFT_THR)
						robust = false;
						
				}else 
					genotype[i-5] = "1/1";
				
				break;
			}
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
