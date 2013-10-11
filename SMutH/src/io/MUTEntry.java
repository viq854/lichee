/**
 * Class: VCFEntry
 * Constructor: VCFEntry(String entry)
 * ----
 * This class represents a particular VCF entry from a 
 * VCF file. Note that all samples are 0-indexed.
 */
package io;


public class MUTEntry extends SNVEntry{
	
    private String EOAtag;
	private double[] AAF;
	
	public MUTEntry(String entry, int numofSamples){
		row = entry;
		
		
		String[] entryParts = entry.split("\t");
		AAF = new double[numofSamples];
		genotype = new String[numofSamples];
		
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
				break;
			default:
				AAF[i-5] = Double.parseDouble(entryParts[i]);
				if (AAF[i-5] <=  VCFConstants.VALIDATION_THR)
					genotype[i-5] = "0/0";
				else 
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
	
	
}
