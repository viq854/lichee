/**
 * Class: FLEntry
 * Constructor: FLEntry(String entry)
 * ----
 * This class represents a particular mutation entry from a 
 * text file in a specific format for Folecular Lymphoma project.
 * #
 */
package util;


public class FLEntry extends SNVEntry{
	
	private char[] call;

	private int[] refCount;
	private int[] altCount;
	
	public FLEntry(String entry, int numofSamples){
		row = entry;
		
		
		String[] entryParts = entry.split("\t");
		call = new char[numofSamples];
		refCount = new int[numofSamples];
		altCount = new int[numofSamples];

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
			default:
				String[] freqParts = entryParts[i].split(":");
				String[] alleleDepths = freqParts[1].split(",");
				call[i-4] = freqParts[0].charAt(0);
				refCount[i-4] = Integer.parseInt(alleleDepths[0]);
				altCount[i-4] = Integer.parseInt(alleleDepths[1]);
				
				
				if (freqParts[0].equals("0") || freqParts[0].equals("Freq")){
					genotype[i-4] = "0/0";
					if (EvidenceOfPresence(i-4))
						robust = false;
						
				}else{
					genotype[i-4] = "1/1";
					if (EvidenceOfAbsence(i-4)){
						//System.out.println(entry);
					
						robust = false;
				}
					}
				
				break;
			}
		}

	}

	
	
	public boolean EvidenceOfPresence(int sample){
		//return (call[sample] == 'F');
		return (getAAF(sample) > 0.1);
				
		/*int a = altCount[sample];
		int d = refCount[sample]+altCount[sample];
		double total = 0.0;
		for (int i = d; i >= a; i--){
			total += getProb(sample, d, i);
		}
		return (total >= Configs.EDIT_PVALUE);*/
	}

	
	public boolean EvidenceOfAbsence(int sample){
		return (call[sample] == '1' && getAAF(sample) < 0.1 && altCount[sample] < 5);
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
	
	
	public double getAAF(int i) {
		return (double)altCount[i]/(refCount[i]+altCount[i]);
	}
	
	public double getLAF(int i) {
		return (double)(refCount[i] < altCount[i]? refCount[i]:altCount[i])/(double)(refCount[i]+altCount[i]);
	}
}
