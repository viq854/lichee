/**
 * Represents an SNV entry
 */
package lineage;

import java.util.ArrayList;
import util.CNVRegion;

public class SNVEntry {
	private String entry;
	private String chrom;
	private int pos;
	private char ref;
	private char alt;	
	private String info;
	private double[] AAF;
	private int[] rc;
	private String[] genotype;
	private boolean robust;
	private boolean inCNVRegion;
	private SNVAnnotation annotation;
	 
	public SNVEntry(String line, int numofSamples){
		entry = line;
		AAF = new double[numofSamples];
		rc = new int[numofSamples];
		genotype = new String[numofSamples];
		robust = true;
		annotation = new SNVAnnotation();
		
		String[] entryParts = entry.split("\t");
		int i = 0;
		chrom = entryParts[i++];
		pos = new Integer(entryParts[i++]).intValue();
		ref = entryParts[i++].charAt(0);
		alt = entryParts[i++].charAt(0);	 
		info = entryParts[i++];
		
		int idx = 0;
		for (; i < entryParts.length; i++) {
			AAF[idx] = Double.parseDouble(entryParts[i]);
			if (AAF[idx] <  Parameters.VALIDATION_THR){
				genotype[idx] = "0/0";
				if (AAF[idx] >=  Parameters.VALIDATION_SOFT_THR) {
					robust = false;
				}
					
			} else { 
				genotype[idx] = "1/1";
			}
			idx++;
		}
	}
	
	public String getChromosome(){
		return chrom;
	}
	
	public int getChromNum(){
		if (chrom.charAt(3) == 'X') return 23;
		if (chrom.charAt(3) == 'Y') return 24;
		return new Integer(chrom.substring(3)).intValue();
	}
	
	
	public int getPosition(){
		return pos;
	}
	
	/**
	 * Returns the reference allele of the entry
	 */
	public char getRefChar(){
		return ref;
	}
	
	/**
	 * Returns the alternate allele of the entry
	 */
	public char getAltChar(){
		return alt;
	}
	
	/**
	 * Returns the particular genotype of a sample from the entry.
	 */
	public String getGenotype(int sample){
		return genotype[sample];
	}
	
	public boolean isRobust(){
		return robust;
	}
	
	/**
	 * Returns the group tag for a sample. The tag is
	 * found by checking whether the genotype is equivalent to "0/0"
	 * for each sample. If so, 0 is appended as the tag
	 * for that sample; otherwise, 1 is appended. At the end,
	 * one will have a binary code of length equal to the number of 
	 * samples for the entry.
	 */
	public String getGroup(){
		String result = "";
		for (int i = 0; i < genotype.length; i++){
			if (genotype[i].equals("0/0")) result += "0";
			else result += "1";
		}
		return result;
	}
	
	public double getAAF(int i) {
		return AAF[i];
	}
	
	public int getReadCount(int i) {
		return rc[i];
	}
	
	public void setAAF(int i, double aaf) {
		AAF[i] = aaf;
	}
	
	public String getEOATag() {
		return info;
	}
	
	public boolean evidenceOfPresence(int sample){
		return (AAF[sample] > Parameters.VALIDATION_SOFT_THR );
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
	
	public SNVAnnotation getAnnotation() {
		return annotation;
	}
	
	public void setAnnotation(SNVAnnotation ann) {
		annotation = ann;
	}
	
	public boolean checkInCNVRegion(ArrayList<CNVRegion> CNVs){
		for (CNVRegion cnv : CNVs){
			int loc = cnv.compareLocation(getChromNum(),pos);
			if (loc == 0){
				
				inCNVRegion = true;
				return true;
			}
			if (loc == -1){
				inCNVRegion = false;
				return false;
			}
		}
		inCNVRegion = false;
		return false;
	}
	
	public boolean isInCNVRegion(){
		return inCNVRegion;
	}
	
	public String toString(){
		return entry;
	}
	
	protected class SNVAnnotation {
		public String annovar;
		public String cosmic;
		public String tcga;
		
		public String toString() {
			return annovar + " " + cosmic + " " + tcga;
		}
	}
}
