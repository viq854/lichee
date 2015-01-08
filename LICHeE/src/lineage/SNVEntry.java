/*
 * Program LICHeE for multi-sample cancer phylogeny reconstruction
 * by Victoria Popic (viq@stanford.edu) 2014
 *
 * MIT License
 *
 * Copyright (c) 2014 Victoria Popic.
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/


/**
 * Represents an SNV entry
 */
package lineage;

import java.util.ArrayList;

import lineage.Parameters.Format;
import util.CNVRegion;

public class SNVEntry {
	
	// Format information
	protected static final int NUM_REQ_FIELDS_MUT = 5;
	protected static final int NUM_REQ_FIELDS_MUTC = 6;
	
	// Input file entry
	private String entryLine;
	
	// Input fields
	private String chromosome;
	private int position;
	private char refAllele;
	private char altAllele;	
	private String info;
	private String presenceProfile;
	
	// VAFs
	private double[] AAF;
	
	private String[] genotype;
	private int[] rc;
	private boolean robust;
	private boolean inCNVRegion;
	private SNVAnnotation annotation;
	 
	public SNVEntry(String line, int numSamples){
		entryLine = line;
		String[] entryParts = entryLine.split("\t");
		int numFields = Parameters.INFORMAT == Format.MUTC ? NUM_REQ_FIELDS_MUTC : NUM_REQ_FIELDS_MUT;
		if(entryParts.length != (numFields + numSamples)) {
			System.err.println("ERROR: Wrong input format: expecting " + numFields + " fields and " + numSamples + " samples" + " in entry: " + entryLine);
			System.exit(1);
		}
		
		// parse fields
		chromosome = entryParts[0].trim();
		try {
			position = Integer.parseInt(entryParts[1].trim());
		} catch (NumberFormatException e) {
			System.err.println("ERROR: Wrong input position format: " + entryParts[1].trim() + " in " + entryLine);
			System.exit(1);
		}
		refAllele = entryParts[2].charAt(0);
		altAllele = entryParts[3].charAt(0);	 
		info = entryParts[4].trim();
		if(Parameters.INFORMAT == Format.MUTC) {
			presenceProfile = entryParts[5].trim();
			if(presenceProfile.length() != numSamples) {
				System.err.println("ERROR: Wrong input presence profile format: " + presenceProfile + " -- the length does not match the number of input samples" + " in entry: " + entryLine);
				System.exit(1);
			}
		}
		
		// parse sample frequency values
		AAF = new double[numSamples];
		genotype = new String[numSamples];
		robust = true;
		for (int i = 0; i < numSamples; i++) {
			try {
				AAF[i] = Double.parseDouble(entryParts[i + numFields]);
			} catch (NumberFormatException e) {
				System.err.println("ERROR: Wrong input VAF value: " + entryParts[i + numFields] + " in entry: " + entryLine);
				System.exit(1);
			}
			if (AAF[i] <  Parameters.VALIDATION_THR){
				genotype[i] = "0/0";
				if (AAF[i] >=  Parameters.VALIDATION_SOFT_THR) {
					robust = false;
				}	
			} else { 
				genotype[i] = "1/1";
			}
		}
		
		rc = new int[numSamples];
		annotation = new SNVAnnotation();	
	}
	
	/** Returns the SNV chromosome */
	public String getChromosome(){
		return chromosome;
	}
	
	/** Returns the SNV chromosome as an integer value */
	public int getChromNum(){
		if (chromosome.charAt(3) == 'X') return 23;
		if (chromosome.charAt(3) == 'Y') return 24;
		try {
			return Integer.parseInt(chromosome);
		} catch (Exception e) {
			System.err.println("ERROR: Wrong input chromosome format " + chromosome + " in entry: " + entryLine);
			System.exit(1);
			return 0;
		}
	}
	
	/** Returns the SNV position  */
	public int getPosition(){
		return position;
	}
	
	/** Returns the reference allele */
	public char getRefChar(){
		return refAllele;
	}
	
	/** Returns the alternate allele */
	public char getAltChar(){
		return altAllele;
	}
	
	/** Returns the particular genotype of a sample from the entry */
	public String getGenotype(int sample){
		return genotype[sample];
	}
	
	/** Returns the VAF in sample i */
	public double getAAF(int i) {
		return AAF[i];
	}
	
	public int getReadCount(int i) {
		return rc[i];
	}
	
	/** Returns the additional info field */
	public String getInfoField() {
		return info;
	}
	
	public boolean evidenceOfPresence(int sample){
		return (AAF[sample] > Parameters.VALIDATION_SOFT_THR );
	}
	
	/** Returns true if the SNV was robustly called in all samples  */
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
	public String getGroupProfile() {
		if(presenceProfile != null) {
			return presenceProfile;
		}
		
		String result = "";
		for (int i = 0; i < genotype.length; i++){
			if (genotype[i].equals("0/0")) result += "0";
			else result += "1";
		}
		return result;
	}
	
	public void updateGroup(String code) {
		for (int i = 0; i < genotype.length; i++) {
			if (genotype[i].equals("0/0") && code.charAt(i) == '1') {
				genotype[i] = "0/1";
			} else if (!genotype[i].equals("0/0") && code.charAt(i) == '0') {
				genotype[i] = "0/0";
			}
		}
	}
	
	// Annotation handling
	
	public SNVAnnotation getAnnotation() {
		return annotation;
	}
	
	public void setAnnotation(SNVAnnotation ann) {
		annotation = ann;
	}
	
	public boolean hasAnnotation() {
		return !annotation.isEmpty();
	}
	
	public boolean checkInCNVRegion(ArrayList<CNVRegion> CNVs){
		for (CNVRegion cnv : CNVs){
			int loc = cnv.compareLocation(getChromNum(),position);
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
		return entryLine;
	}
	
	protected class SNVAnnotation {
		public String annovar;
		public String cosmic;
		public String tcga;
		
		public String toString() {
			return annovar + " " + cosmic + " " + tcga;
		}
		
		public boolean isEmpty() {
			return annovar == null;
		}
	}
}
