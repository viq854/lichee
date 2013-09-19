/**
 * Class: VCFEntry
 * Constructor: VCFEntry(String entry)
 * ----
 * This class represents a particular VCF entry from a 
 * VCF file. Note that all samples are 0-indexed.
 */
package io;


import java.util.ArrayList;

public class VCFEntry {
	
	
	/* Private Instance Variables */
	private String raw;
	private String chrom;
	private String pos;
	private String id;
	private char ref;
	private char alt;
	private Double qual;
	private String filter;
	private String info;
	private String format;
	private ArrayList<String> alleleFreqList;
	
		
	private String[] genotype;
	private int[] refCount;
	private int[] altCount;
	
	/**
	 * Function: VCFEntry(String entry)
	 * Usage: (Constructor)
	 * ----
	 * Creates a new VCFEntry.
	 * 
	 * @param entry	The VCF entry line as a string from the VCF file.
	 * 
	 */
	
	public VCFEntry(String entry, int numeofSamples){
		raw = entry;
		alleleFreqList = new ArrayList<String>();
		String[] entryParts = entry.split("\t");
		for (int i = 0; i < entryParts.length; i++){
			switch(i){
			case 0:
				chrom = entryParts[i];
				break;
			case 1:
				pos = entryParts[i];
				break;
			case 2:
				id = entryParts[i];
				break;
			case 3:
				ref = entryParts[i].charAt(0);
				break;
			case 4:
				alt = entryParts[i].charAt(0);
				break;
			case 5:
				qual = Double.valueOf(entryParts[i]);
				break;
			case 6:
				filter = entryParts[i];
				break;
			case 7:
				info = entryParts[i];
				break;
			case 8:
				format = entryParts[i];
				break;
			default:
				String[] freqParts = entryParts[i].split(":");
				genotype[i-9] = freqParts[0];
				String[] alleleDepths = freqParts[1].split(",");
				refCount[i-9] = Integer.parseInt(alleleDepths[0]);
				altCount[i-9] = Integer.parseInt(alleleDepths[1]);

				break;
			}
		}

	}
	
	
	/**
	 * Function: getChromosome()
	 * Usage: String chrom = entry.getChromosome()
	 * ----
	 * Returns the chromosome of the entry as a string
	 * 
	 * @return	The chrom as a string
	 */
	public String getChromosome(){
		return chrom;
	}
	
	/**
	 * Function: getPosition()
	 * Usage: String pos = entry.getPos()
	 * ----
	 * Returns the position of the entry as a string
	 * @return
	 */
	public String getPosition(){
		return pos;
	}
	
	/**
	 * Function: getRefChar()
	 * Usage: char ref = entry.getRefChar()
	 * ----
	 * Returns the reference allele of the entry
	 * 
	 * @return the reference allele as a char
	 */
	public char getRefChar(){
		return ref;
	}
	
	/**
	 * Function: getAltChar()
	 * Usage: char alt = entry.getAltChar()
	 * ----
	 * Returns the alternate allele of the entry
	 * 
	 * @return the alternate allele as a char
	 */
	public char getAltChar(){
		return alt;
	}
	
	/**
	 * Function: getQuality()
	 * Usage: double qual = entry.getQuality()
	 * ----
	 * Returns the quality of the entry
	 * 
	 * @return the quality as a double
	 */
	public double getQuality(){
		return qual;
	}
	
	/**
	 * Function: getFilter()
	 * Usage: String filter = entry.getFilter()
	 * ----
	 * Returns the filter report of the entry
	 * 
	 * @return the filter as a String
	 */
	public String getFilter(){
		return filter;
	}
	
	/**
	 * Function: getInfo()
	 * Usage: String info = entry.getInfo()
	 * ----
	 * Returns the info column of the entry
	 * 
	 * @return the info column as a String
	 */
	public String getInfo(){
		return info;
	}
	
	/**
	 * Function: getFormat()
	 * Usage: String format = entry.getFormat()
	 * ----
	 * Returns the format column of the entry
	 * 
	 * @return the format column as a String
	 */
	public String getFormat(){
		return format;
	}
	
	
	/**
	 * Function: getAlleleFreq(int sample)
	 * Usage: String freq = entry.getAlleleFreq(sample)
	 * ----
	 * Returns the allele frequency of the particular sample of
	 * the entry
	 * 
	 * @param sample	the particular sample wanted from the entry
	 * @return	the raw allele frequency of the sample of the entry as a String
	 */
	public String getAlleleFreq(int sample){
		return alleleFreqList.get(sample);
	}
	
	/**
	 * Function: getFreqParts(int sample)
	 * Usage: String[] freqArray = entry.getFreqParts(sample)
	 * ----
	 * Returns an array of Strings which represent components
	 * from the frequency parts of the allele frequency. This is the same
	 * as the getAlleleFreq function, only splitting the string by ":"
	 * as a delimiter.
	 * @param sample	the particular sample wanted from the entry
	 * @return	the raw allele frequency of the sample of the entry as a String
	 */
	private String[] getFreqParts(int sample){
		return alleleFreqList.get(sample).split(":");
	}
	
	/**
	 * Function: getAlleleCount(int sample)
	 * Usage: String counts = entry.getAlleleCount(sample)
	 * ----
	 * Returns the counts of the major and minor allele in the format "#,#"
	 * 
	 * @param sample	the particular sample from the entry
	 * @return 	the major and minor allele of the sample of the entry as a String
	 */
	public int getAltCount(int sample){
		return altCount[sample];
	}
	
	
	/**
	 * Function: getReadDepth(int sample)
	 * Usage: int depth = entry.getReadDepth(sample)
	 * ----
	 * Returns the read depth of a particular sample of an entry
	 * @param sample	the particular sample from the entry
	 * @return	the depth of a sample of the entry as an int
	 */
	public int getReadDepth(int sample){
		return refCount[sample]+altCount[sample];
	}
	
	/**
	 * Function: getGenotype(int sample)
	 * Usage: String genotype = entry.getGenotype(sample)
	 * ----
	 * Returns the particular genotype of a sample from the entry.
	 * The genotype is generated by GATK.
	 * @param sample	the particular sample from the entry
	 * @return	the genotype of the sample of the entry as a string
	 */
	public String getGenotype(int sample){
		return genotype[sample];
	}
	
	/**
	 * Function: getGATK(int sample)
	 * Usage: String gatk = entry.getGATK(sample)
	 * ----
	 * Returns the GATK code for a sample. The GATK code is
	 * found by checking whether the genotype is equivalent to "0/0"
	 * for each sample. If so, 0 is appended as the GATK code
	 * for that sample; otherwise, 1 is appended. At the end,
	 * one will have a binary code of length equal to the number of 
	 * samples for the entry.
	 * @return	The GATK code for an entry as a string
	 */
	public String getGATK(){
		String result = "";
		for (int i = 0; i < alleleFreqList.size(); i++){
			if (genotype[i].equals("0/0")) result += "0";
			else result += "1";
		}
		return result;
	}

// -------Old Probability Formula-------
//	public double getConversionProb(int sample){
//		int a = getAlleleCount(sample, 1);
//		int d = getReadDepth(sample);
//		return nCr(d, a) * Math.pow(BASE_ERROR, a) * Math.pow((1 - BASE_ERROR), d - a);
//	}
// -------End-------
	
	/**
	 * Function: getSumProb(int sample)
	 * Usage: double sumProb = entry.getSumProb(sample)
	 * ----
	 * This function returns the probability that a given sample
	 * was called correctly by GATK. This part of the function
	 * handles the sigma and adding up all the individual
	 * probabilities from the read depth to the minor allele
	 * count.
	 * 
	 * @param sample	The particular sample of the entry
	 * @return	The summed probability of a read being correct as a double
	 */
	public double getSumProb(int sample){
		int a = altCount[sample];
		int d = getReadDepth(sample);
		double total = 0.0;
		for (int i = d; i >= a; i--){
			total += getProb(sample, d, i);
		}
		return total;
	}
	
	/**
	 * Function: getProb(int sample, int d, int k)
	 * Usage: double prob = entry.getProb(sample, d, k)
	 * ----
	 * This function calculates the probability of a sample
	 * being called incorrectly. Is used by getSumProb (which
	 * acts as a wrapper).
	 * 
	 * @param sample	The particular sample of the entry
	 * @param d			The allele depth of the sample
	 * @param k			The minor allele count
	 * @return	The probability as a double
	 */
	private double getProb(int sample, int d, int k){
		return nCr(d, k) * Math.pow(VCFConstants.BASE_ERROR / 4.0, k) * Math.pow((1 - VCFConstants.BASE_ERROR), d - k);
	}
	
	/**
	 * Function: nCr(int n, int k)
	 * Usage: double result = entry.nCr(n, k)
	 * ----
	 * This calculates the result of "n choose k". In other words, 
	 * it returns the number of combinations to choose k from 
	 * n objects.
	 * @param n	The number of total objects
	 * @param k	How many objects chosen
	 * @return	n choose k as a double
	 */
	private double nCr(int n, int k) {
	    if (k < 0 || k > n) return 0;
	    if (k > n/2) k = n - k;
	    double denominator = 1.0, numerator = 1.0;
	    for (int i = 1; i <= k; i++) {
	        denominator *= i;
	        numerator *= (n + 1 - i);
	    }
	    return numerator / denominator;
	}

	/**
	 * Function: toString()
	 * Usage: String rawEntry = entry.toString()
	 * ----
	 * Returns the VCF entry in string form.
	 * This is the same argument that was input
	 * into the constructor of the class.
	 * 
	 * @return	The raw VCF entry as a string
	 */
	public String toString(){
		return raw;
	}
	
	public void updateGATK(String code){
		for (int i = 0; i < alleleFreqList.size(); i++){
			String parts[] = alleleFreqList.get(i).split(":");
			String result = "";
			if (parts[0].equals("0/0") && code.charAt(i) =='1') {
				result = "0/1";
			}else if (!parts[0].equals("0/0") && code.charAt(i) =='0') {
				result = "0/0";
			}
			
			if (result != ""){
				for (int j=1; j < parts.length; j++){
					result = result + ":"+parts[j];
				}
				alleleFreqList.set(i, result);
			}
		}
	}

	public String getId() {
		return id;
	}
	
	
	public double getAAF(int i) {
		return (double)altCount[i]/(refCount[i]+altCount[i]);
	}
}
