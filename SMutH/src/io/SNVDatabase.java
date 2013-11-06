/**
 * Class: VCFDatabase
 * Constructor: VCFDatabase(String TESTFILE)
 * Last Edited: September 13, 2012
 * ----
 * This class is a container for VCFEntries.
 * The constructor takes the path name of the VCF file,
 * constructs a VCFEntry for each record, and stores it
 * internally. One can get back VCFEntries using their
 * GATK codes.
 */
package io;


import java.io.*;
import java.util.*;

import SMutH.TreeChecker;

import unmixing.CNVregion;


public class SNVDatabase {
	
	/* Private Instance Variables */
	private ArrayList<SNVEntry> somaticSNVs;
	//private ArrayList<VCFEntry> hgSNVs;
	private ArrayList<String> names; 
	private int allCounter = 0;
	private	int germlineCounter = 0;
	HashMap<String, ArrayList<SNVEntry>> TAG2SNVs;
	HashMap<String, Integer> TAG2RobustSNVNum;
	
	
	/**
	 * Function: VCFDatabase(String TESTFILE)
	 * Usage: (Constructor)
	 * ----
	 * This is the constructor for the VCFDatabase.
	 */
	public SNVDatabase(String inputFile, int normalSample){
		
		somaticSNVs = new ArrayList<SNVEntry>();
		names = new ArrayList<String>();
		System.out.println(inputFile);
		
		String[] fileParts = inputFile.split("\\.");
		
		if (fileParts[fileParts.length-1].equals("vcf"))
			loadVCF(inputFile,normalSample);
		else
			loadMUT(inputFile,normalSample);
			
		generateMap();
		
		System.out.println("Original Mutation Map");
		reportGroups();
		
		// Trying to resolve some conflicts
		int robustGroupSizeThreshold = TreeChecker.getGroupSizeThreshold(TAG2RobustSNVNum.values().size(), TAG2RobustSNVNum.size());	
		System.out.println("Group Size Threshold is "+ robustGroupSizeThreshold);
		
		Set<String> conflicts = new HashSet<String>();
	
		
		String all1s = "", all0s="";
		for (int i = 0; i < names.size(); i++){ all1s += "1";all0s += "0";}
		//SNVs which did not passed filters to be called in any group are conflicts
		conflicts.add(all0s);
		//Germline as valid group
		TAG2SNVs.put(all1s, new ArrayList<SNVEntry>());	
		Set<String> codes = new HashSet<String>(TAG2SNVs.keySet());
			
		//small groups are conflicts
		/*for (String code :codes){
			if (TAG2SNVs.get(code).size() <= groupSizeThreshold){
				conflicts.add(code);
			}
		}*/
		
		//non robust groups are conflicts
		for (String code :codes){
			if (!TAG2RobustSNVNum.containsKey(code) /*|| TAG2RobustSNVNum.get(code).intValue() < robustGroupSizeThreshold*/){
				conflicts.add(code);
			}
		}	
		
		//Edit conflicts to robust groups
		for (String conflict: conflicts){
			codes.remove(conflict);
		}
				
		codes.add(all1s);
		editSNVs(conflicts, codes);
		System.out.println("Mutation Map after editing to solid groups!");
		reportGroups();
		
		
		//Merge conflicts
		editSNVs(conflicts,conflicts);
		System.out.println("Mutation Map after merging conflict groups!");
		reportGroups();
		
	}
	
	public SNVDatabase(String inputFile, String hgFile, int normalSample){
		BufferedReader rd;
		somaticSNVs = new ArrayList<SNVEntry>();
		names = new ArrayList<String>();
		try{
			rd = new BufferedReader(new FileReader(inputFile));
			PrintWriter pw = new PrintWriter(new FileWriter(hgFile));
			String currLine = rd.readLine();
			String lastLine = "";
			while (currLine.substring(0, 1).equals("#")){ lastLine = currLine; currLine = rd.readLine();}
			//System.out.println(lastLine+"\n");
			getNames(lastLine);
			int numSamples = names.size();
			
			while (currLine != null){
				VCFEntry entry = new VCFEntry(currLine,numSamples);
				System.out.println(currLine);
				currLine = rd.readLine();
				
				
				//GATK filter
				//if (!entry.getFilter().equals("PASS") || entry.getQuality() < VCFConstants.MIN_QUAL)
				//	continue;
				allCounter++;
				
				//Makes sure entries are legitimate.
				boolean isLegitimate = true;
				boolean isGermline = false;
				int totalCoverage = 0;
				
				for (int i = 0; i < numSamples; i++){
					/* TO FILTER OUT SNVs with very low coverage in some SAMPLES*/
					if (entry.getGenotype(i).equals("./.") || entry.getReadDepth(i) <= VCFConstants.MIN_COVERAGE){
						isLegitimate = false;
						break;
					}
					totalCoverage += entry.getReadDepth(i);
					
					/*if(!entry.getGenotype(i).equals("0/1") && !entry.getGenotype(i).equals("1/1")){
						isGermline = false;
						break;
					}*/
				}
				//if (!isLegitimate) continue;
				/* TO FILTER OUT SNVs with low average coverage in all SAMPLES*/
				/*if (totalCoverage/numSamples <= VCFConstants.AVG_COVERAGE){
					isLegitimate = false;
				}*/
				if (!isLegitimate) continue;
				
				/* Germline mutations */
				if(!entry.getGenotype(normalSample).equals("0/0")){
					//Heterozygous Germline SNVs
					
					if( entry.getGenotype(normalSample).equals("0/1") && 
							entry.getAAF(normalSample) > VCFConstants.HETEROZYGOUS &&
							!entry.getId().equals(".")){
						//hgSNVs.add(entry);
						//System.out.println("h "+entry);
						
						pw.write(entry.getChromNum()+"\t"+entry.getPosition());
						for (int i = 0; i < names.size(); i++){
							//if (i != VCFConstants.NormalSample)
								pw.write("\t"+entry.getRefCount(i)+"\t"+entry.getAltCount(i));
						}
						pw.write("\n");
						
					}
					isGermline = true;
				}else 
				if (entry.getSumProb(normalSample) < VCFConstants.EDIT_PVALUE ){
					//System.out.println("*"+entry);
					isGermline = true;
				}
				
				if (isGermline) {
					germlineCounter++;
					continue;
				}
				
				//if (!entry.getId().equals("."))
					//continue;
				
				somaticSNVs.add(entry);
				
			}
			rd.close();
			pw.close();
			System.out.println("There are " + allCounter + " SNVs PASS by GATK hard filters. Of those, we pass "+ germlineCounter + " as germline, and "+ somaticSNVs.size() +" as somatic. \n");
		} catch (IOException e){
			System.out.println("File Reading Error!");
		}
		
		generateMap();
	}
	
	private void loadVCF(String inputFile, int normalSample){
		//"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
		try{
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String currLine = rd.readLine();
			String lastLine = "";
			while (currLine.substring(0, 1).equals("#")){ lastLine = currLine; currLine = rd.readLine();}
			//System.out.println(lastLine+"\n");
			getNames(lastLine);
			int numSamples = names.size();
			
			while (currLine != null){
				VCFEntry entry = new VCFEntry(currLine,numSamples);
				//System.out.println(currLine);
				currLine = rd.readLine();
				
				
				//GATK filter
				//if (!entry.getFilter().equals("PASS") || entry.getQuality() < VCFConstants.MIN_QUAL)
					//continue;
				allCounter++;
				
				//Makes sure entries are legitimate.
				boolean isLegitimate = true;
				boolean isGermline = true;
				int totalCoverage = 0;
				
				for (int i = 0; i < numSamples; i++){
					/* TO FILTER OUT SNVs with very low coverage in some SAMPLES*/
					if (entry.getGenotype(i).equals("./.") || entry.getReadDepth(i) <= VCFConstants.MIN_COVERAGE){
						isLegitimate = false;
						break;
					}
					totalCoverage += entry.getReadDepth(i);
					
					/* Germline mutations - soft filtering */
					if(!entry.getGenotype(i).equals("0/1") && !entry.getGenotype(i).equals("1/1")){
						isGermline = false;
						break;
					}
				}
				//if (!isLegitimate) continue;
				/* TO FILTER OUT SNVs with low average coverage in all SAMPLES*/
				/*if (totalCoverage/numSamples <= VCFConstants.AVG_COVERAGE){
					isLegitimate = false;
				}*/
				if (!isLegitimate) continue;
				
				/* Germline mutations - hard filtering */
				/*if(!entry.getGenotype(normalSample).equals("0/0")){
					isGermline = true;
				}else 
				if (entry.getSumProb(normalSample) < VCFConstants.EDIT_PVALUE ){
					//System.out.println("*"+entry);
					isGermline = true;
				}*/
				
				if (isGermline) {
					germlineCounter++;
					continue;
				}
				
				//if (!entry.getId().equals("."))
					//continue;
				
				somaticSNVs.add(entry);
				
			}
			rd.close();
			System.out.println("There are " + allCounter + " SNVs PASS by GATK hard filters. Of those, we pass "+ germlineCounter + " as germline (Heterozygous), and "+ somaticSNVs.size() +" as somatic. \n");
		} catch (IOException e){
			System.out.println("File Reading Error!");
		}
		
		
	}
	
	private void loadMUT(String inputFile, int normalSample){
		//"#CHROM\tPOS\tTAG\tREF\tALT\t";
		try{
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String currLine = rd.readLine();
			String lastLine = "";
			while (currLine.substring(0, 1).equals("#")){ lastLine = currLine; currLine = rd.readLine();}
			//System.out.println(lastLine+"\n");
			getNames(lastLine);
			int numSamples = names.size();
			
			while (currLine != null){
				MUTEntry entry = new MUTEntry(currLine,numSamples);
				//System.out.println(currLine);
				currLine = rd.readLine();
				allCounter++;
				
				//Makes sure entries are legitimate.
				//boolean isLegitimate = false;
				boolean isGermline = true;
				//int totalCoverage = 0;
				
				for (int i = 0; i < numSamples; i++){
					/* TO FILTER OUT SNVs with very low coverage in some SAMPLES*/
					/*if ( entry.getReadDepth(i) <= VCFConstants.MIN_COVERAGE){
						isLegitimate = false;
						break;
					}*/
					/* TO FILTER OUT SNVs that has not been called in any samples SAMPLES*/
					/*if(entry.getGenotype(i).equals("0/1") || entry.getGenotype(i).equals("1/1")){
						isLegitimate = true; 
					}*/
					
					/* Germline mutations - soft filtering */
					if(entry.getGenotype(i).equals("0/0")){
						isGermline = false;
					}
				}

				//if (!isLegitimate) continue;
				
				
				if (isGermline) {
					germlineCounter++;
					continue;
				}
				System.out.println(entry+" *"+entry.getGroup());
				somaticSNVs.add(entry);
				
				
			}
			rd.close();
			System.out.println("There are " + allCounter + " SNVs in validation file. Of those, we pass "+ germlineCounter + " as germline, and "+ somaticSNVs.size() +" as somatic. \n");
		} catch (IOException e){
			System.out.println("File Reading Error!");
		}
		
	}
	
	public String getHeader(){
		
		String header = "Samples = ";//"#There are " + allCounter + " SNVs PASS by GATK hard filters. Of those, we pass "+ germlineCounter + " as germline, and "+ somaticSNVs.size() +" as somatic. \n";
		//header +="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
		
		for(int i =0; i< names.size();i++)
			header += names.get(i)+ "\t";
		header +="#CHROM\tPOS\tREF\tALT\tDEPTH\tAAF\t";

		return header;
	}
	
	private void getNames(String inputLine){
		
		String[] header = inputLine.split("\t");
		if(header[2].equals("ID"))
			names = new ArrayList<String>(Arrays.asList(header).subList(9, header.length));
		else if (header[4].equals("EOA"))
			names = new ArrayList<String>(Arrays.asList(header).subList(5, header.length));
		/*else {
			names = new ArrayList<String>();
			for (char c='A'; c<='Z'; c++)
				names.add(c+"");
		}*/
		System.out.println("There are "+names.size()+" samples!");
	}
	
	
	/**
	 * Function: getEntriesByGATK(String inputCode)
	 * Usage: ArrayList<VCFEntry> entries = db.getEntriesByGATK(inputCode)
	 * ----
	 * Searches the somaticSNVs for all VCFEntries with matching GATK codes
	 * and returns an ArrayList of such entries
	 * 
	 * @param inputCode	The GATK code with which to filter the entries
	 * @return	An ArrayList of VCFEntries with a GATK equivalent to inputCode
	 */
	public ArrayList<SNVEntry> getEntriesByGATK(String inputCode){
		return TAG2SNVs.get(inputCode);
	}
	
	public HashMap<String,ArrayList<SNVEntry>> getTAG2SNVsMap(){
		return TAG2SNVs;
	}
	
	public HashMap<String,ArrayList<SNVEntry>> generateFilteredTAG2SNVsMap(ArrayList<CNVregion> CNVs){
		HashMap<String, ArrayList<SNVEntry>> filteredTAG2SNVs = new HashMap<String, ArrayList<SNVEntry>>();
		
		int groupSizeThreshold = TreeChecker.getGroupSizeThreshold(somaticSNVs.size(), TAG2SNVs.size());
		
		
		String all1s="", all0s="";
		for (int i = 0; i < names.size(); i++) {all1s += "1";all0s +="0";}


		int codeLength = -1;
		int cnvID = 0;
		for (SNVEntry entry: somaticSNVs){
			String code = entry.getGroup();
			//Filter bad groups
			if (code.equals(all1s) || code.equals(all0s) /*|| !robustTAGs.containsKey(code)*/ || TAG2SNVs.get(code).size() <= groupSizeThreshold)
				continue;
		
			//Filter CNV regions
			if (CNVs != null && cnvID < CNVs.size()){
				int loc = CNVs.get(cnvID).compareLocation(entry);
				if (loc == 0) break;
				if (loc == 1) cnvID++;
			}
			
			
			if (codeLength == -1) codeLength = code.length();
			if (!filteredTAG2SNVs.containsKey(code)){
				filteredTAG2SNVs.put(code, new ArrayList<SNVEntry>());				
			}
			filteredTAG2SNVs.get(code).add(entry); 
		}
		
		
		return filteredTAG2SNVs;
	}
	
	public void printInfo(PrintWriter pw, String inputCode){
		double[] sum;
		int count = 0;
		sum = new double[names.size()];
		
		
		for (SNVEntry entry: somaticSNVs){
			if (entry.getGroup().equals(inputCode)){
				for (int i = 0; i < names.size(); i++){
					sum[i] += entry.getAAF(i);
				}
				count++;
			}
		}
		pw.write("#Group: "+inputCode + "\t" + count);
		for (int i = 0; i < names.size(); i++){
			pw.write("\t"+sum[i]/count);
		}
		pw.write("\n");
	}
	
	public void printEntriesByGATK(PrintWriter pw, String inputCode) throws IOException{
		printInfo(pw,inputCode);
		for (SNVEntry entry: somaticSNVs){
			if (entry.getGroup().equals(inputCode)){
				//pw.write(entry+"\n");
				pw.write(entry.getChromosome()+"\t"+entry.getPosition());
				for (int i = 0; i < names.size(); i++){
					if (inputCode.charAt(i) == '1') pw.write("\t"+entry.getAAF(i));
				}
				pw.write("\n");
			}
		}
	}
	
/*	public void printLOH(PrintWriter pw) throws IOException{
		String header ="#CHROM\tPOS";
		
		for(int i =0; i< names.size();i++)
			header += "\t"+names.get(i);
		
		//pw.write(header+"\n");
		
		for (VCFEntry entry: hgSNVs){
			pw.write(entry.getChromosome()+"\t"+entry.getPosition());
			for (int i = 0; i < names.size(); i++){
				//if (i != VCFConstants.NormalSample)
					pw.write("\t"+entry.getLAF(i));
			}
			pw.write("\n");
		}
	}
*/
	
	public String getName(int i){
		return names.get(i);
	}
	
	public int getNumofSamples(){
		return names.size();
	}

	/**
	 * Function: getSortedEntriesByGATK(String inputCode, String destCode)
	 * Usage: ArrayList<VCFEntry> entries = db.getSortedEntriesByGATK(inputCode, destCode)
	 * ----
	 * Returns all VCFEntries with matching GATK codes to inputCode in an ArrayList. 
	 * The entries are sorted by how likely they are to be converted to destCode in 
	 * ascending order.
	 * 
	 * @param inputCode	The GATK code that should match valid entries
	 * @param destCode	The GATK code to convert inputCode into
	 * @param pw 
	 * @return	A sorted ArrayList of VCFEntries
	 */
	public ArrayList<VCFEntry> getSortedEntriesByGATK(String inputCode, String destCode, PrintWriter pw){
		ArrayList<SNVEntry> entryList = getEntriesByGATK(inputCode);
		//index -> if 0 to 1
		Map<Integer, Boolean> codeDifferenceMap = initCodeDiffMap(inputCode, destCode);
		Set<Integer> mismatchIndices = new HashSet<Integer>(codeDifferenceMap.keySet());
		ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap = new ArrayList<MyEntry<Double, ArrayList<VCFEntry>>>();
		for (SNVEntry entry: entryList){
			double total = -1.0;
			for (Integer index: mismatchIndices){
				//boolean is0to1 = codeDifferenceMap.get(index);
				double currProb = ((VCFEntry)entry).getSumProb(index);
				//if (codeDifferenceMap.get(index)) currProb = 1.0 - currProb;
				total = (total < 0.0 ? currProb : total * currProb);
			}
			if (total < 0.0) total = 0.0;
			
			if (entryListContainsKey(probsToEntryMap, total)) {
				ArrayList<VCFEntry> currEntries = entryListGet(probsToEntryMap, total);
				currEntries.add((VCFEntry)entry);
				entryListPut(probsToEntryMap, total, currEntries);
				//probsToEntryMap.put(total, entry);
			} else {
				ArrayList<VCFEntry> newEntries = new ArrayList<VCFEntry>();
				newEntries.add((VCFEntry)entry);
				entryListPut(probsToEntryMap, total, newEntries);
			}
		}
		ArrayList<Double> probs = new ArrayList<Double>(new HashSet<Double>(entryListKeySet(probsToEntryMap)));
		Collections.sort(probs);
		//String divider = "-----\n";
		//String outputStr = "Input: " + inputCode + " Dest: " + destCode + "\n";
		//System.out.println(outputStr);
		//pw.write(divider);
		//pw.write(outputStr);
		//pw.write(divider);
		ArrayList<VCFEntry> sortedList = new ArrayList<VCFEntry>();
		for (int i = 0; i < probs.size(); i++){
			//sortedList.add(probsToEntryMap.get(probs.get(i)));
			ArrayList<VCFEntry> newEntryList = entryListGet(probsToEntryMap, probs.get(i));
			if (newEntryList.size() > 1){
				newEntryList.size();
			}
			for(int j = 0; j < newEntryList.size(); j++){
				sortedList.add(newEntryList.get(j));
				//System.out.println(newEntryList.get(j).toString());
			}
//			if (i >= 25){
//				System.out.println("Entries: " + i);
//			}
			//sortedList.addAll();
			//String probLine = i+1 + ". " + "Prob: " + probs.get(i) + "\n";
			//pw.write(probLine);
			//System.out.println(probLine);
			//VCFEntry entry = probsToEntryMap.get(probs.get(i));
			//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + "\n";
			//pw.write(entryLine);
			//System.out.println(entryLine);
		}
		//pw.close();
		//System.out.println(sortedList.toString());
		return sortedList;
	}
	
	private ArrayList<Double> entryListKeySet(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap) {
		// TODO Auto-generated method stub
		ArrayList<Double> newKeyList = new ArrayList<Double>();
		for (MyEntry<Double, ArrayList<VCFEntry>> entry: probsToEntryMap){
			newKeyList.add(entry.getKey());
		}
		return newKeyList;
	}

	private void entryListPut(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap,
			double total, ArrayList<VCFEntry> currEntries) {
		MyEntry<Double, ArrayList<VCFEntry>> newEntry = new MyEntry<Double, ArrayList<VCFEntry>>(total, currEntries);
		probsToEntryMap.add(newEntry);
		// TODO Auto-generated method stub
		
	}

	private ArrayList<VCFEntry> entryListGet(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap,
			double total) {
		// TODO Auto-generated method stub
		ArrayList<VCFEntry> returnList = new ArrayList<VCFEntry>();
		for (MyEntry<Double, ArrayList<VCFEntry>> currEntry : probsToEntryMap){
			//ArrayList<VCFEntry> currVCFEntries = currEntry.getValue();
			if (currEntry.getKey().doubleValue() == (total)) returnList.addAll(currEntry.getValue());
		}
		return returnList;
	}

	private boolean entryListContainsKey(
			ArrayList<MyEntry<Double, ArrayList<VCFEntry>>> probsToEntryMap,
			double total) {
		for (MyEntry<Double, ArrayList<VCFEntry>> currEntry: probsToEntryMap){
			if (currEntry.getKey().equals(total)) return true;
		}
		return false;
	}

	/**
	 * Function: getValidEntries(String inputCode, String destCode, ArrayList<VCFEntry> validEntries)
	 * Usage: int numValid = db.getValidEntries(String inputCode, String destCode, ArrayList<VCFEntry> validEntries)
	 * ----
	 * Finds all valid entries which have the GATK code equivalent to inputCode and can be converted to
	 * destCode. These valid entries are stored in validEntries, which should be an empty ArrayList that is
	 * passed into the function. The function returns the number of valid entries stored in validEntries.
	 * 
	 * @param inputCode	The GATK code to convert
	 * @param destCode	The GATK code to convert to
	 * @param editDistance 
	 * @param validEntries	The ArrayList of entries in which valid entries are stored
	 * @param failCodes 
	 * @param totalMut 
	 * @param pw 
	 * @return	The number of validEntries found as an int
	 */
	public int getValidEntries(String inputCode, String destCode, ArrayList<SNVEntry> validEntries, 
			ArrayList<SNVEntry> failCodes){
		//ArrayList<VCFEntry> entries = getSortedEntriesByGATK(inputCode, destCode, pw);
		ArrayList<SNVEntry> entries = TAG2SNVs.get(inputCode);
		System.out.println("Edit code " + inputCode + " to dest " + destCode);//+", for " + entries.size() + " entries.");
		//For each mismatch in an entry 
		//boolean is0to1 = checkIf0to1(inputCode, destCode);
		Map<Integer, Boolean> mismatchMap = new HashMap<Integer, Boolean>();
		for (int j = 0; j < inputCode.length(); j++){
			if (inputCode.charAt(j) != destCode.charAt(j)) {
				if (inputCode.charAt(j) == '0') mismatchMap.put(j, true);
				else mismatchMap.put(j, false);
			}
		}
		boolean isAll0to1 = true;
		for (Boolean is0to1 : mismatchMap.values()) if (!is0to1) isAll0to1 = false;
		int counter = 0;
		
		for (int i = 0; i < entries.size(); i++){
			
			//int sampleIndex = mismatchIndex(inputCode, destCode);
			if (isAll0to1){
				//check each probability. If all pass, entries are valid
				boolean canConvertAll = true;
				for (Integer indexKey : mismatchMap.keySet()){
					if(!(entries.get(i)).EvidenceOfPresence(indexKey)){
						canConvertAll = false;
						failCodes.add(entries.get(i));
						break;
					}
				}
				if (canConvertAll){
					counter++;
					validEntries.add(entries.get(i));
				}
			} else {
				//for (int j = entries.size() - 1; j >= 0; j--){
				failCodes.add(entries.get(i));
				//}
			}
		}
//			//If all 0 -> 1 keep, else fail
//			double entryProb = entries.get(i).getSumProb(sampleIndex);
//			if (entryProb < THRESHOLD){
//				counter++;
//				VCFEntry entry = entries.get(i);
//				validEntries.add(entry);
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "T\n";
//				//System.out.println(entryLine);
//				//Check these
//				//pw.write(entryLine);
//			} else {
//				VCFEntry entry = entries.get(i);
//				failCodes.add(entry);
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + inputCode + "\n";
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "F\n";
//				//System.out.println(entryLine);
//				//pw.write(entryLine);
//				//break;
//			}
//		}
//		if (is0to1){
//			//pw.write("--Testing Codes--\n");
//			//ArrayList<VCFEntry> currFailCodes = new ArrayList<VCFEntry>();
//			for (int i = 0; i < entries.size(); i++){
//				//Note this currently only works for edit distance == 1
//				int sampleIndex = mismatchIndex(inputCode, destCode);
//				double entryProb = entries.get(i).getSumProb(sampleIndex);
//				if (entryProb < THRESHOLD){
//					counter++;
//					VCFEntry entry = entries.get(i);
//					validEntries.add(entry);
//					//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "T\n";
//					//System.out.println(entryLine);
//					//Check these
//					//pw.write(entryLine);
//				} else {
//					VCFEntry entry = entries.get(i);
//					failCodes.add(entry);
//					//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + inputCode + "\n";
//					//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + destCode + "F\n";
//					//System.out.println(entryLine);
//					//pw.write(entryLine);
//					//break;
//				}
//			}
//			//if 
//			//pw.write("--Testing Ended--\n");
//		} else {
//			for (int i = entries.size() - 1; i >= 0; i--){
//				VCFEntry entry = entries.get(i);
//				failCodes.add(entry);
//				//String entryLine = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + inputCode + "\t" + inputCode + "\n";
//				//pw.write(entryLine);
////				int sampleIndex = mismatchIndex(inputCode, destCode);
////				double entryProb = entries.get(i).getSumProb(sampleIndex);
////				if (entryProb >= THRESHOLD){
////					counter++;
////					validEntries.add(entries.get(i));
////				} else break;
//			}
//		}
		return counter;
	}
	
	
	/**
	 * Function: initCodeDiffMap(String inputCode, String destCode)
	 * Usage: Map<Integer, Boolean> diffMap = initCodeDiffMap(inputCode, destCode)
	 * ----
	 * Generates a map of mismatched indices between the two input parameters and
	 * whether the mismatch is from 0 to 1. In other words, every key is an index
	 * where the characters between the strings do not match. Every value to a key
	 * is a boolean that tells whether the original string's mismatched character 
	 * was a 0.
	 * @param inputCode	The original GATK code
	 * @param destCode	The GATK code to convert to
	 * @return	The map with the mismatched indices to booleans 
	 */
	private Map<Integer, Boolean> initCodeDiffMap(String inputCode,
			String destCode) {
		Map<Integer, Boolean> result = new HashMap<Integer, Boolean>();
		for (int i = 0; i < inputCode.length(); i++){
			char inputChar = inputCode.charAt(i);
			char destChar = destCode.charAt(i);
			if (inputChar != destChar) {
				if (inputChar == '0') result.put(i, true);
				else result.put(i, false);
			}
		}
		return result;
	}

	
	/**
	 * Function: generateMatrix(String outputFile)
	 * Usage: generateMatrix(outputFile)
	 * ----
	 * Generates the matrix of GATK codes and mutation numbers needed
	 * to build a phylogenetic tree. The matrix is written to 
	 * the file specified by the path outputFile.
	 * @param outputFile	The file which to write the matrix
	 */
	private void generateMap(){
		TAG2SNVs = new HashMap<String, ArrayList<SNVEntry>>();
		TAG2RobustSNVNum = new HashMap<String, Integer>();
		
		for (SNVEntry entry: somaticSNVs){
			String code = entry.getGroup();
			if (!TAG2SNVs.containsKey(code)){
				TAG2SNVs.put(code, new ArrayList<SNVEntry>());				
			}
			TAG2SNVs.get(code).add(entry); 
			if (entry.isRobust()){
				if (!TAG2RobustSNVNum.containsKey(code)){
					TAG2RobustSNVNum.put(code, new Integer(1));				
				}else
					TAG2RobustSNVNum.put(code, new Integer(TAG2RobustSNVNum.get(code).intValue()+1));	
			}
		}
		//Quick-fix to remove the entry for all 1's
		/*String imposCode = "";
		for (int i = 0; i < codeLength; i++){
			imposCode += "1";
		}
		
		if (GATKCounter.containsKey(imposCode)) GATKCounter.remove(imposCode);
		*/						
			
		
	}


	public void generateGATKFile(String outputFile) {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new FileWriter(outputFile));
		} catch (IOException e){
			System.out.println("IOEXCEPTION OCCURRED!");
		}
		String header = "Chromosome\tLocation\tGATK\n";
		pw.print(header);
		for (int i = 0; i < somaticSNVs.size(); i++){
			SNVEntry entry = somaticSNVs.get(i);
			String outputStr = entry.getChromosome() + "\t" + entry.getPosition() + "\t" + entry.getGroup() + "\n";
			pw.print(outputStr);
		}
		pw.close();
	}
	
	public void reportGroups() {
		ArrayList<String> codes = new ArrayList<String>(TAG2SNVs.keySet());
		Collections.sort(codes);
		Collections.reverse(codes);
		
		for (int i = 0; i < codes.size(); i++){
			/*if (TAG2SNVs.get(codes.get(i)).size() == 0){
				continue;
			}*/
			System.out.println(codes.get(i) + ": " + TAG2SNVs.get(codes.get(i)).size());
			
		}	
	}
	
	public void printSNVs(String fileName) {
		ArrayList<String> codes = new ArrayList<String>(TAG2SNVs.keySet());
		Collections.sort(codes);
		Collections.reverse(codes);
		PrintWriter pwv = null;//, pwh=null;
		try{
			pwv = new PrintWriter(new FileWriter(fileName));
			pwv.write(getHeader()+"\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
		for (int i = 0; i < codes.size(); i++){
			if (TAG2SNVs.get(codes.get(i)).size() == 0)
				continue;
			System.out.println(codes.get(i) + ": " + TAG2SNVs.get(codes.get(i)).size());
			try {
				printEntriesByGATK(pwv, codes.get(i));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}	
		pwv.close();
	}
	
	
	
	/******************************* Beginning of editing functions **********************************/	
	
	
	/**
	 * Function: editSNV(Set<ArrayList<Integer>> conflicts, Map<String, Double> mutMap, VCFDatabase vcfDB)
	 * Usage: editSNV(conflicts, mutMap, vcfDB)
	 * ----
	 * @param conflicts	The set of conflicting binary codes in the original music
	 * @param mutMap	The map of codes to number of mutations with that GATK code
	 * @param vcfDB		The VCFDatabse corresponding to the VCF file for this run of matrix building
	 * @param iterCounter Which iteration of editSNV this is
	 * @param testName 
	 */
	public void editSNVs(Set<String> conflicts, Set<String> destCodes) {
		
		System.out.append("CONFLICTS: "+conflicts.toString());
		Map<String, ArrayList<SNVEntry>> conflictToPossMutMap = new HashMap<String, ArrayList<SNVEntry>>();

		
		for (String conflict: conflicts){
			
			ArrayList<SNVEntry> failCodes = null;	
			
			//int mutConverted = 0;
			for (String conflictMatch: destCodes){
				if (getEditDist(conflict, conflictMatch) > VCFConstants.EDIT_DISTANCE || conflict.equals(conflictMatch)) 
					continue;
				ArrayList<SNVEntry> possMutations = new ArrayList<SNVEntry>();
				ArrayList<SNVEntry> currFailCodes = new ArrayList<SNVEntry>();
				//boolean is0to1 = checkIf0to1(conflictStr, conflictMatch);
				int currMutConverted =
					getValidEntries(conflict, conflictMatch, possMutations, currFailCodes);
				if (currMutConverted > 0) System.out.println("For code " + conflict + " and dest " + conflictMatch + " we can convert " +
						currMutConverted + ".");
				if (conflictToPossMutMap.containsKey(conflictMatch)) conflictToPossMutMap.get(conflictMatch).addAll(possMutations);
				else
					conflictToPossMutMap.put(conflictMatch,possMutations);
					
				if (failCodes == null && !currFailCodes.isEmpty()) failCodes = new ArrayList<SNVEntry>(currFailCodes);
				else if (!currFailCodes.isEmpty()) failCodes.retainAll(currFailCodes);
				//mutMap.put(conflictMatch, mutMap.get(conflictMatch) + currMutConverted);
				//mutConverted += currMutConverted;
				//System.out.println("Intermediate Mutation Map");
				//System.out.println(mutMap.toString());
			}
			
			//System.out.println("Total Converted: " + mutConverted);
			/*
			 * 1. Find all possible 1-edit distance changes
			 * 2. See which ones are valid; discard rest
			 * 3. 
			 */

		}
		
		
		//Run set cover algorithm
		ArrayList<SNVEntry> movedEntries = new ArrayList<SNVEntry>();
		for (int i = 0; i < conflictToPossMutMap.size(); i++){
			String conflictMatch = findLargestUncoveredSet(conflictToPossMutMap, movedEntries);
			if (conflictMatch == null) break;
			ArrayList<SNVEntry> matchEntries = conflictToPossMutMap.get(conflictMatch);
			matchEntries.removeAll(movedEntries);
			for (int j = 0; j < matchEntries.size(); j++){
				SNVEntry entry = matchEntries.get(j);
				//pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + conflictStr + "\t" + conflictMatch + "\n");
				TAG2SNVs.get(entry.getGroup()).remove(entry);
				if (TAG2SNVs.get(entry.getGroup()).size() ==0){
					conflicts.remove(entry.getGroup());
					TAG2SNVs.remove(entry.getGroup());
				}
				TAG2SNVs.get(conflictMatch).add(entry);
				entry.updateGroup(conflictMatch);
			}
			
			System.out.println("We converted " + matchEntries.size()+ " to GATK code [" + conflictMatch + "].");
			conflictToPossMutMap.remove(conflictMatch);
			movedEntries.addAll(matchEntries);
		}
		
		

}
	

	
	/**
	 * Function: editSNV(Set<ArrayList<Integer>> conflicts, Map<String, Double> mutMap, VCFDatabase vcfDB)
	 * Usage: editSNV(conflicts, mutMap, vcfDB)
	 * ----
	 * @param conflicts	The set of conflicting binary codes in the original music
	 * @param mutMap	The map of codes to number of mutations with that GATK code
	 * @param vcfDB		The VCFDatabse corresponding to the VCF file for this run of matrix building
	 * @param iterCounter Which iteration of editSNV this is
	 * @param testName 
	 */
	public void editSNV(Set<ArrayList<Integer>> conflicts) {
		System.out.println("Original Mutation Map");
		
		String alls = "";
		for (int i = 0; i < names.size(); i++) alls += "1";
		TAG2SNVs.put(alls, new ArrayList<SNVEntry>());
		
		ArrayList<Integer> all0s = new ArrayList<Integer>();
		for (int i = 0; i < names.size(); i++) all0s.add(new Integer(0));
		conflicts.add(all0s);
		
		
		reportGroups();
		
		ArrayList<String> codes = new ArrayList<String>(TAG2SNVs.keySet());
		Collections.sort(codes);
		Collections.reverse(codes);

		//System.out.append("CONFLICTS: "+conflicts.toString());
		
		for (ArrayList<Integer> conflict: conflicts){
			String conflictStr = "";
			for (int i = 0; i < conflict.size(); i++){
				conflictStr += ("" + conflict.get(i));
			}
			
			ArrayList<SNVEntry> failCodes = null;
			Set<String> allCodes = new HashSet<String>(TAG2SNVs.keySet());
			//Remove all conflicts
			for (ArrayList<Integer> currConflict: conflicts){
				String currConflictStr = "";
				for (int i = 0; i < currConflict.size(); i++){
					currConflictStr += ("" + currConflict.get(i));
				}
				allCodes.remove(currConflictStr);
			}

			//System.out.println(editDistMap.toString());
			allCodes = filterEditDistance(conflictStr, allCodes);
			//get edit distances between conflict str and all codes
			//get all possible by iterating through all int from 1 to desired dist
			Set<String> conflictMatches = allCodes;
			System.out.println("For code "+conflictStr+" POSSIBLE target CODES");
			System.out.println(allCodes.toString());
			
			
			ArrayList<String> conflictMatchesList = new ArrayList<String>(conflictMatches);
			//int mutConverted = 0;
			Map<String, ArrayList<SNVEntry>> conflictToPossMutMap = new HashMap<String, ArrayList<SNVEntry>>();
			for (String conflictMatch: conflictMatchesList){
				ArrayList<SNVEntry> possMutations = new ArrayList<SNVEntry>();
				ArrayList<SNVEntry> currFailCodes = new ArrayList<SNVEntry>();
				//boolean is0to1 = checkIf0to1(conflictStr, conflictMatch);
				int currMutConverted =
					getValidEntries(conflictStr, conflictMatch, possMutations, currFailCodes);
				if (currMutConverted > 0) System.out.println("For code " + conflictStr + " and dest " + conflictMatch + " we can convert " +
						currMutConverted + ".");
				conflictToPossMutMap.put(conflictMatch, possMutations);
				if (failCodes == null && !currFailCodes.isEmpty()) failCodes = new ArrayList<SNVEntry>(currFailCodes);
				else if (!currFailCodes.isEmpty()) failCodes.retainAll(currFailCodes);
				//mutMap.put(conflictMatch, mutMap.get(conflictMatch) + currMutConverted);
				//mutConverted += currMutConverted;
				//System.out.println("Intermediate Mutation Map");
				//System.out.println(mutMap.toString());
			}
			//Run set cover algorithm
			ArrayList<SNVEntry> movedEntries = new ArrayList<SNVEntry>();
			for (int i = 0; i < conflictMatchesList.size(); i++){
				String conflictMatch = findLargestUncoveredSet(conflictToPossMutMap, movedEntries);
				if (conflictMatch == null) break;
				ArrayList<SNVEntry> matchEntries = conflictToPossMutMap.get(conflictMatch);
				matchEntries.removeAll(movedEntries);
				int numMoved = matchEntries.size();
				for (int j = 0; j < matchEntries.size(); j++){
					SNVEntry entry = matchEntries.get(j);
					//pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + conflictStr + "\t" + conflictMatch + "\n");
					entry.updateGroup(conflictMatch);
					TAG2SNVs.get(conflictStr).remove(entry);
					TAG2SNVs.get(conflictMatch).add(entry);
				}
				
				//System.out.println("We converted " + numMoved + " entries with GATK code [" + conflictStr + "] to GATK code [" + conflictMatch + "].");
				conflictToPossMutMap.remove(conflictMatch);
				movedEntries.addAll(matchEntries);
			}
			//System.out.println("Total Converted: " + mutConverted);
			/*
			 * 1. Find all possible 1-edit distance changes
			 * 2. See which ones are valid; discard rest
			 * 3. 
			 */
			if (failCodes != null){
				failCodes.removeAll(movedEntries);
				for (int j = 0; j < failCodes.size(); j++){
					SNVEntry entry = failCodes.get(j);
					//pw.write("FAILED\n");
					//pw.write(entry.getChromosome() + "\t" + entry.getPosition() + "\t" + conflictStr + "\t" + conflictStr + "\n");
				}
			}
		}
		System.out.println("Final Mutation Map after merging");
		reportGroups();
		

}
	
	/**
	 * Function: findLargestUncovered(Map<String, ArrayList<SNVEntry>> conflictToPossMutMap, ArrayList<SNVEntry> movedEntries)
	 * Usage: String largestConflict = findLargestUncoveredSet(conflictToPossMutMap, movedEntries)
	 * ----
	 * This algorithm finds the set with the most "uncovered" objects (or entries) in this case.
	 * This is the brunt of the set cover algorithm. 
	 * 
	 * @param conflictToPossMutMap	The map of a conflict to all possible codes it could be moved to
	 * @param movedEntries			All entries that have been moved already
	 * @return	The binary code where the maximal number of entries can be moved
	 */
	private  String findLargestUncoveredSet(
			Map<String, ArrayList<SNVEntry>> conflictToPossMutMap,
			ArrayList<SNVEntry> movedEntries) {
		String bestMatch = null;
		int maxEntries = 0;
		ArrayList<String> conflictMatchesList = new ArrayList<String>(conflictToPossMutMap.keySet());
		for (String conflictMatch: conflictMatchesList){
			ArrayList<SNVEntry> matchEntries = conflictToPossMutMap.get(conflictMatch);
			matchEntries.removeAll(movedEntries);
			if (matchEntries.size() !=0 && matchEntries.size()+TAG2SNVs.get(conflictMatch).size() > maxEntries) {
				bestMatch = conflictMatch;
				maxEntries = matchEntries.size()+TAG2SNVs.get(conflictMatch).size();
			}
		}
		return bestMatch;
	}

	private  Set<String> filterEditDistance(String conflictCode, Set<String> allCodes) {
		// TODO Auto-generated method stub
		Set<String> resultCodes = new HashSet<String>();
		for (String code : allCodes){
			if (getEditDist(conflictCode, code) <= VCFConstants.EDIT_DISTANCE && !conflictCode.equals(code))
				resultCodes.add(code);
		}
		return resultCodes;
	}

	private  int getEditDist(String conflictStr, String code) {
		// TODO Auto-generated method stub
		int dist = 0;
		for (int i = 0; i < conflictStr.length(); i++)
			if (conflictStr.charAt(i) != code.charAt(i)) dist++;
		return dist;
	}

	public int getNumRobustSNVs(String group){
		return TAG2RobustSNVNum.get(group);
	}
/******************************* END of editing functions **********************************/	
	
	
	
	
	//Taken from StackOverflow
	//http://stackoverflow.com/questions/3110547/java-how-to-create-new-entry-key-value
	final class MyEntry<K, V> implements Map.Entry<K, V> {
	    private final K key;
	    private V value;

	    public MyEntry(K key, V value) {
	        this.key = key;
	        this.value = value;
	    }

	    @Override
	    public K getKey() {
	        return key;
	    }

	    @Override
	    public V getValue() {
	        return value;
	    }

	    @Override
	    public V setValue(V value) {
	        V old = this.value;
	        this.value = value;
	        return old;
	    }
	}
	
	

}
