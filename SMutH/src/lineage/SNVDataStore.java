package lineage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import util.CNVRegion;
import util.Configs;
import util.MUTEntry;
import util.SNVAnnotation;
import util.SNVEntry;

public class SNVDataStore {
	/** Input samples */
	private int numSamples;
	private ArrayList<String> sampleNames;
	/** List of all input somatic SNVs */
	protected ArrayList<SNVEntry> somaticSNVs;
	/** List of SNVs with unreliable VAFs in at least one sample */
	protected ArrayList<SNVEntry> unreliableSNVs;
	/** Map of sample profile tags to a list of SNVs with this profile */
	private HashMap<String, ArrayList<SNVEntry>> tag2SNVs;
	private int normalSample;

	public SNVDataStore(String snvInputFile, int normalSampleId) {
		normalSample = normalSampleId;
		somaticSNVs = new ArrayList<SNVEntry>();
		tag2SNVs = new HashMap<String, ArrayList<SNVEntry>>();
		unreliableSNVs = new ArrayList<SNVEntry>();
		
		switch (Configs.INFORMAT){		
			case MUT :
				loadMUT(snvInputFile);
				break;
			case VCF :
				//loadVCF(snvInputFile,normalSample);
				break;
			case SIM :
				loadSIM(snvInputFile);
				break;
			case FL:
				break;
			default:
				break;
		}
		
		for (SNVEntry entry: somaticSNVs) {
			String tag = entry.getGroup();
			if(!entry.isRobust()) {
				unreliableSNVs.add(entry);
				continue;
			}
			// group robust SNVs only
			if (!tag2SNVs.containsKey(tag)){
				tag2SNVs.put(tag, new ArrayList<SNVEntry>());	
			}
			tag2SNVs.get(tag).add(entry); 
		}
		reportSNVGroups();
		
		// handle mutations in groups with insufficient support as unreliable
		ArrayList<String> smallGroups = new ArrayList<String>();
		for(String tag : tag2SNVs.keySet()) {
			if(tag2SNVs.get(tag).size() < Configs.ROBUSTGROUP_SIZE_THR) {
				for(SNVEntry entry : tag2SNVs.get(tag)) {
					unreliableSNVs.add(entry);
				}
				smallGroups.add(tag);
			}
		}
		for(String group : smallGroups) {
			tag2SNVs.remove(group);
		}
		
		assignUnreliableSNVs();
		filterGroups();
		reportSNVGroups();
	}
	
	/**
	 * Checks the validity of the SNV VAFs and coverage
	 */
	private boolean isValidSNV(SNVEntry snv) {
		for(int i = 0; i < numSamples; i++) {
			// check if the VAF i too high
			if(snv.getAAF(i) > Configs.MAX_ALLOWED_VAF) {
				return false;
			}
			// check if coverage is too low
		}
		return true;
	}
	
	private void detectSystematicNoise(SNVEntry snv) {
		if(snv.getAAF(normalSample) < Configs.SYST_ERR_THR) {
			return;
		}
		// detected a high level of noise in the normal
		ArrayList<Double> vafs = new ArrayList<Double>();
		for(int i = 0; i < numSamples; i++) {
			if(snv.getAAF(i) < Configs.VALIDATION_THR) {
				vafs.add(snv.getAAF(i));
			}
		}
		
		if(vafs.size() == 0) return;
		
		Collections.sort(vafs);
		double noise = vafs.get(vafs.size()/2);
		System.out.println("High noise in SNV " + snv.getGroup() +" normal sample: " + snv + "\n Subtracting: " + noise);
		
		String group = "";
		for(int i = 0; i < numSamples; i++) {
			if(snv.getAAF(i) >= noise) {
				snv.setAAF(i, (snv.getAAF(i)-noise));
				System.out.print(snv.getAAF(i) + " ");
				if(snv.getAAF(i) >= Configs.VALIDATION_THR) {
					group += '1';
				} else {
					group += '0';
				}
			} else {
				group += '0';
			}
		}
		snv.updateGroup(group);
		System.out.println("Result " + snv.getGroup());
	}
	
	private void filterGroups() {
		// apply minimum size and robust size constraint
		ArrayList<String> filteredOut = new ArrayList<String>();
		for(String tag : tag2SNVs.keySet()) {
			if(tag2SNVs.get(tag).size() < Configs.GROUP_SIZE_THR) {
				filteredOut.add(tag);
				continue;
			}
			int numRobust = 0;
			for(SNVEntry entry : tag2SNVs.get(tag)) {
				if(entry.isRobust()) {
					numRobust++;
					if(numRobust >= Configs.GROUP_ROBUST_NUM_THR) {
						break;
					}
				}
			}
			if(numRobust < Configs.GROUP_ROBUST_NUM_THR) {
				filteredOut.add(tag);
				continue;
			}
		}
		for(String group : filteredOut) {
			tag2SNVs.remove(group);
		}
	}
	
	private void assignUnreliableSNVs() {
		if (unreliableSNVs.size() == 0) return;
		
		ArrayList<String> targetTags = new ArrayList<String>(tag2SNVs.keySet());
		String all1s = "", all0s = "";
		for (int i = 0; i < numSamples; i++) { 
			all1s += "1";
			all0s += "0";
		}
		if(!tag2SNVs.containsKey(all0s)) {
			targetTags.add(all0s);
		}
		targetTags.add(all1s);

		System.out.println("Ambiguous profile SNV assignment (" + unreliableSNVs.size() + " total): ");
		ArrayList<SNVEntry> toRemove = new ArrayList<SNVEntry>();
		for(SNVEntry snv : unreliableSNVs) {
			String tag = snv.getGroup();	
			double bestDistToTarget = 0;
			String bestTarget = "";
				
			for(String target : targetTags) {				
				if(getHammingDist(tag, target) > Configs.EDIT_DISTANCE || !canConvert(snv, target)) continue;
				// if the target is germline, move to germline regardless of distance
				if(target.equals(all1s)) {
					toRemove.add(snv);
					System.out.println("**Removed as germline: " + snv);
					break;
				}				
				if(target.equals(all0s)) {
					double dist = distToZero(snv);
					if(dist > bestDistToTarget) {
						bestDistToTarget = dist;
						bestTarget = target;
					}
					continue;
				}
				// compare to each mutation in the target group
				for(SNVEntry targetSNV : tag2SNVs.get(target)) {
					// compute the distance to the target snv
					double dist = distToTargetVAF(snv, targetSNV);
					if(dist > bestDistToTarget) {
						bestDistToTarget = dist;
						bestTarget = target;
					}
				}
			}
			if(bestTarget.equals(all0s)) {
				snv.updateGroup(bestTarget);
				System.out.println("Assigned " + tag + " to " + bestTarget + " with dist " + bestDistToTarget + " (" + snv + ")");
				continue;
			}
			
			if(bestDistToTarget != 0 && bestDistToTarget >= getHammingWeight(bestTarget)*Configs.MIN_VAF_TARGET_RATIO_PER_SAMPLE) {
				// found a valid match
				snv.updateGroup(bestTarget);
				toRemove.add(snv);
				//System.out.println(bestTarget);
				tag2SNVs.get(bestTarget).add(snv);
				System.out.println("Assigned " + tag + " to " + bestTarget + " with dist " + bestDistToTarget + " (" + snv + ")");
				continue;
			}
			
			// if there is a matching robust group tag, remove the snv
			//if((!tag.equals(all0s)) && tag2SNVs.containsKey(tag)) {
				//toRemove.add(snv);
				//System.out.println("**Removed " + tag + " as conflict" + ((bestDistToTarget != 0) ? " (best dist = " + bestDistToTarget + " to target " + bestTarget + "): " : " ") + snv);
			//} else {
				System.out.println("**Could not assign " + tag + ", no conflict" + ((bestDistToTarget != 0) ? " (best dist = " + bestDistToTarget + " to target " + bestTarget + "): " : " ") + snv);
			//}
		}
		unreliableSNVs.removeAll(toRemove);
		
		// remaining snvs had no suitable matches and potentially represent true branches
		// we minimize the number of additional nodes by merging the groups
		HashMap<String, ArrayList<SNVEntry>> unreliableGroups = mergeUnreliableSNVs(unreliableSNVs);
		
		for(String tag : unreliableGroups.keySet()) {
			if(tag.equals(all0s)) continue;
			//if(!allowGroup(snv)) continue;	
			if(!tag2SNVs.containsKey(tag)) {
				tag2SNVs.put(tag, unreliableGroups.get(tag));
				System.out.println("New group: " + tag);
			}
		}
		tag2SNVs.remove(all0s);
	}
	
	/*private boolean allowGroup(SNVEntry snv) {
		String tag = snv.getGroup();
		// if there are ambiguous 1s, don't allow it
		for(int i = 0; i < tag.length(); i++) {
			if(snv.getAAF(i) < (Configs.VALIDATION_THR + 0.01) && snv.getAAF(i) > (Configs.VALIDATION_THR - 0.01)) {
				return false;
			}
		}
		return true;
	}*/
	
	private HashMap<String, ArrayList<SNVEntry>> mergeUnreliableSNVs(ArrayList<SNVEntry> snvs) {
		HashMap<String, ArrayList<SNVEntry>> groups = new HashMap<String, ArrayList<SNVEntry>>();
		HashMap<String, ArrayList<SNVEntry>> adj = new HashMap<String, ArrayList<SNVEntry>>();
		
		// generate all possible target profiles 
		ArrayList<String> targets = generateAllPossibleTargets(snvs);
		for(String t : targets) {
			//System.out.println("Target " + t);
			if(!adj.keySet().contains(t)) {
				adj.put(t, new ArrayList<SNVEntry>());
			}
		}
		
		// min vertex cover
		for(SNVEntry entry : snvs) {
			for(String target : adj.keySet()) {
				if(canConvert(entry, target)) {
					//System.out.println("Can convert " + entry.getGroup() + " to " + target);
					adj.get(target).add(entry);
				}
			}
		}
		while(snvs.size() > 0) {
			// find the largest set
			int maxSize = 0;
			String maxSet = "";
			for(String target : adj.keySet()) {
				if(adj.get(target).size() > maxSize) {
					maxSize = adj.get(target).size();
					maxSet = target;
				}
			}
			System.out.println("Max size = " + maxSize);
			if(maxSize == 1) break;
			
			groups.put(maxSet, adj.get(maxSet));
			for(SNVEntry entry : adj.get(maxSet)) {
				System.out.println("Converted " + entry.getGroup() + " to " + maxSet);
				entry.updateGroup(maxSet);
				for(String target : adj.keySet()) {
					if(target.equals(maxSet)) continue;
					ArrayList<SNVEntry> l = adj.get(target);
					if(l.contains(entry)) {
						l.remove(entry);
					}
				}
				snvs.remove(entry);
			}
			adj.remove(maxSet);
		}
		
		// the remaining targets are supported by only 1 SNV
		// decide the groups of the remaining SNVs based on their proximity to the thresholds
		for(SNVEntry snv : snvs) {
			String tag = "";
			for(int i = 0; i < numSamples; i++) {
				if(snv.getGroup().charAt(i) == '0' && snv.EvidenceOfPresence(i)) {
					double delta0 = snv.getAAF(i) - Configs.VALIDATION_SOFT_THR;
					double delta1 = Configs.VALIDATION_THR - snv.getAAF(i);
					tag += (delta1 < delta0) ? '1' : snv.getGroup().charAt(i);
				} else {
					tag += snv.getGroup().charAt(i);
				}
			}
			System.out.println("Converted " + snv.getGroup() + " to " + tag);
			if(!groups.containsKey(tag)) {
				groups.put(tag, new ArrayList<SNVEntry>());
			}
			groups.get(tag).add(snv);
		}
		
		return groups;
	}
	
	private ArrayList<String> generateAllPossibleTargets(ArrayList<SNVEntry> snvs) {
		ArrayList<String> targets = new ArrayList<String>();
		for(SNVEntry snv: snvs) {
			extendTarget("", snv, targets);
		}
		return targets;
	}
	
	private void extendTarget(String partialProfile, SNVEntry snv, ArrayList<String> targets) {
		int sample = partialProfile.length();
		if(sample == numSamples) {
			//System.out.println(snv.getGroup() + " to " + partialProfile);
			targets.add(partialProfile);
			return;
		}
		extendTarget(partialProfile + ""+ snv.getGroup().charAt(sample), snv, targets);
		if(snv.getGroup().charAt(sample) == '0' && snv.EvidenceOfPresence(sample)) {
			extendTarget(partialProfile + '1', snv, targets);
		}
	}
	
	private double distToZero(SNVEntry snv) {
		String tag = snv.getGroup();
		double dist = 0;
		for(int i = 0; i < numSamples; i++) {
			if(tag.charAt(i) == '1' || snv.EvidenceOfPresence(i)) {
				dist += 0.0001/snv.getAAF(i);
			}
		}
		return dist;
	}
	
	private double distToTargetVAF(SNVEntry snv, SNVEntry targetSNV) {
		String tag = snv.getGroup();
		String targetTag = targetSNV.getGroup();
		double dist = 0;
		for(int i = 0; i < numSamples; i++) {
			if(targetTag.charAt(i) == '0' && tag.charAt(i) == '0' && !snv.EvidenceOfPresence(i)) continue;
			double min = targetSNV.getAAF(i) < snv.getAAF(i) ? targetSNV.getAAF(i) : snv.getAAF(i);
			double max = targetSNV.getAAF(i) > snv.getAAF(i) ? targetSNV.getAAF(i) : snv.getAAF(i);
			if(min == 0) {
				min = 0.0001;
			}
			dist += min/max;
		}
		return dist;
	}
	
	private boolean canConvert(SNVEntry snv, String target) {
		String tag = snv.getGroup();		
		for (int i = 0; i < tag.length(); i++) {
			if (tag.charAt(i) != target.charAt(i)) {
				if(tag.charAt(i) == '1') {
					return false;
				}
				if(!snv.EvidenceOfPresence(i)) {
					return false;
				}
			}
		}
		return true;
	}
	
	private int getHammingDist(String s1, String s2) {
		int dist = 0;
		for (int i = 0; i < s1.length(); i++){
			if (s1.charAt(i) != s2.charAt(i)) dist++;
		}
		return dist;
	}
	
	private int getHammingWeight(String tag) {
		int w = 0;
		for(int i = 0; i < tag.length(); i++) {
			if(tag.charAt(i) == '1') {
				w++;
			}
		}
		return w;
	}
	
	/****** Accessors *****/
	
	public int getNumSamples() {
		return numSamples;
	}
	
	public String getSampleName(int i) {
		return sampleNames.get(i);
	}
	
	public HashMap<String, ArrayList<SNVEntry>> getSomaticGroups() {
		return tag2SNVs;
	}
	
	
	/** A group is robust if it contains at least 1 robust mutation */
	public boolean isRobustGroup(String groupTag) {
		int numRobust = 0;
		for(SNVEntry entry : tag2SNVs.get(groupTag)) {
			if(entry.isRobust()) {
				numRobust++;
				if(numRobust >= Configs.ROBUSTGROUP_SIZE_THR) {
					return true;
				}
			}
		}
		return false;
	}
	
	/****** IO ******/
	
	public void reportSNVGroups() {
		ArrayList<String> tags = new ArrayList<String>(tag2SNVs.keySet());
		Collections.sort(tags);
		Collections.reverse(tags);
		
		System.out.println("SNV Map:");
		System.out.println("Group\t#SNVs\t#Robust\tMean VAF(Robust Mean VAF)");	
		for (String tag : tags){
			if (tag2SNVs.get(tag).size() == 0) {
				continue;
			}
			
			double[] mean = new double[numSamples];	
			double[] robustMean = new double[numSamples];	
			int robustCount = 0;
			for (SNVEntry entry : tag2SNVs.get(tag)) {
				for (int i = 0; i < numSamples; i++){
					mean[i] += entry.getAAF(i);
					if(entry.isRobust()) {
						robustMean[i] += entry.getAAF(i);
					}
				}
				if(entry.isRobust()) {
					robustCount++;
				}
			}
			for (int i = 0; i < numSamples; i++) {
				mean[i] = mean[i]/(tag2SNVs.get(tag).size());
				if(robustCount > 0) {
					robustMean[i] = robustMean[i]/robustCount;
				}
			}
			
			System.out.print(tag + "\t" + tag2SNVs.get(tag).size()+ "\t" + robustCount + "\t");
			for (int i = 0; i < numSamples; i++) {
				if(robustCount != 0) {
					System.out.printf(" %.3f(%.3f)", mean[i], robustMean[i]);
				} else {
					System.out.printf(" %.3f(%.3f)", mean[i], 0f);
				}
			}
			System.out.println();
		}	
	}
	
	private void setSampleNames(String headerLine) {	
		String[] header = headerLine.split("\t");
		switch (Configs.INFORMAT){
			case VCF :
				sampleNames = new ArrayList<String>(Arrays.asList(header).subList(9, header.length));
				break;
			case MUT :
				sampleNames = new ArrayList<String>(Arrays.asList(header).subList(5, header.length));
				break;
			case FL :
				sampleNames = new ArrayList<String>(Arrays.asList(header).subList(4, header.length));
				break;
			case SIM :
				sampleNames = new ArrayList<String>(Arrays.asList(header).subList(1, header.length));
				break;			
		}
	}
	
	private void loadMUT(String inputFile) {
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String lastLine = "";
			String currLine = rd.readLine();
			while (currLine.substring(0, 1).equals("#")) { 
				lastLine = currLine; 
				currLine = rd.readLine();
			}
			setSampleNames(lastLine);
			numSamples = sampleNames.size();
			System.out.println("There are " + numSamples + " samples!");
			int totalSNVCounter = 0;
			int germlineCounter = 0;
			System.out.println("SNV Entries:");
			while (currLine != null){
				MUTEntry entry = new MUTEntry(currLine, numSamples);
				currLine = rd.readLine();
				totalSNVCounter++;
				
				if(entry.getGenotype(normalSample).equals("1/1")) {
					System.out.println("**Filtered as germline: " + entry);
					germlineCounter++;
					continue;
				}
				if(!isValidSNV(entry)) {
					System.out.println("**Filtered as untrustable: " + entry);
					continue;
				}
				System.out.println(entry + "\t" + entry.getGroup() +"\t" + entry.isRobust());
				detectSystematicNoise(entry);
				somaticSNVs.add(entry);	
			}
			rd.close();
			System.out.println("There are " + totalSNVCounter + " SNVs in validation file. Of those, we pass "+ germlineCounter + " as germline, and "+ somaticSNVs.size() +" as somatic. \n");
		} catch (IOException e){
			System.out.println("File Reading Error!");
		}
		
	}
	
	private void loadSIM(String inputFile){
		System.out.println(normalSample);
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String currLine = rd.readLine();
			setSampleNames(currLine);
			numSamples = sampleNames.size();
			int totalSNVCounter = 0;
			int germlineCounter = 0;
			
			currLine = rd.readLine();
			while (currLine != null){
				MUTEntry entry = new MUTEntry(currLine, numSamples);
				currLine = rd.readLine();	
				
				// filter out germline mutations
				if(entry.getGenotype(normalSample).equals("1/1")) {
					System.out.println("**Filtered as germline: " + entry);
					germlineCounter++;
					continue;
				} 				
				System.out.println(entry + "\t" + entry.getGroup() +"\t" + entry.isRobust());
				somaticSNVs.add(entry);	
			}
			rd.close();
			System.out.println("There are " + totalSNVCounter + " SNVs in validation file. Of those, we pass "+ germlineCounter + " as germline, and "+ somaticSNVs.size() +" as somatic. \n");
		} catch (IOException e){
			System.out.println("File Reading Error!");
		}
	}
	
	private ArrayList<CNVRegion> loadCNVs(String inputCNVFile){
		ArrayList<CNVRegion> CNVs = new ArrayList<CNVRegion>();
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputCNVFile));
			String currLine = rd.readLine();
			if (currLine.charAt(0) == '#') {
				currLine = rd.readLine();             
			}
			while (currLine != null){
				CNVRegion c= new CNVRegion(currLine);
				CNVs.add(c);
				currLine = rd.readLine();
			}
			rd.close();
		} catch (IOException e) {
			System.out.println("CNV input file Reading Error!");
		}
		return CNVs;
	}
	
	private ArrayList<SNVAnnotation> loadAnnovarFunction(String inputAnnFile){
		ArrayList<SNVAnnotation> anns = new ArrayList<SNVAnnotation>();
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputAnnFile));
			String currLine = rd.readLine();
			while (currLine != null){
				SNVAnnotation ann = new SNVAnnotation();
				String[] entryParts = currLine.split("\t");
				ann.codingInfo = entryParts[0];
				ann.geneInfo = entryParts[1];
				anns.add(ann);
				currLine = rd.readLine();
			}
			rd.close();
		} catch (IOException e) {
			System.out.println("ANN input file Reading Error!");
		}
		return anns;
	}
	
	// returns COSMIC info by SNV position
	private HashMap<Integer, String> loadCOSMIC(String inputCOSMICFile){
		HashMap<Integer, String> cosmicDB = new HashMap<Integer, String>();
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputCOSMICFile));
			String currLine = rd.readLine();
			while (currLine != null){
				String[] entryParts = currLine.split("\t");
				cosmicDB.put(Integer.parseInt(entryParts[3]), entryParts[1]);
				currLine = rd.readLine();
			}
			rd.close();
		} catch (IOException e) {
			System.out.println("COSMIC input file Reading Error!");
		}
		return cosmicDB;
	}
	
	private HashMap<Integer, ArrayList<Integer>> loadTCGA(String inputTCGAFile){
		HashMap<Integer, ArrayList<Integer>> tcgaDB = new HashMap<Integer, ArrayList<Integer>>();
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputTCGAFile));
			String currLine = rd.readLine();
			while (currLine != null){
				String[] entryParts = currLine.split("\t");
				int chr = Integer.parseInt(entryParts[0].trim());
				int pos = Integer.parseInt(entryParts[1].trim());
				if(tcgaDB.containsKey(chr)) {
					tcgaDB.get(chr).add(pos);
				} else {
					tcgaDB.put(chr, new ArrayList<Integer>());
					tcgaDB.get(chr).add(pos);
				}
				currLine = rd.readLine();
			}
			rd.close();
		} catch (IOException e) {
			System.out.println("TCGA input file Reading Error!");
		}
		return tcgaDB;
	}
	
	public void annotateSNVs(String inputCNVFile, String inputAnnFile, String cosmicFile, String tcgaFile) {
		ArrayList<CNVRegion> cnvs = null;
		if(inputCNVFile != null) {
			cnvs = loadCNVs(inputCNVFile);
		}
		ArrayList<SNVAnnotation> anns = null;
		if(inputAnnFile != null) {
			anns = loadAnnovarFunction(inputAnnFile);
		}
		HashMap<Integer, String> cosmicDB = null;
		if(cosmicFile != null) {
			cosmicDB = loadCOSMIC(cosmicFile);
		}
		HashMap<Integer, ArrayList<Integer>> tcgaDB = null;
		if(cosmicFile != null) {
			tcgaDB = loadTCGA(tcgaFile);
		}
		for (int i = 0; i < somaticSNVs.size(); i++) {
			SNVEntry entry = somaticSNVs.get(i);
			if (cnvs != null){
				entry.checkInCNVRegion(cnvs);
			}
			if(anns != null) {
				entry.setAnnotation(anns.get(i));
			}
			if(cosmicDB != null) {
				if(cosmicDB.containsKey(entry.getPosition())) {
					entry.getAnnotation().cosmic = cosmicDB.get(entry.getPosition());
				}
			}
			if(tcgaDB != null) {
				if(tcgaDB.containsKey(entry.getChromNum()) && tcgaDB.get(entry.getChromNum()).contains(entry.getPosition())) {
					entry.getAnnotation().tcga = "+TCGA+";
				}
			}
		}
	}
}
