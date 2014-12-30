package lineage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import lineage.Parameters.Format;
import util.CNVRegion;

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
	
	private static Logger logger = LineageEngine.logger;

	public SNVDataStore(String snvInputFile, int normalSampleId) {
		normalSample = normalSampleId;
		somaticSNVs = new ArrayList<SNVEntry>();
		tag2SNVs = new HashMap<String, ArrayList<SNVEntry>>();
		unreliableSNVs = new ArrayList<SNVEntry>();
		
		if(Parameters.INFORMAT == Format.MUT || Parameters.INFORMAT == Format.MUTC) {
			loadMUT(snvInputFile);
		} else {
			//loadVCF(snvInputFile,normalSample);
			System.err.println("ERROR: VCF format is not supported");
			System.exit(1);
		}
		
		for (SNVEntry entry: somaticSNVs) {
			String tag = entry.getGroupProfile();
			if(Parameters.INFORMAT != Parameters.Format.MUTC && !entry.isRobust()) {
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
		if(Parameters.INFORMAT != Parameters.Format.MUTC) {
			ArrayList<String> smallGroups = new ArrayList<String>();
			for(String tag : tag2SNVs.keySet()) {
				if(tag2SNVs.get(tag).size() < Parameters.ROBUSTGROUP_SIZE_THR) {
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
	}
	
	/**
	 * Checks the validity of the SNV VAFs and coverage
	 */
	private boolean isValidSNV(SNVEntry snv) {
		for(int i = 0; i < numSamples; i++) {
			// check if the VAF i too high
			if(snv.getAAF(i) > Parameters.MAX_ALLOWED_VAF) {
				return false;
			}
			// check if coverage is too low
			//if(snv.getTotalReadCount(i) < Configs.MIN_ALLOWED_COVERAGE) {
				//System.out.println("LOW coverage " + snv);
				//return false;
			//}
		}
		return true;
	}
	
	private boolean isAbsentSNV(SNVEntry snv) {
		String tag = snv.getGroupProfile();
		for(int i = 0; i < numSamples; i++) {
			if(tag.charAt(i) == '1') {
				return false;
			}
		}
		return true;
	}
	
	private void filterGroups() {
		// apply minimum size and robust size constraint
		ArrayList<String> filteredOut = new ArrayList<String>();
		for(String tag : tag2SNVs.keySet()) {
			if(tag2SNVs.get(tag).size() < Parameters.GROUP_SIZE_THR) {
				filteredOut.add(tag);
				continue;
			}
			boolean zeros = true;
			for(int i = 0; i < numSamples; i++) {
				if(tag.charAt(i) == '1') {
					zeros = false;
					break;
				}
			}
			if(zeros) {
				filteredOut.add(tag);
				continue;
			}
			
			int numRobust = 0;
			for(SNVEntry entry : tag2SNVs.get(tag)) {
				if(entry.isRobust()) {
					numRobust++;
					if(numRobust >= Parameters.GROUP_ROBUST_NUM_THR) {
						break;
					}
				}
			}
			if(numRobust < Parameters.GROUP_ROBUST_NUM_THR) {
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

		logger.log(Level.FINE, "Ambiguous profile SNV assignment (" + unreliableSNVs.size() + " total): ");
		ArrayList<SNVEntry> toRemove = new ArrayList<SNVEntry>();
		for(SNVEntry snv : unreliableSNVs) {
			String tag = snv.getGroupProfile();	
			double bestDistToTarget = 0;
			String bestTarget = "";
				
			for(String target : targetTags) {				
				if(!canConvert(snv, target)) continue;
				// if the target is germline, move to germline regardless of distance
				if(target.equals(all1s)) {
					toRemove.add(snv);
					logger.log(Level.FINE, "**Removed as germline: " + snv);
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
				logger.log(Level.FINE, "Assigned " + tag + " to " + bestTarget + " with dist " + bestDistToTarget + " (" + snv + ")");
				continue;
			}
			
			if(bestDistToTarget != 0 && bestDistToTarget >= getHammingWeight(bestTarget)*Parameters.MIN_VAF_TARGET_RATIO_PER_SAMPLE) {
				// found a valid match
				snv.updateGroup(bestTarget);
				toRemove.add(snv);
				tag2SNVs.get(bestTarget).add(snv);
				logger.log(Level.FINE, "Assigned " + tag + " to " + bestTarget + " with dist " + bestDistToTarget + " (" + snv + ")");
				continue;
			}
			
			logger.log(Level.FINE, "**Could not assign " + tag + ", no conflict" + ((bestDistToTarget != 0) ? " (best dist = " + bestDistToTarget + " to target " + bestTarget + "): " : " ") + snv);

		}
		unreliableSNVs.removeAll(toRemove);
		
		// remaining snvs had no suitable matches and potentially represent true branches
		// we minimize the number of additional nodes by merging the groups
		HashMap<String, ArrayList<SNVEntry>> unreliableGroups = mergeUnreliableSNVs(unreliableSNVs);
		
		for(String tag : unreliableGroups.keySet()) {
			if(tag.equals(all0s)) continue;
			if(!tag2SNVs.containsKey(tag)) {
				tag2SNVs.put(tag, unreliableGroups.get(tag));
				logger.log(Level.FINE, "New group: " + tag);
			}
		}
		tag2SNVs.remove(all0s);
	}
	
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
			if(maxSize == 1) break;
			
			groups.put(maxSet, adj.get(maxSet));
			for(SNVEntry entry : adj.get(maxSet)) {
				logger.log(Level.FINE, "Converted " + entry.getGroupProfile() + " to " + maxSet);
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
				if(snv.getGroupProfile().charAt(i) == '0' && snv.evidenceOfPresence(i)) {
					double delta0 = snv.getAAF(i) - Parameters.VALIDATION_SOFT_THR;
					double delta1 = Parameters.VALIDATION_THR - snv.getAAF(i);
					tag += (delta1 < delta0) ? '1' : snv.getGroupProfile().charAt(i);
				} else {
					tag += snv.getGroupProfile().charAt(i);
				}
			}
			logger.log(Level.FINE, "Converted " + snv.getGroupProfile() + " to " + tag);
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
			targets.add(partialProfile);
			return;
		}
		extendTarget(partialProfile + ""+ snv.getGroupProfile().charAt(sample), snv, targets);
		if(snv.getGroupProfile().charAt(sample) == '0' && snv.evidenceOfPresence(sample)) {
			extendTarget(partialProfile + '1', snv, targets);
		}
	}
	
	private double distToZero(SNVEntry snv) {
		String tag = snv.getGroupProfile();
		double dist = 0;
		for(int i = 0; i < numSamples; i++) {
			if(tag.charAt(i) == '1' || snv.evidenceOfPresence(i)) {
				dist += 0.0001/snv.getAAF(i);
			}
		}
		return dist;
	}
	
	private double distToTargetVAF(SNVEntry snv, SNVEntry targetSNV) {
		String tag = snv.getGroupProfile();
		String targetTag = targetSNV.getGroupProfile();
		double dist = 0;
		for(int i = 0; i < numSamples; i++) {
			if(targetTag.charAt(i) == '0' && tag.charAt(i) == '0' && !snv.evidenceOfPresence(i)) continue;
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
		String tag = snv.getGroupProfile();		
		for (int i = 0; i < tag.length(); i++) {
			if (tag.charAt(i) != target.charAt(i)) {
				if(tag.charAt(i) == '1') {
					return false;
				}
				if(!snv.evidenceOfPresence(i)) {
					return false;
				}
			}
		}
		return true;
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
				if(numRobust >= Parameters.ROBUSTGROUP_SIZE_THR) {
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
		
		logger.log(Level.FINE, "SNV Map:");
		logger.log(Level.FINE, "Group\t#SNVs\t#Robust\tMean VAF(Robust Mean VAF)");	
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
			
			String m = tag + "\t" + tag2SNVs.get(tag).size()+ "\t" + robustCount + "\t";
			DecimalFormat df = new DecimalFormat("#.##");
			for (int i = 0; i < numSamples; i++) {
				if(robustCount != 0) {
					m += df.format(mean[i]) + "(" + df.format(robustMean[i]) + ") \t";
				} else {
					m += df.format(mean[i]) + "(0) \t";
				}
			}
			logger.log(Level.FINE, m);
		}	
	}
	
	private void parseHeader(String headerLine) {	
		String[] headerParts = headerLine.split("\t");
		int numRequiredFields = SNVEntry.NUM_REQ_FIELDS_MUT;
		if(Parameters.INFORMAT == Format.MUTC) {
			numRequiredFields = SNVEntry.NUM_REQ_FIELDS_MUTC;
		}
		if(headerParts.length <= numRequiredFields + 1) {
			System.err.println("ERROR: Wrong input header format: must have " + numRequiredFields + " required fields and at least 2 samples"); 
			System.exit(1);
		}
		sampleNames = new ArrayList<String>(Arrays.asList(headerParts).subList(numRequiredFields, headerParts.length));
	}
	
	private void loadMUT(String inputFile) {
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String currLine = rd.readLine();
			if(currLine == null) {
				System.err.println("ERROR: empty input file"); 
				System.exit(1);
			}
			parseHeader(currLine);
			numSamples = sampleNames.size();
			logger.log(Level.FINE, "There are " + numSamples + " samples!");
			
			currLine = rd.readLine();
			int totalSNVCounter = 0;
			int germlineCounter = 0;
			while (currLine != null){
				SNVEntry entry = new SNVEntry(currLine, numSamples);
				currLine = rd.readLine();
				totalSNVCounter++;
				
				if(entry.getGenotype(normalSample).equals("1/1")) {
					logger.log(Level.INFO, "**Filtered as germline: " + entry);
					germlineCounter++;
					continue;
				}
				if(!isValidSNV(entry)) {
					logger.log(Level.INFO, "**Filtered due to VAFs > allowed MAX: " + entry);
					continue;
				}
				if(entry.isRobust() && isAbsentSNV(entry)) {
					logger.log(Level.INFO, "**Filtered as absent in all samples: " + entry);
				}
				somaticSNVs.add(entry);	
			}
			rd.close();
			logger.log(Level.INFO, "There are " + totalSNVCounter + " SNVs in validation file. Of those, we pass "+ germlineCounter + " as germline, and "+ somaticSNVs.size() +" as somatic. \n");
		} catch (IOException e){
			System.err.println("ERROR: could not read file: " + inputFile);
			System.exit(1);
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
			System.err.println("CNV input file Reading Error!"+ e);
		}
		return CNVs;
	}
	
	private ArrayList<String> loadAnnovarFunction(String inputAnnFile){
		ArrayList<String> anns = new ArrayList<String>();
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputAnnFile));
			String currLine = rd.readLine();
			while (currLine != null){
				String[] entryParts = currLine.split("\t");
				String ann = "";
				ann += entryParts[0] + " ";
				ann += entryParts[1];
				anns.add(ann);
				currLine = rd.readLine();
			}
			rd.close();
		} catch (IOException e) {
			System.err.println("ANN input file Reading Error!");
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
			System.err.println("COSMIC input file Reading Error!");
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
			System.err.println("TCGA input file Reading Error!");
		}
		return tcgaDB;
	}
	
	public void annotateSNVs(String inputCNVFile, String inputAnnFile, String cosmicFile, String tcgaFile) {
		ArrayList<CNVRegion> cnvs = null;
		if(inputCNVFile != null) {
			cnvs = loadCNVs(inputCNVFile);
		}
		ArrayList<String> anns = null;
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
				entry.getAnnotation().annovar = anns.get(i);
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
