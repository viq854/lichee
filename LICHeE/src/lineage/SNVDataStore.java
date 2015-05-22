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

import lineage.AAFClusterer.Cluster;
import lineage.Parameters.Format;
import util.CNVRegion;

/** 
 * Handles input SNV data loading, storing, filtering, and partitioning 
 * based on sample presence-absence profiles
 */
public class SNVDataStore {
	/** Number of input samples */
	private int numSamples;
	/** Input sample names */
	private ArrayList<String> sampleNames;
	/** List of all input somatic SNVs */
	protected ArrayList<SNVEntry> somaticSNVs;
	/** List of SNVs with ambiguous VAFs in at least one sample */
	protected ArrayList<SNVEntry> ambiguousSNVs;
	/** Map of sample profile tags to a list of SNVs with this profile */
	private HashMap<String, ArrayList<SNVEntry>> tag2SNVs;
	private HashMap<String, ArrayList<Cluster>> tag2Clusters;
	/** Index of the normal sample in the input sample list*/
	private int normalSample;
	private static Logger logger = LineageEngine.logger;

	public SNVDataStore(String snvInputFile, String clusterInputFile, int normalSampleId) {
		normalSample = normalSampleId;
		somaticSNVs = new ArrayList<SNVEntry>();
		tag2SNVs = new HashMap<String, ArrayList<SNVEntry>>();
		tag2Clusters = new HashMap<String, ArrayList<Cluster>>();
		ambiguousSNVs = new ArrayList<SNVEntry>();	
		// load and process input SNVs
		if(clusterInputFile == null) {
			loadSNVFile(snvInputFile);		
		} else {
			loadSNVFileWithClusters(snvInputFile, clusterInputFile);
		}
		
		logger.fine("Sample id -> name map:");
		for(int i = 0; i < getNumSamples(); i++) {
			logger.fine(i + ": " + getSampleName(i));
		}
		logger.fine("Initial robust groups:");
		reportSNVGroups();
		if(Parameters.INPUT_FORMAT == Parameters.Format.SNV_WITH_PROFILE || clusterInputFile != null) return;
		
		// handle mutations from small groups as ambiguous
		ArrayList<String> smallGroups = new ArrayList<String>();
		for(String tag : tag2SNVs.keySet()) {
			if(tag2SNVs.get(tag).size() < Parameters.MIN_GROUP_PROFILE_SUPPORT) {
				for(SNVEntry entry : tag2SNVs.get(tag)) {
					ambiguousSNVs.add(entry);
				}
				smallGroups.add(tag);
			}
		}
		for(String group : smallGroups) {
			tag2SNVs.remove(group);
		}
		// assign ambiguous SNVs to existing groups or create new groups
		assignAmbiguousSNVs();
		// process resulting groups
		filterGroups();
		logger.fine("Final SNV groups:");
		reportSNVGroups();
	}
	
	/**
	 * Assign ambiguous SNVs to existing groups; 
	 * create new groups for SNVs with no suitable matches
	 */
	private void assignAmbiguousSNVs() {
		if (ambiguousSNVs.size() == 0) return;
		
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

		logger.log(Level.FINE, "Ambiguous profile SNV assignment (" + ambiguousSNVs.size() + " total): ");
		ArrayList<SNVEntry> toRemove = new ArrayList<SNVEntry>();
		for(SNVEntry snv : ambiguousSNVs) {	
			double bestDistToTarget = 0;
			String bestTarget = "";
				
			for(String target : targetTags) {				
				if(!canConvert(snv, target)) continue;
				// if the target is germline, move to germline regardless of distance
				if(target.equals(all1s)) {
					toRemove.add(snv);
					logger.log(Level.FINE, "**Removed as germline:\n" + snv);
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
				logger.log(Level.FINE, "Assigned " + snv.getAmbigProfile() + " to " + bestTarget + " with dist " + bestDistToTarget + ": " + snv);
				continue;
			}
			
			if(bestDistToTarget != 0 && bestDistToTarget >= getHammingWeight(bestTarget)*Parameters.MIN_VAF_TARGET_RATIO_PER_SAMPLE) {
				// found a valid match
				snv.updateGroup(bestTarget);
				toRemove.add(snv);
				tag2SNVs.get(bestTarget).add(snv);
				logger.log(Level.FINE, "Assigned " + snv.getAmbigProfile() + " to " + bestTarget + " with dist " + bestDistToTarget + ": " + snv);
				continue;
			}
			
			logger.log(Level.FINE, "**No existing candidate group found for " + snv.getAmbigProfile() + "" + ((bestDistToTarget != 0) ? " (best dist = " + bestDistToTarget + " to target " + bestTarget + "): " : ": ") + snv);

		}
		ambiguousSNVs.removeAll(toRemove);
		
		// remaining snvs had no suitable matches and potentially represent true branches
		// we minimize the number of additional nodes by merging the groups
		HashMap<String, ArrayList<SNVEntry>> ambiguousGroups = mergeAmbiguousSNVs(ambiguousSNVs);
		
		for(String tag : ambiguousGroups.keySet()) {
			if(tag.equals(all0s)) continue;
			if(!tag2SNVs.containsKey(tag)) {
				tag2SNVs.put(tag, ambiguousGroups.get(tag));
				logger.log(Level.FINE, "Created new group: " + tag);
			}
		}
		tag2SNVs.remove(all0s);
	}
	
	/**
	 * Finds the minimum number of groups that can incorporate the ambiguous input SNVs
	 * applying the greedy set cover algorithm
	 */
	private HashMap<String, ArrayList<SNVEntry>> mergeAmbiguousSNVs(ArrayList<SNVEntry> snvs) {
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
				logger.log(Level.FINE, "Assigned " + entry.getAmbigProfile() + " to " + maxSet + ": " + entry);
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
				if(snv.getProfile().charAt(i) == '0' && snv.evidenceOfPresence(i)) {
					double delta0 = snv.getVAF(i) - Parameters.MAX_VAF_ABSENT;
					double delta1 = Parameters.MIN_VAF_PRESENT - snv.getVAF(i);
					tag += (delta1 < delta0) ? '1' : snv.getProfile().charAt(i);
				} else {
					tag += snv.getProfile().charAt(i);
				}
			}
			logger.log(Level.FINE, "Assigned " + snv.getAmbigProfile() + " to " + tag + ": " + snv);
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
		extendTarget(partialProfile + ""+ snv.getProfile().charAt(sample), snv, targets);
		if(snv.getProfile().charAt(sample) == '0' && snv.evidenceOfPresence(sample)) {
			extendTarget(partialProfile + '1', snv, targets);
		}
	}
	
	private double distToZero(SNVEntry snv) {
		String tag = snv.getProfile();
		double dist = 0;
		for(int i = 0; i < numSamples; i++) {
			if(tag.charAt(i) == '1' || snv.evidenceOfPresence(i)) {
				dist += 0.0001/snv.getVAF(i);
			}
		}
		return dist;
	}
	
	private double distToTargetVAF(SNVEntry snv, SNVEntry targetSNV) {
		String tag = snv.getProfile();
		String targetTag = targetSNV.getProfile();
		double dist = 0;
		for(int i = 0; i < numSamples; i++) {
			if(targetTag.charAt(i) == '0' && tag.charAt(i) == '0' && !snv.evidenceOfPresence(i)) continue;
			double min = targetSNV.getVAF(i) < snv.getVAF(i) ? targetSNV.getVAF(i) : snv.getVAF(i);
			double max = targetSNV.getVAF(i) > snv.getVAF(i) ? targetSNV.getVAF(i) : snv.getVAF(i);
			if(min == 0) {
				min = 0.0001;
			}
			dist += min/max;
		}
		return dist;
	}
	
	private boolean canConvert(SNVEntry snv, String target) {
		String tag = snv.getProfile();		
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
	
	/**
	 * Group filtering based on:
	 * - group size
	 * - absence in all samples
	 * - minimum number of robust SNVs
	 */
	private void filterGroups() {
		// apply minimum size and robust size constraint
		ArrayList<String> filteredOut = new ArrayList<String>();
		for(String tag : tag2SNVs.keySet()) {
			if(isAbsent(tag) || (tag2SNVs.get(tag).size() < Parameters.MIN_SNVS_PER_GROUP)) {
				filteredOut.add(tag);
				continue;
			}
	
			if(Parameters.MIN_ROBUST_SNVS_PER_GROUP > 0) {
				int numRobust = 0;
				for(SNVEntry entry : tag2SNVs.get(tag)) {
					if(entry.isRobust()) {
						numRobust++;
						if(numRobust >= Parameters.MIN_ROBUST_SNVS_PER_GROUP) {
							break;
						}
					}
				}
				if(numRobust < Parameters.MIN_ROBUST_SNVS_PER_GROUP) {
					filteredOut.add(tag);
				}
			}
		}
		for(String group : filteredOut) {
			tag2SNVs.remove(group);
		}
	}
	
	/**
	 * SNV filtering based on:
	 * - germline SNV
	 * - VAF higher than max threshold
	 * - absence in all samples
	 * 
	 * SNVs that pass the above filters are stored based
	 * on their robustness
	 */
	private void processSNVEntry(SNVEntry entry) {
		if(Parameters.INPUT_FORMAT != Format.SNV_WITH_PROFILE && entry.isPresent(normalSample)) {
			logger.log(Level.INFO, "**Filtered as germline: \n" + entry);
			return;
		}
		if(!hasValidVAFs(entry)) {
			logger.log(Level.INFO, "**Filtered due to VAFs > allowed MAX: \n" + entry);
			return;
		}
		if(entry.isRobust() && isAbsent(entry.getProfile())) {
			logger.log(Level.INFO, "**Filtered as robustly absent in all samples: \n" + entry);
			return;
		}
		somaticSNVs.add(entry);	
		
		if(entry.isRobust()) {
			String tag = entry.getProfile();
			if (!tag2SNVs.containsKey(tag)){
				tag2SNVs.put(tag, new ArrayList<SNVEntry>());	
			}
			tag2SNVs.get(tag).add(entry); 
		} else {
			ambiguousSNVs.add(entry);
		}
	}
	
	/**
	 * Checks for VAFs that are above a certain threshold
	 */
	private boolean hasValidVAFs(SNVEntry snv) {
		for(int i = 0; i < numSamples; i++) {
			// check if the VAF i too high
			if(snv.getVAF(i) > Parameters.MAX_ALLOWED_VAF) {
				return false;
			}
		}
		return true;
	}
	
	/** 
	 * Checks if the SNV is absent in all the samples
	 */
	private boolean isAbsent(String tag) {
		for(int i = 0; i < numSamples; i++) {
			if(tag.charAt(i) == '1') {
				return false;
			}
		}
		return true;
	}
	
	/****** Accessors *****/
	
	public int getNumSamples() {
		return numSamples;
	}
	
	public String getSampleName(int i) {
		return sampleNames.get(i);
	}
	
	public ArrayList<String> getSampleNames() {
		return sampleNames;
	}
	
	public HashMap<String, ArrayList<SNVEntry>> getSomaticGroups() {
		return tag2SNVs;
	}
	
	public HashMap<String, ArrayList<Cluster>> getClusters() {
		return tag2Clusters;
	}
	
	
	/** A group is robust if it contains sufficient robust mutations */
	public boolean isRobustGroup(String groupTag) {
		int numRobust = 0;
		for(SNVEntry entry : tag2SNVs.get(groupTag)) {
			if(entry.isRobust()) {
				numRobust++;
				if(numRobust >= Parameters.MIN_GROUP_PROFILE_SUPPORT) {
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

		logger.log(Level.FINE, "Profile\t#SNVs\t#Robust\tMean VAF(Robust Mean VAF)");	
		for (String tag : tags){
			if (tag2SNVs.get(tag).size() == 0) {
				continue;
			}
			
			double[] mean = new double[numSamples];	
			double[] robustMean = new double[numSamples];	
			int robustCount = 0;
			for (SNVEntry entry : tag2SNVs.get(tag)) {
				for (int i = 0; i < numSamples; i++){
					mean[i] += entry.getVAF(i);
					if(entry.isRobust()) {
						robustMean[i] += entry.getVAF(i);
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
	
	////// SNV File I/O //////
	
	private void loadSNVFile(String inputFile) {
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String currLine = rd.readLine();
			if(currLine == null) {
				returnInputFileFormatError("Empty file", null); 
			}
			parseHeader(currLine);
			
			currLine = rd.readLine();
			int totalSNVCounter = 0;
			while (currLine != null) {
				SNVEntry entry = parseSNVEntry(currLine, totalSNVCounter+1);
				totalSNVCounter++;
				processSNVEntry(entry);
				currLine = rd.readLine();
			}
			rd.close();
			logger.log(Level.INFO, "There are " + totalSNVCounter + " SNVs in the input file. \nAfter pre-processing, the input consists of " + somaticSNVs.size() +" somatic SNVs. \n");
		} catch (IOException e){
			returnInputFileFormatError("Could not read file: " + inputFile, null);
		}
		
	}
	
	private void loadSNVFileWithClusters(String inputFile, String clustersFile) {
		// load the input SNV file 
		try {
			BufferedReader rd = new BufferedReader(new FileReader(inputFile));
			String currLine = rd.readLine();
			if(currLine == null) {
				returnInputFileFormatError("Empty file", null); 
			}
			parseHeader(currLine);
			
			currLine = rd.readLine();
			int totalSNVCounter = 0;
			while (currLine != null) {
				SNVEntry entry = parseSNVEntry(currLine, totalSNVCounter+1);
				totalSNVCounter++;
				somaticSNVs.add(entry);	// no filtering
				currLine = rd.readLine();
			}
			rd.close();
			logger.log(Level.INFO, "There are " + totalSNVCounter + " SNVs in the input file.\n");
		} catch (IOException e){
			returnInputFileFormatError("Could not read file: " + inputFile, null);
		}
		// load the clusters file
		try {
			BufferedReader rd = new BufferedReader(new FileReader(clustersFile));
			String currLine = rd.readLine();
			if(currLine == null) {
				returnInputFileFormatError("Empty file", null); 
			}
			int numClusters = 0;
			AAFClusterer clusterer = new AAFClusterer();
			while (currLine != null) {
				double[] vafs = new double[numSamples];
				ArrayList<Integer> members = new ArrayList<Integer>();
				String profile = parseSNVCluster(currLine, vafs, members);
				int numSetSamples = getHammingWeight(profile);
				double[] centroid = new double[numSetSamples];
				int idx = 0;
				for(int i = 0; i < numSamples; i++) {
					if(profile.charAt(i) == '1') {
						centroid[idx] = vafs[i];
						idx++;
					}
				}
				
				// new cluster
				Cluster c = clusterer.new Cluster(centroid, numClusters);
				numClusters++;
				c.setStdDev(new double[numSamples]);
				int startId = 0;
				if(tag2SNVs.containsKey(profile)) {
					startId = tag2SNVs.get(profile).size();
				} else {
					tag2SNVs.put(profile, new ArrayList<SNVEntry>());
				}
				for(int i = 0; i < members.size(); i++) {
					c.addMember(startId + i);
					SNVEntry snv = somaticSNVs.get(members.get(i)-1);
					snv.presenceProfile = profile;
					snv.isRobust = true;
					tag2SNVs.get(profile).add(snv);
				}
				if(!tag2Clusters.containsKey(profile)) {
					tag2Clusters.put(profile, new ArrayList<Cluster>());
				}
				tag2Clusters.get(profile).add(c);
				currLine = rd.readLine();
			}
			rd.close();
			logger.log(Level.INFO, "There are " + numClusters + " clusters in the input file.\n");
		} catch (IOException e){
			returnInputFileFormatError("Could not read file: " + clustersFile, null);
		}
	}
	
	// Format information
	protected static final int NUM_REQ_FIELDS_SNV_FILE = 3;
	protected static final int NUM_REQ_FIELDS_SNV_W_PROFILE_FILE = 4;
	public static int getNumRequiredFields() {
		return Parameters.INPUT_FORMAT == Format.SNV ? NUM_REQ_FIELDS_SNV_FILE : NUM_REQ_FIELDS_SNV_W_PROFILE_FILE;
	}
		
	public static void returnInputFileFormatError(String desc, String entry) {
		System.err.println("[Wrong input file format] " + desc);
		if(entry != null) {
			System.err.println("Line: " + entry);
		}
		System.exit(1);
	}
	
	private void parseHeader(String headerLine) {	
		int numRequiredFields = getNumRequiredFields();
		String[] headerParts = headerLine.split("\t");
		if(headerParts.length < numRequiredFields + 2) {
			returnInputFileFormatError("The header must have " + numRequiredFields + " required fields and at least 2 samples (separated by tabs)", headerLine); 
		}
		sampleNames = new ArrayList<String>(Arrays.asList(headerParts).subList(numRequiredFields, headerParts.length));
		numSamples = sampleNames.size();
		logger.log(Level.FINE, "Input file contains " + numSamples + " samples!");
	}
	
	/** Returns the SNV chromosome as an integer value */
	private int parseChrNum(String chrString) {
		if(chrString.contains("chr")) {
			chrString = chrString.substring(3, chrString.length());
		}
		if(chrString.equals("x")) return 23;
		if(chrString.equals("y")) return 24;
		try {
			return Integer.parseInt(chrString);
		} catch (Exception e) {
			return -1;
		}
	}
	
	private SNVEntry parseSNVEntry(String line, int lineId) {
		SNVEntry entry = new SNVEntry(line, lineId);
		
		// parse the entry fields
		String[] entryParts = line.split("\t");
		int numFields = getNumRequiredFields();
		if(entryParts.length != (numFields + numSamples)) {
			returnInputFileFormatError("Expecting " + numFields + " fields and " + numSamples + " samples based on the header", line);
		}
		
		entry.chr = parseChrNum(entryParts[0].trim().toLowerCase());
		if(entry.chr == -1) {
			returnInputFileFormatError("Chromosome " + entryParts[0], line);
		}
		try {
			entry.position = Integer.parseInt(entryParts[1].trim());
		} catch (NumberFormatException e) {
			returnInputFileFormatError("Position " + entryParts[1].trim(), line);
		}	 
		entry.description = entryParts[2].trim();
		if(Parameters.INPUT_FORMAT == Format.SNV_WITH_PROFILE) {
			entry.presenceProfile = entryParts[3].trim();
			if(entry.presenceProfile.length() != numSamples) {
				returnInputFileFormatError("Presence profile " + entry.presenceProfile + " length does not match the number of input samples", line);
			}
		} else {
			entry.presenceProfile = "";
			entry.ambigPresenceProfile = "";
		}
		
		// parse per sample VAF values
		entry.VAF = new double[numSamples];
		entry.isRobust = true;
		for(int i = 0; i < numSamples; i++) {
			try {
				entry.VAF[i] = Double.parseDouble(entryParts[numFields + i]);
			} catch (NumberFormatException e) {
				returnInputFileFormatError("VAF value " + entryParts[i + numFields], line);
			}
			if(Parameters.INPUT_FORMAT == Format.SNV_WITH_PROFILE) continue;
			if (entry.VAF[i] <  Parameters.MIN_VAF_PRESENT){
				if (entry.VAF[i] >=  Parameters.MAX_VAF_ABSENT) {
					entry.isRobust = false;
					entry.ambigPresenceProfile += "*";
				} else {
					entry.ambigPresenceProfile += "0";
				}
				entry.presenceProfile += "0";
			} else { 
				entry.presenceProfile += "1";
				entry.ambigPresenceProfile += "1";
			}
		}
		
		return entry;
	}
	
	private String parseSNVCluster(String line, double[] centroid, ArrayList<Integer> members) {	
		// parse the entry fields
		String[] entryParts = line.split("\t");
		int numFields = 2;
		if(entryParts.length != (numFields + numSamples)) {
			returnInputFileFormatError("Expecting " + numFields + " fields and " + numSamples + " samples based on the SSNV file", line);
		}
		
		String profile = entryParts[0].trim();
		if(profile.length() != numSamples) {
			returnInputFileFormatError("Cluster presence profile " + profile + " length does not match the number of input samples", line);
		}
		for(int i = 0; i < numSamples; i++) {
			try {
				centroid[i] = Double.parseDouble(entryParts[i + 1]);
			} catch (NumberFormatException e) {
				returnInputFileFormatError("Cluster centroid VAF value " + entryParts[i + 1], line);
			}
		}
		
		String[] snvEntries = entryParts[numSamples + 1].split(",");
		if(snvEntries.length == 0) {
			returnInputFileFormatError("Cluster SSNV list is empty or malformated (SNVs should be separated by commas)", line);
		}
		
		for(String s : snvEntries) {
			try {
				int snvId = Integer.parseInt(s.trim());
				members.add(snvId);
			} catch (NumberFormatException e) {
				returnInputFileFormatError("SSNV member entry " + s.trim() + " is not a valid number ", line);
			}
		}
		return profile;
	}
	
	////// Additional SNV annotations ////
	
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
				entry.addAnnotation(anns.get(i));
			}
			if(cosmicDB != null) {
				if(cosmicDB.containsKey(entry.getPosition())) {
					entry.addAnnotation(cosmicDB.get(entry.getPosition()));
				}
			}
			if(tcgaDB != null) {
				if(tcgaDB.containsKey(entry.getChromosome()) && tcgaDB.get(entry.getChromosome()).contains(entry.getPosition())) {
					entry.addAnnotation("+TCGA+");
				}
			}
		}
	}
}
