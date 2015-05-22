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

import util.CNVRegion;

public class SNVEntry {
	
	/** SNV entry info */
	protected int id;
	protected int chr;
	protected int position;
	protected String description;
	protected String presenceProfile;
	protected String ambigPresenceProfile;
	protected double[] VAF;
	protected boolean isRobust;
	protected String snvEntryString;
	protected boolean inCNVRegion;
	protected String annotation;
	 
	public SNVEntry(String line, int lineId) {
		id = lineId;
		snvEntryString = line;
		annotation = "";
	}
	
	public int getId() {
		return id;
	}
	
	/** Returns the SNV chromosome */
	public int getChromosome() {
		return chr;
	}
	
	/** Returns the SNV position  */
	public int getPosition() {
		return position;
	}
	
	/** Returns the description field */
	public String getDescription() {
		return description;
	}
	
	/** Returns the VAF in sample i */
	public double getVAF(int i) {
		return VAF[i];
	}
	
	/** Returns true if the SNV was robustly called in all samples  */
	public boolean isRobust() {
		return isRobust;
	}
	
	/** Returns the sample presence-absence profile */
	public String getProfile() {
		return presenceProfile;
	}
	
	public String getAmbigProfile() {
		return ambigPresenceProfile;
	}
	
	/** Sets the sample presence-absence profile */
	public void updateGroup(String code) {
		presenceProfile = code;
	}
	
	/** Returns true if profile[sampleId] == 1 */
	public boolean isPresent(int sampleId) {
		return presenceProfile.charAt(sampleId) == '1';
	}
	
	public boolean evidenceOfPresence(int sample){
		return (VAF[sample] > Parameters.MAX_VAF_ABSENT );
	}
	
	public String getAnnotation() {
		return annotation;
	}
	
	public void addAnnotation(String ann) {
		annotation += ann + "; ";
	}

	public boolean checkInCNVRegion(ArrayList<CNVRegion> CNVs) {
		for (CNVRegion cnv : CNVs){
			int loc = cnv.compareLocation(chr, position);
			if (loc == 0) {
				inCNVRegion = true;
				return true;
			}
			if (loc == -1) {
				inCNVRegion = false;
				return false;
			}
		}
		inCNVRegion = false;
		return false;
	}
	
	public boolean isInCNVRegion() {
		return inCNVRegion;
	}
	
	public String toString() {
		return snvEntryString;
	}
}
