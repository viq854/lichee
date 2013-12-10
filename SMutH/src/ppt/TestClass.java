package ppt;

import java.util.ArrayList;

import util.*;



//import java.util.ArrayList;

public class TestClass {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//SNVDatabase db1 = new SNVDatabase("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_1/Patient_1.validation.txt", 0);
		//SNVDatabase db2 = new SNVDatabase("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_2/Patient_2.validation.txt", 0);
		//SNVDatabase db3 = new SNVDatabase("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_3/Patient_3.validation.txt", 1);
		//SNVDatabase db4 = new SNVDatabase("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_4/Patient_4.validation.txt", 1);
		SNVDatabase db5 = new SNVDatabase("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_5/Patient_5.validation.txt", 1);
		//SNVDatabase db6 = new SNVDatabase("/Users/rahelehs/Work/BreastCancer/patients_vcfs/full_vcfs/Patient_6/Patient_6.validation.txt", 0);
		/*int count = 0;
		
		for (SNVEntry entry :db1.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='1' && entry.getAAF(6) < 0.1 && entry.getAAF(6) > 0.01){
				count++;
				System.out.println(entry.getAAF(6));
			}
			
		}
		
		for (SNVEntry entry :db2.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='1'  && entry.getAAF(4) < 0.1  && entry.getAAF(4) > 0.01){
				count++;
				System.out.println(entry.getAAF(4));
			}
			
		}
		
		for (SNVEntry entry :db3.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='1' && entry.getAAF(9) < 0.1  && entry.getAAF(9) > 0.01){
				count++;
				System.out.println(entry.getAAF(9));
			}
			
		}
		for (SNVEntry entry :db4.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='1' && entry.getAAF(6) < 0.1 && entry.getAAF(6) >0.01){
				count++;
				System.out.println(entry.getAAF(6));
			}
			
		}
		
		
		for (SNVEntry entry :db5.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='1' && entry.getAAF(10) < 0.1&& entry.getAAF(10) > 0.01){
				count++;
				System.out.println(entry.getAAF(10));
			}
			
		}
		
		for (SNVEntry entry :db6.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='1' && entry.getAAF(11) < 0.1  && entry.getAAF(11) > 0.01){
				count++;
				System.out.println(entry.getAAF(11));
			}
			
		}
		
		for (SNVEntry entry :db1.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='0' && entry.getAAF(6) > 0.01){
				System.out.println(entry.getAAF(6));
			}
			
		}
		
		for (SNVEntry entry :db2.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='0' && entry.getAAF(4) > 0.01){
				System.out.println(entry.getAAF(4));
			}
			
		}
		
		for (SNVEntry entry :db3.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='0' && entry.getAAF(9) > 0.01){
				System.out.println(entry.getAAF(9));
			}
			
		}
		for (SNVEntry entry :db4.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='0' && entry.getAAF(6) > 0.01){
				System.out.println(entry.getAAF(6));
			}
			
		}
		
		
		for (SNVEntry entry :db5.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='0' && entry.getAAF(10) > 0.01){
				System.out.println(entry.getAAF(10));
			}
			
		}
		
		for (SNVEntry entry :db6.somaticSNVs){
			String tag = ((MUTEntry) entry).getEOATag();
			if (tag.charAt(tag.length()-1) =='0' && entry.getAAF(11) > 0.01){
				System.out.println(entry.getAAF(11));
			}
			
		}
		
		System.out.println("Hello, world!"+count);*/
		//VCFDatabase db = new VCFDatabase("testCases/simulation_vcfs/tree_4_03.raw.vcf");
		//db.generateMatrix("output");
		//ArrayList<VCFEntry> test = db.getSortedEntriesByGATK("0101", "0100");
//		for (VCFEntry entry: test){
//			System.out.println(entry.toString());
//		}
	}

}
