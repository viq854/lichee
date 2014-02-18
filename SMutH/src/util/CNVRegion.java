package util;

public class CNVRegion {

	private int chr;
	private int startPos;
	private int endPos;
	
	//a sequence of ploidy numbers in different samples. e.g. tag 221 means there is loss in the third sample.
	private String ploidy;
	
	private boolean whole;

	
	
	public CNVRegion(String entry){
		String[] parts = entry.split("\t");

		if (parts[0].charAt(3) == 'X') chr=23;
		else if (parts[0].charAt(3) == 'Y') chr=24;
		else chr=Integer.parseInt(parts[0].substring(3));
		
		if (parts.length == 2){
			ploidy = parts[1];
			whole = true;
		}else{
			startPos = Integer.parseInt(parts[1]);
			endPos = Integer.parseInt(parts[2]);
			ploidy = parts[3];
			whole = false;
		}
		
	}
	
	public CNVRegion(int chrom, int start, int end){
		chr = chrom;
		startPos = start;
		endPos = end;
		whole = false;
	}
	
	public CNVRegion(int chrom){
		chr = chrom;
		whole = true;
	}

	
	public String toString(){
		String str = chr +"\t"+ startPos+"\t"+endPos;
		return str;
	}

	
	/**
	 * -1 before, 0 inside, 1 after
	 * @param entry
	 * @return
	 */
	public int compareLocation(int SNVChr, int SNVPos){
		if (SNVChr > chr )
			return 1;
				
		if (SNVChr < chr  )
			return -1;
	
		if (whole)
			return 0;
		
		if (SNVPos < startPos)
			return -1;
			
		if (SNVPos > endPos)
			return 1;
		
		return 0;


	}

}
