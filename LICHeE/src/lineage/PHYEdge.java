package lineage;

/**
 * Edge in the phylogenetic constraint network
 */
public class PHYEdge {
	protected PHYNode from;
	protected PHYNode to;
	
	public PHYEdge(PHYNode from, PHYNode to) {
		this.from = from;
		this.to = to;
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof PHYEdge)) {
			return false;
		}
		PHYEdge e = (PHYEdge) o;
		if(this.from.equals(e.from) && this.to.equals(e.to)) {
			return true;
		} else {
			return false;
		}
	}
}