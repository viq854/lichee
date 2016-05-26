package util;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Toolkit;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.sql.Timestamp;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.swing.Icon;
import javax.swing.JPanel;

import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import lineage.PHYNetwork;
import lineage.PHYNode;
import lineage.PHYTree;
import lineage.SNVEntry;


public class LineageDisplayConfig {

	// colors
	static class ColorPalette {
		public static Color[] palette = {new Color(255,255,255), new Color(200,155,82), new Color(100,100,117), new Color(194,110,144), 
				new Color(255,101,101), new Color(212,196,184), new Color(3,115,129), new Color(23,118,191), new Color(109,113,63), 
				new Color(146,185,132), new Color(207,46,28), 
				new Color(47,75,124), new Color(254,224,0), new Color(205,92,92), new Color(0,128,128), new Color(106,90,205),
				new Color(136,232,190), new Color(33,57,61), new Color(60,144,163), new Color(233,180,0), new Color(191,93,71),
				new Color(40,95,107), new Color(112,143,0), new Color(130,157,186), new Color(211,51,63), new Color(147,53,69)};
		
		public static Color genColorDiff(int i, int n) {
			return Color.getHSBColor((float)i/(float)n, 0.75f, 1.0f);
		}
		
		public static Color pastelify(Color c) {
			 return new Color((c.getRed() + 255)/2, (c.getGreen() + 255)/2, (c.getBlue() + 255)/2);
		}
		
		public static Color genColorRand(Random rand) {
		    return new Color((rand.nextInt(256) + 255)/2, (rand.nextInt(256) + 255)/2, (rand.nextInt(256) + 255)/2);
		}
		
		public static Color getEdgeColor(Color c) {
			return new Color(c.getRed(), c.getGreen(), c.getBlue(), 0x77);
		}
		
		public static Color getBorderColor(Color c) {
			return c.darker();
		}
		
		public static Color getBorderHighlighColor(Color c) {
			return HIGHLIGHT_COLOR;
		}
		
		public static Color getNodeHighlighColor(Color c) {
			return c;
		}
		
		public static Color makeTransparent(Color c) {
			return new Color(c.getRed(), c.getGreen(), c.getBlue(), 0x77);
		}
		
		public static Color SELECTED_NODE_COLOR = Color.yellow;
		public static Color INTERNAL_NODE_BORDER_COLOR = Color.BLACK;
		public static Color LEAF_NODE_BORDER_COLOR = Color.BLACK;
		public static Color INTERNAL_NODE_COLOR = Color.white;
		public static Color LEAF_NODE_COLOR = Color.white;
		public static Color LABEL_COLOR = Color.BLACK;
		public static Color EDGE_COLOR = Color.BLACK;
		public static Color HIGHLIGHT_COLOR = new Color(186,85,211);
	}
	
	// fonts
	public static Font GUI_FONT = new Font("Comic Sans MS", Font.BOLD, 16);
	public static Font CONSOLE_FONT = new Font("Comic Sans MS", Font.PLAIN, 14);
	public static Font INTERNAL_NODE_FONT = new Font("Comic Sans MS", Font.BOLD, 20);
	public static Font LEAF_NODE_FONT = new Font("Comic Sans MS", Font.BOLD, 20);
	public static Font EDGE_FONT = new Font("Comic Sans MS", Font.BOLD, 18);
	public static Font PLAIN_FONT = new Font("Arial", Font.PLAIN, 20);
	public static Font DOT_FONT = new Font("Arial", Font.BOLD, 24);
		
	// shapes
	public static Shape LEAF_NODE_SHAPE = new Rectangle2D.Float(-40, -40, 80, 80);
	public static Shape getInternalNodeShape(int size) {
        if(size <= 10) size = 15;
        int w = 4*size;
        if(w > 80) w = 100;
		return new Ellipse2D.Float(-w/2, -w/2, w, w);
	}
	
	protected boolean beautify;
	protected boolean subclone_path_highligth;
	protected boolean subclone_sample_highligth;
	protected Color[] nodeColors;
	protected boolean[] nodeHighlight;
	protected boolean[] nodeSelected;
	Icon[] nodeIcons;
	boolean[] nodeIconPainted;
	
	
	protected DirectedGraph<Integer, Integer> displayTree;
	protected PHYTree spanningTree;
	protected PHYNetwork constrNetwork;
	protected HashMap<Integer, PHYNode> nodeInfo;
	protected ArrayList<String> sampleNames;
	protected Random randGen;
	
	
	public LineageDisplayConfig(PHYNetwork network, PHYTree tree, ArrayList<String> samples, boolean colorful) {
		constrNetwork = network;
		spanningTree = tree;
		sampleNames = samples;
		beautify = colorful;
		computeDisplayTree();
		
		int numNodes = displayTree.getVertexCount() + 1;
		subclone_path_highligth = false;
		subclone_sample_highligth = false;
		nodeHighlight = new boolean[numNodes];
		nodeSelected = new boolean[numNodes];
		nodeIcons = new Icon[numNodes];
		nodeIconPainted = new boolean[numNodes];
		if(beautify) {
			int numColors = numNodes;
			randGen = new Random(85354597);
			nodeColors = new Color[numColors];
			for(int i = 0; i < numColors; i++) {
				nodeColors[i] = (i < ColorPalette.palette.length) ? ColorPalette.pastelify(ColorPalette.palette[i]) : ColorPalette.genColorRand(randGen); 
			}
		}
	}
	
	public void computeDisplayTree() {			
		displayTree = new DirectedSparseGraph<Integer, Integer>();	
		nodeInfo = new HashMap<Integer, PHYNode>();
		int edgeId = 0;
		for (PHYNode n : spanningTree.treeEdges.keySet()) {
			displayTree.addVertex(n.getNodeId());
			nodeInfo.put(n.getNodeId(), n);
			for(PHYNode n2 : spanningTree.treeEdges.get(n)) {
				if(!displayTree.containsVertex(n2.getNodeId())) {
					displayTree.addVertex(n2.getNodeId());
					nodeInfo.put(n2.getNodeId(), n2);
				}
				displayTree.addEdge(edgeId, n.getNodeId(), n2.getNodeId(), EdgeType.DIRECTED);
				edgeId++;
			}
		}
		
		// add sample leaves
		for(int i = 1; i < constrNetwork.numSamples; i++) {
			PHYNode n = new PHYNode(0, i, constrNetwork.numNodes + i);
			displayTree.addVertex(n.getNodeId());
			nodeInfo.put(n.getNodeId(), n);
			// find a parent in the closest higher level		 
			boolean found = false;
			ArrayList<PHYNode> parents = new ArrayList<PHYNode>();
			ArrayList<PHYNode> sameLevelParents = new ArrayList<PHYNode>();
			for(int j = n.getLevel() + 1; j <= constrNetwork.numSamples; j++) {
				ArrayList<PHYNode> fromLevelNodes = constrNetwork.nodes.get(j);
				if(fromLevelNodes == null) continue;
				for(PHYNode n2 : fromLevelNodes) {
					if(n2.getAAF(i) > 0) {
						boolean addEdge = true;
						for(PHYNode p : parents) {
							if(spanningTree.isDescendent(n2, p)) {
								addEdge = false;
								break;
							}
						}
						if(addEdge) {
							sameLevelParents.add(n2);
							parents.add(n2);
							found = true;
						}
					}
				}
				// remove nodes that are in same level that are connected
				ArrayList<PHYNode> toRemove = new ArrayList<PHYNode>();
				for(PHYNode n1 : sameLevelParents) {
					for(PHYNode n2 : sameLevelParents) {
						if(spanningTree.isDescendent(n1, n2)) {
							toRemove.add(n1);
						}
					}
				}
				sameLevelParents.removeAll(toRemove);
				
				for(PHYNode n2 : sameLevelParents) {
					displayTree.addEdge(edgeId, n2.getNodeId(), n.getNodeId());
					edgeId++;
				}
				sameLevelParents.clear();
			}
			if(!found) {
				displayTree.addEdge(edgeId, 0, n.getNodeId());
				edgeId++;
			}
		}				
	}
	
	private String edge2DOT(Integer e) {
		String s = "";
		Integer n1 = displayTree.getSource(e);
		Integer n2 = displayTree.getDest(e);
		s += n1 + " -> " + n2;
		s += "[";
		s += " label=\"" + getEdgeLabel(n1, n2) + "\"";
		s += " fontname=\"" + DOT_FONT.getFontName() + "\"";
		s += " fontsize=" + getEdgeFont().getSize() + "";
		s += "];\n";
		return s;
	}
	
	private String node2DOT(Integer v, String tempImageDir, long timestamp) {
		String s = v + " [";
		if(nodeInfo.get(v).isRoot()) {
			s += " shape=plaintext";
			s += " label=\"" + getInternalNodeLabel(v) + "\"";
		} else if(nodeInfo.get(v).isLeaf()) {
			if(beautify) {
				Icon icon = getNodeIcon(v);
				BufferedImage image = new BufferedImage(icon.getIconWidth(), icon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
				icon.paintIcon(new JPanel(), image.getGraphics(), 0, 0);
				File outputfile = new File(tempImageDir + "/img_" + timestamp + "_" + v + ".png");
				try {
					ImageIO.write(image, "png", outputfile);
				} catch (IOException e) {
					e.printStackTrace();
				}
				s += "image=\"" + outputfile.getAbsolutePath() + "\"";
			} 
			s += " shape=square";
			s += " label=\"" + getExternalNodeLabel(v) + "\"";
			s += " labelloc=b";
		} else {
			s += " shape=circle";
			s += " label=\"" + getInternalNodeLabel(v) + "\"";
		}
		s += " fontname=\"" + DOT_FONT.getFontName() + "\"";
		s += " fontsize=" + getNodeFont(v).getSize() + "";
		Color c = getNodeFillColor(v);
		s += " style=filled fillcolor=\"" + String.format("#%02x%02x%02x", c.getRed(), c.getGreen(), c.getBlue()) + "\"";
		c = Color.black; //getNodeBorderColor(v);
		s += " color=\"" + String.format("#%02x%02x%02x", c.getRed(), c.getGreen(), c.getBlue()) + "\"";
		s += " width=" + getNodeShape(v).getBounds().getWidth()/Toolkit.getDefaultToolkit().getScreenResolution() + "";
		s += " heigth=" + getNodeShape(v).getBounds().getHeight()/Toolkit.getDefaultToolkit().getScreenResolution() + "";
		s += "];\n";
	    
		return s;
	}
	
	public String toDOT(String tempImageDirParent) {
		String t = "";
		t += "digraph G { \n";
		t += "size =\"10,10\"\n";
		t += "forcelabels=true\n";
		for(Integer e : displayTree.getEdges()) {
			t += edge2DOT(e);
		}
		
		String tempDirName = "/lichee_dot_img_temp";
		long ts = 0;
		if(beautify) {
			(new File(tempImageDirParent + tempDirName)).mkdirs();
			ts = new Timestamp(Calendar.getInstance().getTime().getTime()).getTime();
		}
		String rank = "{ rank = sink; ";
		for(Integer v : displayTree.getVertices()) {
			t += node2DOT(v, tempImageDirParent+tempDirName, ts);
			if(nodeInfo.get(v).isLeaf()) {
				rank += v + "; ";
			}
		}
		rank += "} \n";
		t += rank;
		t += "}";
		return t;
	}
	
	public Color getNodeLabelColor(Integer nid) {
		return LineageDisplayConfig.ColorPalette.LABEL_COLOR;
	}
	
	public Color getNodeFillColor(Integer nid) {
		PHYNode n = nodeInfo.get(nid);
		if(n.isRoot()) return Color.WHITE;
		if(n.isLeaf()) {
			return LineageDisplayConfig.ColorPalette.LEAF_NODE_COLOR;
		}
		Color c = LineageDisplayConfig.ColorPalette.INTERNAL_NODE_COLOR;
		if(beautify) {
			c = nodeColors[n.getNodeId()];
		}
		if(subclone_path_highligth && nodeHighlight[nid]) {
			return ColorPalette.getNodeHighlighColor(c);
		} 
		return c;
	}
	
	public Color getNodeBorderColor(Integer nid) {
		PHYNode n = nodeInfo.get(nid);
		if(n.isRoot()) return Color.WHITE;
		if(n.isLeaf()) {
			return LineageDisplayConfig.ColorPalette.LEAF_NODE_BORDER_COLOR;
		}
		// internal
		Color c = LineageDisplayConfig.ColorPalette.INTERNAL_NODE_BORDER_COLOR;
		if(beautify) {
			c = ColorPalette.getBorderColor(nodeColors[n.getNodeId()]);
		}
		if(subclone_path_highligth && nodeHighlight[nid]) {
			return ColorPalette.getBorderHighlighColor(c);
		} 
		return c;
	}
	
	public Color getEdgeColor(Integer nid1, Integer nid2) {
		PHYNode n1 = nodeInfo.get(nid1);
		Color c = LineageDisplayConfig.ColorPalette.EDGE_COLOR;
		if(beautify) {
			c = ColorPalette.getEdgeColor(nodeColors[n1.getNodeId()]);
			if(n1.isRoot()) {
				c = ColorPalette.getEdgeColor(Color.BLACK);
			}
		} 
		if(subclone_path_highligth && nodeHighlight[nid1] && nodeHighlight[nid2]) {
			return ColorPalette.getBorderHighlighColor(c);
		} 
		return c;
		
	}
	
	public Font getNodeFont(Integer nid) {
		if(!beautify) return PLAIN_FONT;
		PHYNode n = nodeInfo.get(nid);
		if(n.isRoot() || n.isLeaf()) {
			return LEAF_NODE_FONT;
		} else {
			return INTERNAL_NODE_FONT;
		}
	}
	
	public Font getEdgeFont() {
		if(beautify) {
			return LineageDisplayConfig.EDGE_FONT;
		} else  {
			return LineageDisplayConfig.PLAIN_FONT;
		}
	}
	
	public Shape getNodeShape(Integer nid) {
		PHYNode n = nodeInfo.get(nid);
		if(n.isLeaf()) {
			return LineageDisplayConfig.LEAF_NODE_SHAPE;
		}
		return getInternalNodeShape(n.getSize()); 
	}
	
	// only leaf nodes are currently labeled 
	public String getExternalNodeLabel(Integer nid) {
		PHYNode n = nodeInfo.get(nid);
		if(n.isLeaf()) {
			return sampleNames.get(n.getLeafSampleId());
		} else {
			return "";
		}
	}
	
	public String getInternalNodeLabel(Integer nid) {
		PHYNode n = nodeInfo.get(nid);
		if(n.isRoot()) {
			return "GL";
		} else if (!n.isLeaf()) {
			return n.getSize() + "";
		}
		return "";
	}
	
	// only leaf edges are currently labeled 
	public String getEdgeLabel(Integer nid1, Integer nid2) {
		PHYNode n1 = nodeInfo.get(nid1);
		PHYNode n2 = nodeInfo.get(nid2);
		String label = "";
		if(beautify) {
			DecimalFormat df = new DecimalFormat("#.##");
			if(n2.isLeaf()) {
				label += df.format(n1.getAAF(n2.getLeafSampleId()));
			}
		}
		return label;
	}
	
	public String getConsoleText(Integer nid) {
		PHYNode n = nodeInfo.get(nid);
		if(n.isRoot()) return "";
		if(n.isLeaf()) { // sample 
			return spanningTree.getLineage(n.getLeafSampleId(), sampleNames.get(n.getLeafSampleId()));		                	
        } else {
        	String s = "***Cluster " + n.getNodeId() + "***\n";
    		s += "Binary sample presence profile: " + n.getSNVGroup().getTag() + "\n";
    		int[] sampleIds = n.getSNVGroup().getSampleIds();
    		s+= "Samples (in order):";
    		for(int i = 0; i < n.getSNVGroup().getNumSamples(); i++) {
    			s+= " " + sampleNames.get(sampleIds[i]);
    		}
    		s += "\n";
    		s += n.getCluster().toString();
        	s += "\nMutations (chr:pos, info):\n";
         	ArrayList<SNVEntry> snvs = n.getSNVs(n.getSNVGroup().getSNVs());         	
         	for(SNVEntry snv : snvs) {
         		s += snv.getChromosome() + ":";
         		s += snv.getPosition() + ".........";
         		s += snv.getDescription() + "\n";
         	}
         	return s;
         }
	}
	
	public Icon getNodeIcon(final Integer nid) {
		final PHYNode n = nodeInfo.get(nid);
		final Shape s = getNodeShape(nid);
		final int w = (int) Math.floor(s.getBounds2D().getWidth());
		final int h = (int) Math.floor(s.getBounds2D().getHeight());
		if(!n.isLeaf()) {
			nodeIcons[nid] = new Icon() {
				public int getIconHeight() {
					return h;
				}
				public int getIconWidth() {
					return w;
				}
				public void paintIcon(Component c, Graphics g, int x, int y) {
					g.setFont(getNodeFont(nid));
					g.setColor(getNodeLabelColor(nid));
					String label = getInternalNodeLabel(nid);
					FontMetrics fm = g.getFontMetrics();
				    int x1 = (int) (w - fm.stringWidth(label))/2 + x;
				    int y1 = (int) (fm.getAscent() + (h - (fm.getAscent() + fm.getDescent())) / 2) + y;
				    g.drawString(label, x1, y1);
				}
			};
			return nodeIcons[nid];
		}
		
		final int sampleId = n.getLeafSampleId();
		ArrayList<ArrayList<PHYNode>> lineages = new ArrayList<ArrayList<PHYNode>>();
		spanningTree.getLineageClusters(new ArrayList<PHYNode>(), lineages, spanningTree.getRoot(), n.getLeafSampleId());
		final ArrayList<ArrayList<PHYNode>> paths = lineages;
		
		nodeIcons[nid] = new Icon() {
			public int getIconHeight() {
				return h;
			}
			public int getIconWidth() {
				return w;
			}
			public void paintIcon(Component c, Graphics g, int x, int y) {
				if(!beautify) {
					return;
				}
				
				Rectangle[] boxes = new Rectangle[constrNetwork.numNodes];
				int[] boxOffsetX = new int[constrNetwork.numNodes];
				int[] boxOffsetY = new int[constrNetwork.numNodes];
				int offsetY = 0;
				int X = x+1;
				int Y = y+1;
				int hBOX = (int) h-2;
				int wBOX = (int) w-2;
				int aBOX = hBOX*wBOX;
				ArrayList<Rectangle> highlight = new ArrayList<Rectangle>();
				for(int i = 0; i < paths.size(); i++) {
					ArrayList<PHYNode> path = paths.get(i);
					for(int j = 0; j < path.size(); j++) {
						PHYNode cluster = path.get(j);
						
						Rectangle r = boxes[cluster.getNodeId()];
						if(r == null) { // must create a new box for the node
							if(j == 0) { // no parent (always vertical)
								int H = (int) Math.floor(2*cluster.getAAF(sampleId)*hBOX);
								r = new Rectangle(X, Y + offsetY, wBOX, H);
								offsetY += H;
							} else {
								int parentId = path.get(j-1).getNodeId();
								Rectangle parent = boxes[parentId];
								int H = parent.height;
								int W = parent.width;
								int xOffset = 0;
								int yOffset = 0;
								if(j % 2 == 0) { // vertical
									H = (int) Math.floor(2*cluster.getAAF(n.getLeafSampleId())*aBOX/W);
									if(parent.height < H) { // possible due to overflow
										H = parent.height; 
									} 
									yOffset = boxOffsetY[parentId];
									boxOffsetY[parentId] += H;
								} else { // horizontal
									W = (int) Math.floor(2*cluster.getAAF(n.getLeafSampleId())*aBOX/H);
									if(parent.width < W) { // possible due to VAF overflow
										W = parent.width; 
									} 
									xOffset = boxOffsetX[parentId];
									boxOffsetX[parentId] += W;
								}
								// truncate if cumulative size larger than box(due to VAF overflow)
								if(parent.x + xOffset + W >= X + wBOX) {
									W = X + wBOX - (parent.x + xOffset);
								}
								if(parent.y + yOffset + H >= Y + hBOX) {
									H = Y + hBOX - (parent.y + yOffset);
								}	
								r = new Rectangle(parent.x + xOffset, parent.y + yOffset, W, H);
							}
							boxes[cluster.getNodeId()] = r;
							
							// draw
							g.setColor(nodeColors[cluster.getNodeId()]);
							g.fillRect(r.x, r.y, r.width, r.height);
							g.setColor(Color.black);
							g.drawLine(r.x, r.y, r.x, r.y + r.height);
							g.drawLine(r.x, r.y, r.x + r.width, r.y);
							g.drawLine(r.x + r.width, r.y, r.x + r.width, r.y + r.height);
							g.drawLine(r.x, r.y + r.height, r.x + r.width, r.y + r.height);
							
							if(subclone_sample_highligth && nodeSelected[cluster.getNodeId()]) {
								highlight.add(r);
							}
						}
					}
				}
				if(subclone_sample_highligth) {
					g.setColor(ColorPalette.HIGHLIGHT_COLOR);
					for(int i = 0; i < highlight.size(); i++) {
						Rectangle r = highlight.get(i);
						g.fillRect(r.x, r.y, r.width, r.height);
						if(r.width == 0) {
							// possible when the sibling covers the entire parent (due to overflows)
							g.fillRect(r.x-2, r.y, r.width+2, r.height);
						}if(r.height == 0) {
							g.fillRect(r.x, r.y-2, r.width, r.height+2);
						}
					}
				}
			}
		};
		return nodeIcons[nid];
	}
}
