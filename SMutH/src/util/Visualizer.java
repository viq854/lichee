package util;

import ppt.*;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Paint;
import java.awt.Panel;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.TextArea;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JToggleButton;

import lineage.PHYNode;
import lineage.PHYTree;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.functors.ConstantTransformer;

import edu.uci.ics.jung.algorithms.layout.DAGLayout;
import edu.uci.ics.jung.algorithms.layout.FRLayout;
import edu.uci.ics.jung.algorithms.layout.PolarPoint;
import edu.uci.ics.jung.algorithms.layout.RadialTreeLayout;
import edu.uci.ics.jung.algorithms.layout.TreeLayout;
import edu.uci.ics.jung.graph.DelegateTree;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.Forest;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.visualization.Layer;
import edu.uci.ics.jung.visualization.VisualizationServer;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import edu.uci.ics.jung.visualization.control.ModalGraphMouse;
import edu.uci.ics.jung.visualization.decorators.AbstractEdgeShapeTransformer;
import edu.uci.ics.jung.visualization.decorators.EdgeShape;
import edu.uci.ics.jung.visualization.layout.LayoutTransition;
import edu.uci.ics.jung.visualization.picking.PickedState;
import edu.uci.ics.jung.visualization.renderers.Renderer.VertexLabel.Position;
import edu.uci.ics.jung.visualization.util.Animator;
import edu.uci.ics.jung.graph.util.Context;

public class Visualizer {

	
	/**
	 * Displays the constraint network (used by the lineage package)
	 */
	public static void showNetwork(DirectedGraph<Integer, Integer> g, HashMap<Integer, String> nodeLabels) {
		JFrame frame = new JFrame("Constraint Network");		
		FRLayout<Integer, Integer> dagLayout = new FRLayout<Integer, Integer>(g);
		dagLayout.setRepulsionMultiplier(0.4);
		VisualizationViewer<Integer, Integer> visServer = new VisualizationViewer<Integer, Integer>(dagLayout);
		visServer.setPreferredSize(new Dimension(600, 500));
		
		// mouse
		DefaultModalGraphMouse<Integer, Integer> graphMouse = new DefaultModalGraphMouse<Integer, Integer>();
		graphMouse.setMode(ModalGraphMouse.Mode.TRANSFORMING);
		visServer.setGraphMouse(graphMouse);
		
		// node labels
		final HashMap<Integer, String> nodeLabelsFinal = new HashMap<Integer, String>(nodeLabels);
		Transformer<Integer, String> lt = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) {
					return nodeLabelsFinal.get(num);
				}
				else return null;
			}
		};
		
		// node colors
		Transformer<Integer, Paint> ct = new Transformer<Integer, Paint>() {
			public Paint transform(Integer num){
				return new Color(178, 34, 34);
			}
		};
		
		@SuppressWarnings({ "unchecked", "rawtypes" })
		Transformer<Integer,Font> vt = new ConstantTransformer(new Font("Helvetica", Font.BOLD, 12));
		visServer.getRenderContext().setVertexFontTransformer(vt);
		visServer.getRenderContext().setVertexFillPaintTransformer(ct);
		visServer.getRenderContext().setVertexLabelTransformer(lt);
		Transformer<Context<Graph<Integer,Integer>,Integer>,Shape> et = new EdgeShape.Line<Integer,Integer>();		
		visServer.getRenderContext().setEdgeShapeTransformer(et);
		
		// content
		Container content = frame.getContentPane();
		content.add(visServer);
		JPanel controls = new JPanel();
		content.add(controls, BorderLayout.SOUTH);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		
		Dimension size = frame.getSize();
	    BufferedImage image = (BufferedImage)frame.createImage(size.width, size.height);
	    Graphics gr = image.getGraphics();
	    frame.paint(gr);
	    gr.dispose();
	}
	
	/**
	 * Displays the lineage tree
	 * The tree nodes can be interacted with to obtain more information
	 */
	public static void showLineageTreeBreakdown(DirectedGraph<Integer, Integer> g, final HashMap<Integer, String> nodeLabels, 
			final HashMap<String, ArrayList<SNVEntry>> snvsByTag, 
			String fileOutputName, final HashMap<Integer, PHYNode> nodeInfo, final PHYTree t) {	
		
		// create a color for each internal node
		Random rand = new Random();
		final HashMap<Integer, Color> nodeColors = new HashMap<Integer, Color>();
		for(Integer n : nodeLabels.keySet()) {
			if(n < 0) {
				// leaf
				nodeColors.put(n, new Color(0, 191, 255));
			} else {
				nodeColors.put(n, new Color(rand.nextFloat(), rand.nextFloat(), rand.nextFloat()).brighter());
			}
		}
		
		
		JFrame frame = new JFrame("Best Lineage Tree");
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		TreeLayout<Integer, Integer> treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree,100,70);
		VisualizationViewer<Integer, Integer> visServer = new VisualizationViewer<Integer, Integer>(treeLayout);
		visServer.setPreferredSize(new Dimension(1500, 600));
		
		DefaultModalGraphMouse<Integer, Integer> graphMouse = new DefaultModalGraphMouse<Integer, Integer>();
		graphMouse.setMode(ModalGraphMouse.Mode.PICKING);
		visServer.setGraphMouse(graphMouse);
		
		// node click listeners
		JPanel console = new JPanel();
		final TextArea info = new TextArea();
		info.setEditable(false);
		console.add(info);
		final PickedState<Integer> pickedState = visServer.getPickedVertexState();
		pickedState.addItemListener(new ItemListener() { 
		    @Override
		    public void itemStateChanged(ItemEvent e) {
		    Object o = e.getItem();
		        if (o instanceof Integer) {
		            Integer node = (Integer) o;
		            if (pickedState.isPicked(node)) {
		            	PHYNode n = nodeInfo.get(node);
		                if(node < 0) { // sample 
		                	info.setText(t.getLineage(n.getLeafSampleId(), nodeLabels.get(node)));
		                	
		                } else {
		                	String s = n.getLongLabel();
		                	s += "\nSNVs: \n";
		                	ArrayList<SNVEntry> snvs;
		                	if(snvsByTag == null) {
		                		snvs = n.getSNVs(n.getSNVGroup().getSNVs());
		                	} else {
		                		snvs = n.getSNVs(snvsByTag.get(n.getSNVGroup().getTag()));
		                	}
		                	for(SNVEntry snv : snvs) {
		                		s += snv.getChromosome() + " ";
		                		s += snv.getPosition() + " ";
		                		s += snv.alt + "/" + snv.ref;
		                		s += "\n";
		                	}
		                	info.setText(s);
		                }
		            }
		        }
		    }
		});
		
		// node labels
		Transformer<Integer, String> vlt = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabels != null && nodeLabels.containsKey(num)) {
					return nodeLabels.get(num);
				}
				else return null;
			}
		};
		// node colors
		Transformer<Integer, Paint> vpt = new Transformer<Integer, Paint>() {
			public Paint transform(Integer num){
				if (nodeLabels != null && nodeLabels.containsKey(num)) {
					return nodeColors.get(num);
				} 
				else return Color.ORANGE;
			}
		};
		// node shapes
		Transformer<Integer, Shape> vts = new Transformer<Integer, Shape>() {
			public Shape transform(Integer num){
				if(num < 0) {
					Rectangle2D r = new Rectangle2D.Float(-10, -10, 20, 20);
					Rectangle2D r2 = new Rectangle2D.Float(-10, -10, 5, 5);
			
					Area a = new Area(r);
					
					a.add(new Area(r2));
					return a;
				} 
				return new Ellipse2D.Float(-10, -10, 20, 20);
			}
		};
		
		Transformer<Context<Graph<Integer,Integer>,Integer>,Shape> et = new EdgeShape.Line<Integer,Integer>();		
		@SuppressWarnings({ "unchecked", "rawtypes" })
		Transformer<Integer,Font> vt = new ConstantTransformer(new Font("Helvetica", Font.BOLD, 12));
		
		visServer.getRenderer().getVertexLabelRenderer().setPosition(Position.S);
		visServer.getRenderContext().setVertexShapeTransformer(vts);
		visServer.getRenderContext().setVertexFontTransformer(vt);
		visServer.getRenderContext().setEdgeShapeTransformer(et);
		visServer.getRenderContext().setVertexLabelTransformer(vlt);
		visServer.getRenderContext().setVertexFillPaintTransformer(vpt);
		
		// frame content
		Container content = frame.getContentPane();
		content.add(visServer);
		content.add(console, BorderLayout.SOUTH);
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		
		Dimension size = frame.getSize();
	    BufferedImage image = (BufferedImage)frame.createImage(size.width, size.height);
	    Graphics gr = image.getGraphics();
	    frame.paint(gr);
	    gr.dispose();
      
	    if(fileOutputName != null) {
	    	try {
	    		ImageIO.write(image, "jpg", new File(fileOutputName));
	    	}
	    	catch (IOException e) {
	    		e.printStackTrace();
	    		System.err.println("Failed to save tree to the file: " + fileOutputName);
	    	}
	    }
	}
	
	/**
	 * Displays the lineage tree
	 * The tree nodes can be interacted with to obtain more information
	 */
	public static void showLineageTree(DirectedGraph<Integer, Integer> g, final HashMap<Integer, String> nodeLabels, 
			final HashMap<String, ArrayList<SNVEntry>> snvsByTag, 
			String fileOutputName, final HashMap<Integer, PHYNode> nodeInfo, final PHYTree t) {	
		JFrame frame = new JFrame("Best Lineage Tree");
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		TreeLayout<Integer, Integer> treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree,100,70);
		VisualizationViewer<Integer, Integer> visServer = new VisualizationViewer<Integer, Integer>(treeLayout);
		visServer.setPreferredSize(new Dimension(1500, 600));
		
		DefaultModalGraphMouse<Integer, Integer> graphMouse = new DefaultModalGraphMouse<Integer, Integer>();
		graphMouse.setMode(ModalGraphMouse.Mode.PICKING);
		visServer.setGraphMouse(graphMouse);
		
		// node click listeners
		JPanel console = new JPanel();
		final TextArea info = new TextArea();
		info.setEditable(false);
		console.add(info);
		final PickedState<Integer> pickedState = visServer.getPickedVertexState();
		pickedState.addItemListener(new ItemListener() { 
		    @Override
		    public void itemStateChanged(ItemEvent e) {
		    Object o = e.getItem();
		        if (o instanceof Integer) {
		            Integer node = (Integer) o;
		            if (pickedState.isPicked(node)) {
		            	PHYNode n = nodeInfo.get(node);
		                if(node < 0) { // sample 
		                	info.setText(t.getLineage(n.getLeafSampleId(), nodeLabels.get(node)));
		                	
		                } else {
		                	String s = n.getLongLabel();
		                	s += "\nSNVs: \n";
		                	ArrayList<SNVEntry> snvs;
		                	if(snvsByTag == null) {
		                		snvs = n.getSNVs(n.getSNVGroup().getSNVs());
		                	} else {
		                		snvs = n.getSNVs(snvsByTag.get(n.getSNVGroup().getTag()));
		                	}
		                	for(SNVEntry snv : snvs) {
		                		s += snv.getChromosome() + " ";
		                		s += snv.getPosition() + " ";
		                		s += snv.alt + "/" + snv.ref;
		                		s += "\n";
		                	}
		                	info.setText(s);
		                }
		            }
		        }
		    }
		});
		
		// node labels
		Transformer<Integer, String> vlt = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabels != null && nodeLabels.containsKey(num)) {
					return nodeLabels.get(num);
				}
				else return null;
			}
		};
		// node colors
		Transformer<Integer, Paint> vpt = new Transformer<Integer, Paint>() {
			public Paint transform(Integer num){
				if (nodeLabels != null && nodeLabels.containsKey(num)) {
					if(num < 0) {
						return new Color(0, 191, 255);
					} else {
						return new Color(60, 179, 113);
					}
				} 
				else return Color.ORANGE;
			}
		};
		// node shapes
		Transformer<Integer, Shape> vts = new Transformer<Integer, Shape>() {
			public Shape transform(Integer num){
				if(num < 0) {
					return new Rectangle2D.Float(-10, -10, 20, 20);
				} 
				return new Ellipse2D.Float(-10, -10, 20, 20);
			}
		};
		
		Transformer<Context<Graph<Integer,Integer>,Integer>,Shape> et = new EdgeShape.Line<Integer,Integer>();		
		@SuppressWarnings({ "unchecked", "rawtypes" })
		Transformer<Integer,Font> vt = new ConstantTransformer(new Font("Helvetica", Font.BOLD, 12));
		
		visServer.getRenderer().getVertexLabelRenderer().setPosition(Position.S);
		visServer.getRenderContext().setVertexShapeTransformer(vts);
		visServer.getRenderContext().setVertexFontTransformer(vt);
		visServer.getRenderContext().setEdgeShapeTransformer(et);
		visServer.getRenderContext().setVertexLabelTransformer(vlt);
		visServer.getRenderContext().setVertexFillPaintTransformer(vpt);
		
		// frame content
		Container content = frame.getContentPane();
		content.add(visServer);
		content.add(console, BorderLayout.SOUTH);
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		
		Dimension size = frame.getSize();
	    BufferedImage image = (BufferedImage)frame.createImage(size.width, size.height);
	    Graphics gr = image.getGraphics();
	    frame.paint(gr);
	    gr.dispose();
      
	    if(fileOutputName != null) {
	    	try {
	    		ImageIO.write(image, "jpg", new File(fileOutputName));
	    	}
	    	catch (IOException e) {
	    		e.printStackTrace();
	    		System.err.println("Failed to save tree to the file: " + fileOutputName);
	    	}
	    }
	}
	
	public static void TreeVisualizer(DirectedGraph<Integer, Integer> g, HashMap<Integer, String> nodeLabels) {
		
		VisualizationViewer<Integer, Integer> visServer;
		TreeLayout<Integer, Integer> treeLayout;
		RadialTreeLayout<Integer, Integer> radialLayout;
		
		final HashMap<Integer, String> nodeLabelsFinal = new HashMap<Integer, String>(nodeLabels);
		
		JFrame frame = new JFrame(Configs.testName+" Tree View");
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		//FRLayout<Integer, Integer> tree = new FRLayout<Integer, Integer>(g);//DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree,100,70);
		visServer = new VisualizationViewer<Integer, Integer>(treeLayout);
		visServer.setPreferredSize(new Dimension(1500, 600));
		
		//Creating graph mouse
		DefaultModalGraphMouse graphMouse = new DefaultModalGraphMouse();
		graphMouse.setMode(ModalGraphMouse.Mode.TRANSFORMING);
		visServer.setGraphMouse(graphMouse);
		
		Transformer<Integer, String> vlt = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) {
					return nodeLabelsFinal.get(num);
				}
				else return null;
			}
		};
		
		Transformer<Integer, Paint> vpt = new Transformer<Integer, Paint>() {
			public Paint transform(Integer num){
				if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) {
					if(num < 0) {
						return new Color(0, 191, 255);
					} else {
						return new Color(60, 179, 113);
					}
				} 
				else return Color.ORANGE;
			}
		};
		
		Transformer<Integer, Shape> vts = new Transformer<Integer, Shape>() {
			public Shape transform(Integer num){
				if(num < 0) {
					return new Rectangle2D.Float(-10, -10, 20, 20);
				} 
				return new Ellipse2D.Float(-10, -10, 20, 20);
			}
		};
		
		Transformer<Context<Graph<Integer,Integer>,Integer>,Shape> et = new EdgeShape.Line<Integer,Integer>();		
		Transformer<Integer,Font> vt = new ConstantTransformer(new Font("Helvetica", Font.BOLD, 12));
		
		visServer.getRenderer().getVertexLabelRenderer().setPosition(Position.S);
		visServer.getRenderContext().setVertexShapeTransformer(vts);
		visServer.getRenderContext().setVertexFontTransformer(vt);
		visServer.getRenderContext().setEdgeShapeTransformer(et);
		visServer.getRenderContext().setVertexLabelTransformer(vlt);
		visServer.getRenderContext().setVertexFillPaintTransformer(vpt);
		//visServer.getRenderContext().s
		
		
		
		Container content = frame.getContentPane();
		content.add(visServer);
		
		JPanel controls = new JPanel();
		content.add(controls, BorderLayout.SOUTH);
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		
		Dimension size = frame.getSize();
	    BufferedImage image = (BufferedImage)frame.createImage(size.width, size.height);
	    Graphics gr = image.getGraphics();
	    frame.paint(gr);
	    gr.dispose();
      
	    
	    /**
	     * Raheleh: saving the tree
	     */
	    try {
	    	ImageIO.write(image, "jpg", new File(Configs.path+Configs.testName+".network.jpg"));
	    }
	    catch (IOException e) {
	    	e.printStackTrace();
	    }
	}
	
	
	/**
	 * Deprecated Function
	 * Creates the actual graphics and sets up visualizing the tree
	 * 
	 * Uses a combination of JUNG and Java's Swing libraries to display
	 * tree. This might be made into it's own class if necessary.
	 * 
	 * @param g	A specific instance of a DirectedGraph built already
	 * @param hashMap 
	 */
	public static void TreeVisualizer(DirectedGraph<Integer, Integer> g, HashMap<Integer, String> nodeLabels, HashMap<Integer, String> edgeLabels, SNVDatabase db) {
		
		
		VisualizationViewer<Integer, Integer> visServer;
		TreeLayout<Integer, Integer> treeLayout;
		RadialTreeLayout<Integer, Integer> radialLayout;
		VisualizationServer.Paintable rings;
		
		final HashMap<Integer, String> nodeLabelsFinal;
		final HashMap<Integer, String> edgeLabelsFinal;
		
		if (nodeLabels == null) nodeLabelsFinal = new HashMap<Integer, String>();
		else nodeLabelsFinal = new HashMap<Integer, String>(nodeLabels);
		
		nodeLabelsFinal.put(0, "germline");
		for(int i=0; i < db.getNumofSamples(); i++){
			nodeLabelsFinal.put(-i-1, db.getName(i));
		}
		
		
		if (edgeLabels == null) edgeLabelsFinal = new HashMap<Integer, String>();
		else edgeLabelsFinal = new HashMap<Integer, String>(edgeLabels);
		
		JFrame frame = new JFrame(Configs.testName+" Tree View");
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree,100,70);
		//treeLayout.setSize(new Dimension(600, 600));
		radialLayout = new RadialTreeLayout<Integer, Integer>(tree);
		radialLayout.setSize(new Dimension(700, 700));
//		BasicVisualizationServer<Integer, Integer> visServer =
//			new BasicVisualizationServer<Integer, Integer>(treeLayout);
		visServer = new VisualizationViewer<Integer, Integer>(treeLayout);
		visServer.setPreferredSize(new Dimension(600, 500));
		
		//Setting up transformers for JUNG
		
		Transformer<Integer, String> PhyVertexLabelTransformer = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) return nodeLabelsFinal.get(num);
				/*if (num < 0)
					return db.getName(-1*num - 1);
				else if (num == 0)
					return "Root";*/
				//else return "N-" + num.toString();
				else return null;//"N-" + num.toString();
			}
		};
		
		Transformer<Integer, String> PhyEdgeLabelTransformer = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (edgeLabelsFinal != null && edgeLabelsFinal.containsKey(num)) return edgeLabelsFinal.get(num); 
				else
					return null;
			}
		};
		
//		Transformer<Integer, EdgeShape<Integer, Integer>> PhyEdgeShapeTransformer = new Transformer<Integer, EdgeShape<Integer, Integer>>(){
//			public Line<Integer, Integer> transformer(Integer num){
//				if (num >= 0) return (new EdgeShape.Line<Integer, Integer>());
//				else return (new EdgeShape.BentLine<Integer,Integer>());
//			}
//
//			@Override
//			public EdgeShape<Integer, Integer> transform(Integer num) {
//				if (num >= 0) return (new EdgeShape.Line<Integer, Integer>());
//				else return (new EdgeShape.BentLine<Integer,Integer>());
//			}
//		}
		
		float dash[] = {5.0f};
		final Stroke edgeStroke = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,
				BasicStroke.JOIN_MITER, 10.0f, dash, 0.0f);
		Transformer<Integer, Stroke> PhyEdgeStrokeTransformer =
			new Transformer<Integer, Stroke>() {
				public Stroke transform(Integer num) {
					if (num < 0) return edgeStroke;
					else return null;
				}
			};
			
		Transformer<Integer, Paint> PhyVertexPaintTransformer =
			new Transformer<Integer, Paint>() {
				public Paint transform(Integer num){
					if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) return Color.BLUE;
					if (num < 0) return Color.GREEN;
					else return Color.RED;
				}
			};
		
		//visServer.getRenderContext().setVertexLabelTransformer(new ToStringLabeller<Integer>());
		//visServer.setBackground(Color.WHITE);
		visServer.getRenderContext().setVertexLabelTransformer(PhyVertexLabelTransformer);
		visServer.getRenderContext().setVertexFillPaintTransformer(PhyVertexPaintTransformer);
		visServer.getRenderContext().setEdgeShapeTransformer(new EdgeShape.Line<Integer, Integer>());
		visServer.getRenderContext().setEdgeLabelTransformer(PhyEdgeLabelTransformer);
		visServer.getRenderContext().setEdgeStrokeTransformer(PhyEdgeStrokeTransformer);
		
		//Creating graph mouse
		DefaultModalGraphMouse graphMouse = new DefaultModalGraphMouse();
		graphMouse.setMode(ModalGraphMouse.Mode.TRANSFORMING);
		visServer.setGraphMouse(graphMouse);
		
		Container content = frame.getContentPane();
		//final GraphZoomScrollPane panel = new GraphZoomScrollPane(visServer);
		content.add(visServer);
		
		JToggleButton radial = new JToggleButton("Radial");
		
		
		JPanel controls = new JPanel();
		//controls.setBackground(Color.WHITE);
		controls.add(radial);
		content.add(controls, BorderLayout.SOUTH);
		
		//JFrame frame = new JFrame("Simple Tree View");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		frame.getContentPane().add(content);
		frame.pack();
		frame.setVisible(true);
		
		Dimension size = frame.getSize();
	      //BufferedImage image = new BufferedImage(size.width, size.height, BufferedImage.TYPE_INT_RGB);
	      BufferedImage image = (BufferedImage)frame.createImage(size.width, size.height);
	      Graphics gr = image.getGraphics();
	      frame.paint(gr);
	      gr.dispose();
	      
	      try
	      {
	        ImageIO.write(image, "jpg", new File(Configs.path+Configs.testName+".tree.jpg"));
	      }
	      catch (IOException e)
	      {
	        e.printStackTrace();
	      }
	}

}
