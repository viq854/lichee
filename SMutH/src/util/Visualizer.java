package util;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.TextArea;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import lineage.PHYNetwork;
import lineage.PHYNode;
import lineage.PHYTree;
import lineage.SNVEntry;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.functors.ConstantTransformer;

import edu.uci.ics.jung.algorithms.layout.FRLayout;
import edu.uci.ics.jung.algorithms.layout.TreeLayout;
import edu.uci.ics.jung.graph.DelegateTree;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.Forest;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Context;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import edu.uci.ics.jung.visualization.control.ModalGraphMouse;
import edu.uci.ics.jung.visualization.decorators.EdgeShape;
import edu.uci.ics.jung.visualization.picking.PickedState;
import edu.uci.ics.jung.visualization.renderers.Renderer.VertexLabel.Position;

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
		graphMouse.setMode(ModalGraphMouse.Mode.PICKING);
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
		                		s += snv.getAltChar() + "/" + snv.getRefChar();
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
			String fileOutputName, final HashMap<Integer, PHYNode> nodeInfo, final PHYTree t, final PHYNetwork net, final String[] sampleNames) {	
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		TreeLayout<Integer, Integer> treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree,100,70);
		final VisualizationViewer<Integer, Integer> visServer = new VisualizationViewer<Integer, Integer>(treeLayout);
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
		                		if(n.getSNVGroup() == null) return;
		                		snvs = n.getSNVs(n.getSNVGroup().getSNVs());
		                	} else {
		                		if(n.getSNVGroup() == null) return;
		                		snvs = n.getSNVs(snvsByTag.get(n.getSNVGroup().getTag()));
		                	}
		                	for(SNVEntry snv : snvs) {
		                		s += snv.getChromosome() + " ";
		                		s += snv.getPosition() + " ";
		                		s += snv.getAltChar() + "/" + snv.getRefChar() + " ";
		                		s += ((snv.isInCNVRegion()) ? "*cnv* " : "");
		                		s += ((snv.getAnnotation() != null) ? snv.getAnnotation() : "");
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
		final JFrame frame = new JFrame("Top Lineage Tree");
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		
		JButton kill = new JButton();
		kill.setText("Kill");
		kill.setBackground(Color.red);
		kill.setOpaque(true);
		kill.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent evt) {
				Object[] pickedNodes = new HashSet<Integer>(visServer.getPickedVertexState().getPicked()).toArray();
				if(pickedNodes.length != 1) {
					JOptionPane.showMessageDialog(frame, "One node needs to be selected.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				int selectedNode = (Integer)pickedNodes[0];
				if(selectedNode < 0) {
					JOptionPane.showMessageDialog(frame, "Sample nodes cannot be removed.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				PHYNode n = nodeInfo.get(selectedNode);
				if(n.isRoot()) {
					JOptionPane.showMessageDialog(frame, "Cannot remove the root germline node.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				PHYNetwork constrNetwork = net.removeNode(n);
				ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees(); 
				constrNetwork.evaluateLineageTrees();
				if(spanningTrees.size() > 0) {
					System.out.println("Best tree error score: " + spanningTrees.get(0).getErrorScore());
					constrNetwork.displayTree(spanningTrees.get(0), sampleNames, null, null);
				} else {
					JOptionPane.showMessageDialog(frame, "No valid lineage tree was found.", "Message", JOptionPane.PLAIN_MESSAGE);
				}
			}
		});
		
		JButton collapse = new JButton();
		collapse.setText("Collapse");
		collapse.setBackground(Color.blue);
		collapse.setOpaque(true);
		collapse.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent evt) {
				Object[] pickedNodes = new HashSet<Integer>(visServer.getPickedVertexState().getPicked()).toArray();
				if(pickedNodes.length != 2) {
					JOptionPane.showMessageDialog(frame, "Two nodes need to be selected.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				PHYNode n1 = nodeInfo.get((Integer)pickedNodes[0]);
				PHYNode n2 = nodeInfo.get((Integer)pickedNodes[1]);
				
				// can collapse only clusters
				if(!n1.getSNVGroup().equals(n2.getSNVGroup())) {
					JOptionPane.showMessageDialog(frame, "Cannot collapse nodes that are not clusters of the same group.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				
				PHYNetwork constrNetwork = net.collapseClusterNodes(n1, n2);
				ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees();  
				constrNetwork.evaluateLineageTrees();
				if(spanningTrees.size() > 0) {
					constrNetwork.displayTree(spanningTrees.get(0), sampleNames, null, null);
				} else {
					JOptionPane.showMessageDialog(frame, "No valid lineage tree was found.", "Message", JOptionPane.PLAIN_MESSAGE);
				}
			}
		});
		
		GridBagConstraints c = new GridBagConstraints();
		//c.weighty = 0.8;
		c.fill = GridBagConstraints.BOTH;
		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 1.0;
		c.weighty = 1.0;
		c.ipady = 100;
		//c.anchor = GridBagConstraints.PAGE_START;
		panel.add(visServer, c);

		c.ipady = 0;
		c.fill = GridBagConstraints.CENTER;
		c.weighty = 0;
		c.weightx = 0;
		c.gridx = 0;
		c.gridy = 1;
		//int h = console.getPreferredSize().height;
		//console.setMinimumSize(new Dimension(console.getPreferredSize().width, h/4));
		//console.setPreferredSize(new Dimension(console.getPreferredSize().width, h/4));
		//console.setMaximumSize(new Dimension(console.getPreferredSize().width, h/4));
		panel.add(console, c);
		
		c.fill = GridBagConstraints.CENTER;
		c.gridx = 0;
		c.gridy = 2;
		kill.setMargin(new Insets(0, 0, 1, 0));
		panel.add(kill, c);
		
		c.fill = GridBagConstraints.CENTER;
		c.gridx = 0;
		c.gridy = 3;
		panel.add(collapse, c);
		
		frame.setContentPane(panel);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		
		
		
		//Dimension size = frame.getSize();
	    //BufferedImage image = (BufferedImage)frame.createImage(size.width, size.height);
	    //Graphics gr = image.getGraphics();
	    //frame.paint(gr);
	    //gr.dispose();
      
	    if(fileOutputName != null) {
	    	//try {
	    		//ImageIO.write(image, "jpg", new File(fileOutputName));
	    	//}
	    	//catch (IOException e) {
	    		//e.printStackTrace();
	    		//System.err.println("Failed to save tree to the file: " + fileOutputName);
	    //	}
	    }
	}
}
