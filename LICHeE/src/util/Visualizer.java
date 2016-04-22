package util;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.UIManager;

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
import edu.uci.ics.jung.visualization.decorators.PickableVertexPaintTransformer;
import edu.uci.ics.jung.visualization.picking.PickedState;
import edu.uci.ics.jung.visualization.renderers.Renderer.VertexLabel.Position;
import de.erichseifert.vectorgraphics2d.PDFGraphics2D;

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

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void showLineageTree(DirectedGraph<Integer, Integer> g, final HashMap<Integer, String> nodeLabels, 
			final HashMap<String, ArrayList<SNVEntry>> snvsByTag, 
			String fileOutputName, final HashMap<Integer, PHYNode> nodeInfo, 
			final PHYTree t, final PHYNetwork net, final ArrayList<String> sampleNames) {	
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(g);
		tree.setRoot(0);
		TreeLayout<Integer, Integer> treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree, 100, 70);
		final VisualizationViewer<Integer, Integer> visServer = new VisualizationViewer<Integer, Integer>(treeLayout);	
		DefaultModalGraphMouse<Integer, Integer> graphMouse = new DefaultModalGraphMouse<Integer, Integer>();
		graphMouse.setMode(ModalGraphMouse.Mode.PICKING);
		visServer.setGraphMouse(graphMouse);
		
		try {
			UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
		} catch (Exception e1) {
			System.exit(0);
		} 
		
		JPanel licheePanel = new JPanel();
		final JPanel graphPanel = new javax.swing.JPanel();
		final JTextArea consoleTextArea = new JTextArea();
		javax.swing.JScrollPane consolePanel = new javax.swing.JScrollPane(consoleTextArea);
		consolePanel.setFont(new java.awt.Font("Arial Unicode MS", 0, 13));
        JButton removeButton = new javax.swing.JButton();
        JButton collapseButton = new javax.swing.JButton();
        JButton snapshotButton = new javax.swing.JButton();

        // frame content
     	final JFrame frame = new JFrame("Top Lineage Tree");
     	frame.setBackground(Color.WHITE);
     	frame.setContentPane(licheePanel);
     		
        licheePanel.setBackground(new java.awt.Color(255, 255, 255));
        graphPanel.setBackground(new java.awt.Color(211, 255, 248));
        graphPanel.setLayout(new BorderLayout());
        graphPanel.add(visServer, BorderLayout.PAGE_START);
        visServer.setBackground(graphPanel.getBackground());
        
        removeButton.setBackground(new java.awt.Color(255, 255, 255));
        removeButton.setFont(new java.awt.Font("Arial Unicode MS", 0, 13));
        removeButton.setText("Remove");
        collapseButton.setBackground(new java.awt.Color(255, 255, 255));
        collapseButton.setFont(new java.awt.Font("Arial Unicode MS", 0, 13));
        collapseButton.setText("Collapse");
        snapshotButton.setBackground(new java.awt.Color(255, 255, 255));
        snapshotButton.setFont(new java.awt.Font("Arial Unicode MS", 0, 13));
        snapshotButton.setText("Snapshot");
        
        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(licheePanel);
        licheePanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(graphPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(consolePanel, javax.swing.GroupLayout.PREFERRED_SIZE, 500, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(removeButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(collapseButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(snapshotButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(451, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(graphPanel)//, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(removeButton, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(collapseButton, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(snapshotButton, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 47, Short.MAX_VALUE))
                    .addComponent(consolePanel))
                .addContainerGap())
        );
		
		// node click listeners
		consoleTextArea.setEditable(false);
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
		                	consoleTextArea.setText(t.getLineage(n.getLeafSampleId(), nodeLabels.get(node)));		                	
		                } else {
		                	String s = n.getLongLabel();
		                	s += "\nSSNV members:\n";
		                	ArrayList<SNVEntry> snvs;
		                	if(snvsByTag == null) {
		                		if(n.getSNVGroup() == null) {
		                			consoleTextArea.setText(""); 
		                			return;
		                		}
		                		snvs = n.getSNVs(n.getSNVGroup().getSNVs());
		                	} else {
		                		if(n.getSNVGroup() == null) {
		                			consoleTextArea.setText(""); 
		                			return;
		                		}
		                		snvs = n.getSNVs(snvsByTag.get(n.getSNVGroup().getTag()));
		                	}
		                	for(SNVEntry snv : snvs) {
		                		s += snv.getChromosome() + " ";
		                		s += snv.getPosition() + " ";
		                		s += snv.getDescription() + "\n";
		                	}
		                	consoleTextArea.setText(s);
		                }
		            }
		        }
		    }
		});

		// graph node/edge properties
		visServer.getRenderer().getVertexLabelRenderer().setPosition(Position.S);
		visServer.getRenderContext().setVertexShapeTransformer(new Transformer<Integer, Shape>() {
			public Shape transform(Integer num) {
				if(num < 0) {
					return new Rectangle2D.Float(-10, -10, 20, 20);
				}
				PHYNode n = nodeInfo.get(num);
                int size = n.getSize();
                if(size == 0) size = 3;
                int w = 10*size;
                if(w > 50) w = 50;
				return new Ellipse2D.Float(-w/2, -w/2, w, w); } });
		visServer.getRenderContext().setVertexFontTransformer(new ConstantTransformer(new Font("Arial Unicode MS", Font.PLAIN, 12)));
		
		visServer.getRenderContext().setEdgeShapeTransformer(new EdgeShape.QuadCurve<Integer, Integer>());
		visServer.getRenderContext().setVertexLabelTransformer(new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabels != null && nodeLabels.containsKey(num)) {
					return nodeLabels.get(num);
				}
				else return null;
			}});
		visServer.getRenderContext().setVertexFillPaintTransformer(new PickableVertexPaintTransformer<Integer>(pickedState, Color.white, Color.yellow));
		visServer.getRenderContext().setEdgeDrawPaintTransformer(new Transformer<Integer, Paint>() {
			public Paint transform(Integer s) {
				return Color.BLACK;
			}
		});
		
		// remove node button
		removeButton.addActionListener(new ActionListener() {
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
		
		// collapse button
		collapseButton.addActionListener(new ActionListener() {
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
	
		// snapshot button
		snapshotButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent evt) {
				//BufferedImage bufImage = new BufferedImage(graphPanel.getSize().width, graphPanel.getSize().height, BufferedImage.TYPE_INT_RGB);
				//graphPanel.paint(bufImage.createGraphics());
				PDFGraphics2D bufImage = new PDFGraphics2D(0.0, 0.0, graphPanel.getSize().width, graphPanel.getSize().height);
				graphPanel.paint(bufImage);				
				
				final JFileChooser fc = new JFileChooser();
				int returnVal = fc.showSaveDialog(graphPanel);
		        if (returnVal == JFileChooser.APPROVE_OPTION) {
		            File imageFile = fc.getSelectedFile();
		            try{
		            	imageFile.createNewFile();
		            	//ImageIO.write(bufImage, "jpeg", imageFile);
		            	Files.write(imageFile.toPath(), bufImage.getBytes());
		            }catch(Exception ex){}
		        }
			}
		});
		
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
}
