package util;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Paint;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.swing.BorderFactory;
import javax.swing.GroupLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.border.TitledBorder;

import lineage.PHYNetwork;
import lineage.PHYNode;
import lineage.PHYTree;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.functors.ConstantTransformer;

import de.erichseifert.vectorgraphics2d.PDFGraphics2D;
import edu.uci.ics.jung.algorithms.layout.FRLayout;
import edu.uci.ics.jung.algorithms.layout.Layout;
import edu.uci.ics.jung.algorithms.layout.TreeLayout;
import edu.uci.ics.jung.graph.DelegateTree;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.Forest;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Context;
import edu.uci.ics.jung.visualization.Layer;
import edu.uci.ics.jung.visualization.RenderContext;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.CrossoverScalingControl;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import edu.uci.ics.jung.visualization.control.ModalGraphMouse;
import edu.uci.ics.jung.visualization.control.ScalingControl;
import edu.uci.ics.jung.visualization.decorators.ConstantDirectionalEdgeValueTransformer;
import edu.uci.ics.jung.visualization.decorators.EdgeShape;
import edu.uci.ics.jung.visualization.picking.PickedState;
import edu.uci.ics.jung.visualization.renderers.BasicVertexRenderer;
import edu.uci.ics.jung.visualization.renderers.DefaultVertexLabelRenderer;
import edu.uci.ics.jung.visualization.renderers.Renderer.VertexLabel.Position;
import edu.uci.ics.jung.visualization.transform.shape.GraphicsDecorator;

public class Visualizer {

	/**
	 * Displays the constraint network
	 */
	public static void showNetwork(DirectedGraph<Integer, Integer> g, HashMap<Integer, String> nodeLabels) {
		JFrame frame = new JFrame("LICHeE Evolutionary Constraint Network Viewer");		
		FRLayout<Integer, Integer> dagLayout = new FRLayout<Integer, Integer>(g);
		dagLayout.setRepulsionMultiplier(0.4);
		VisualizationViewer<Integer, Integer> visServer = new VisualizationViewer<Integer, Integer>(dagLayout);
		visServer.setPreferredSize(new Dimension(600, 500));
		DefaultModalGraphMouse<Integer, Integer> graphMouse = new DefaultModalGraphMouse<Integer, Integer>();
		graphMouse.setMode(ModalGraphMouse.Mode.PICKING);
		visServer.setGraphMouse(graphMouse);
		final HashMap<Integer, String> nodeLabelsFinal = new HashMap<Integer, String>(nodeLabels);
		Transformer<Integer, String> lt = new Transformer<Integer, String>(){
			public String transform(Integer num) {
				if (nodeLabelsFinal != null && nodeLabelsFinal.containsKey(num)) {
					return nodeLabelsFinal.get(num);
				}
				else return null;
			}
		};
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
	
	@SuppressWarnings({"serial"})
	public static VisualizationViewer<Integer, Integer> setupJungViewer(final DelegateTree<Integer, Integer> tree, final LineageDisplayConfig config) {
		Layout<Integer, Integer> treeLayout = new TreeLayout<Integer, Integer>((Forest<Integer, Integer>) tree, 90, 130);
		VisualizationViewer<Integer, Integer> treeView = new VisualizationViewer<Integer, Integer>(treeLayout);
		DefaultModalGraphMouse<Integer, Integer> graphMouse = new DefaultModalGraphMouse<Integer, Integer>();
		graphMouse.setMode(ModalGraphMouse.Mode.PICKING);
		treeView.setGraphMouse(graphMouse);
		treeView.setBackground(Color.WHITE);
		treeView.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder()));
		ScalingControl scaler = new CrossoverScalingControl();
		scaler.scale(treeView, 1 / 1.5f, treeView.getLocation());
		
		// ---- vertices -----
		treeView.getRenderer().getVertexLabelRenderer().setPosition(Position.S);
		treeView.getRenderContext().setVertexLabelRenderer(new DefaultVertexLabelRenderer(Color.BLACK) {
            @Override
            public <V> Component getVertexLabelRendererComponent(
                JComponent vv, Object value, Font font, 
                boolean isSelected, V vertex) { // color
                	super.getVertexLabelRendererComponent(vv, value, font, isSelected, vertex);
                	if(!isSelected) {
                		setForeground(config.getNodeLabelColor((Integer) vertex));
                	}
                	return this;
        }}); 
		treeView.getRenderer().setVertexRenderer(new BasicVertexRenderer<Integer, Integer>() {
			@Override
			protected void paintIconForVertex(RenderContext<Integer,Integer> rc, Integer v, Layout<Integer,Integer> layout) { 
		        GraphicsDecorator g = rc.getGraphicsContext(); 
		        boolean vertexHit = true; 
		        Shape shape = rc.getVertexShapeTransformer().transform(v); 
		        Point2D p = layout.transform(v); 
		        p = rc.getMultiLayerTransformer().transform(Layer.LAYOUT, p); 
		        float x = (float)p.getX(); 
		        float y = (float)p.getY(); 
		        AffineTransform xform = AffineTransform.getTranslateInstance(x,y); 
		        shape = xform.createTransformedShape(shape); 
		        vertexHit = vertexHit(rc, shape);  
		        if (vertexHit) { 
		        	if(rc.getVertexIconTransformer() != null) { 
		        		Icon icon = rc.getVertexIconTransformer().transform(v); 
		        		if(icon != null) { 
		        			paintShapeForVertex(rc, v, shape);
		        			g.draw(icon, rc.getScreenDevice(), shape, (int)x, (int)y); 
		        		} else { 
		        			paintShapeForVertex(rc, v, shape); 
		        		} 
		        	} else { 
		        		paintShapeForVertex(rc, v, shape); 
		        	} 
		        } 
		    } 
		});
		treeView.getRenderContext().setVertexFontTransformer(new Transformer<Integer, Font>() {
			public Font transform(Integer num) { // font
				return config.getNodeFont(num);
		}});	
		treeView.getRenderContext().setVertexLabelTransformer(new Transformer<Integer, String>(){
			public String transform(Integer num) { // label
				return config.getExternalNodeLabel(num);
		 }});
		treeView.getRenderContext().setVertexShapeTransformer(new Transformer<Integer, Shape>() {
			public Shape transform(Integer num) { // shape
				return config.getNodeShape(num);
		} });
		final PickedState<Integer> pickedState = treeView.getPickedVertexState();
		treeView.getRenderContext().setVertexIconTransformer(new Transformer<Integer, Icon> () {
			public Icon transform(final Integer num) {
				return config.getNodeIcon(num);
		}});
		treeView.getRenderContext().setVertexFillPaintTransformer(new Transformer<Integer, Paint>() {
			public Paint transform(Integer s) { // colors
				if(pickedState.isPicked(s)) {
					return LineageDisplayConfig.ColorPalette.SELECTED_NODE_COLOR;
				}
				return config.getNodeFillColor(s);
			}
			
		});
		treeView.getRenderContext().setVertexDrawPaintTransformer(new Transformer<Integer, Paint>() {
			public Paint transform(Integer s) {
				return config.getNodeBorderColor(s);
			}
		});
		treeView.getRenderContext().setVertexStrokeTransformer(new Transformer<Integer, Stroke>() {
			public Stroke transform(Integer s) {
				return new BasicStroke(3.0f);
			}
		});
		
		// ----- edges -----
		if(config.beautify) {
			treeView.getRenderContext().setEdgeShapeTransformer(new EdgeShape.CubicCurve<Integer, Integer>());
		} else {
			treeView.getRenderContext().setEdgeShapeTransformer(new EdgeShape.QuadCurve<Integer, Integer>());
		}
		//treeView.getRenderContext().setEdgeShapeTransformer(new EdgeShape.Line<Integer, Integer>());
		treeView.getRenderContext().getEdgeLabelRenderer().setRotateEdgeLabels(false);
		treeView.getRenderContext().setLabelOffset(20);
		treeView.getRenderContext().setEdgeLabelClosenessTransformer(new ConstantDirectionalEdgeValueTransformer<Integer, Integer>(0.5d, 0.5d));	
		treeView.getRenderContext().setArrowFillPaintTransformer(treeView.getRenderContext().getEdgeDrawPaintTransformer());
		treeView.getRenderContext().setArrowDrawPaintTransformer(treeView.getRenderContext().getEdgeDrawPaintTransformer());
		treeView.getRenderContext().setEdgeDrawPaintTransformer(new Transformer<Integer, Paint>() {
			public Paint transform(Integer s) { // color
				return config.getEdgeColor(tree.getSource(s), tree.getDest(s));
			}
		});
		treeView.getRenderContext().setEdgeStrokeTransformer(new Transformer<Integer, Stroke>() {
			public Stroke transform(Integer s) { // edge width
				if(config.beautify) {
					return new BasicStroke(14.0f);
				} else {
					return new BasicStroke(2.0f);
				}
			}
		});
		treeView.getRenderContext().setEdgeFontTransformer(new Transformer<Integer, Font>() {
			public Font transform(Integer num) { // font
				return config.getEdgeFont();
		}});
		treeView.getRenderContext().setEdgeLabelTransformer(new Transformer<Integer, String>(){
			public String transform(Integer s) { // text
				return config.getEdgeLabel(tree.getSource(s), tree.getDest(s));
		 }});
		return treeView;
	}
	
	@SuppressWarnings("serial")
	public static JButton setupButton(final String name, final boolean emph) {
		return new JButton() {
			@Override
		    protected void paintComponent(Graphics g){		
				setFont(LineageDisplayConfig.GUI_FONT);
				setContentAreaFilled(false);
				if(!isEnabled()) {
					setText("Saving...");
				} else {
					setText(name);
				}
				if(emph) {
					setForeground(Color.yellow);
					if(!isEnabled()) {
						setForeground(Color.GRAY);
					}
				} else {
					setForeground(Color.WHITE);
				}
		        Graphics2D g2 = (Graphics2D)g.create();
		        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		        Color c = new Color(87, 82, 126);
		        g2.setPaint(c.darker());
		        g2.fillOval(2, 2, getWidth()-2, getHeight()-2);
		        if (getModel().isRollover()) {
		            g2.setPaint(c.brighter());
		        } else if (getModel().isPressed()) {
		        	g2.setPaint(c.brighter());
		        } else {
		        	g2.setPaint(c);
		        }
		        g2.fillOval(4, 4, getWidth()-8, getHeight()-8);
		        g2.dispose();
		        super.paintComponent(g);
		    }
		};
	}
	
	/**
	 * Launch the interactive lineage tree viewer
	 */
	public static void showLineageTree(final LineageDisplayConfig config) {	
		try {
			UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
		} catch (Exception e1) {
			System.exit(0);
		} 
		
		// main frame
     	final JFrame frame = new JFrame("LICHeE: Lineage Tree Viewer");
		JPanel licheePanel = new JPanel();
		frame.setContentPane(licheePanel);
		frame.setBackground(Color.WHITE);
		licheePanel.setBackground(new Color(78,84,98).brighter());//232,162,104).darker());//(194,110,144));
		
		// tree view
		final JPanel graphPanel = new JPanel();
		graphPanel.setBackground(Color.WHITE);
		graphPanel.setLayout(new GridLayout());
		DelegateTree<Integer, Integer> tree = new DelegateTree<Integer, Integer>(config.displayTree);
		tree.setRoot(0);
		final VisualizationViewer<Integer, Integer> treeView = setupJungViewer(tree, config);
		graphPanel.add(treeView, BorderLayout.PAGE_START);
        
		// buttons
        JButton snapshotButton = setupButton("Snapshot", true);
		JButton removeButton = setupButton("Remove", false);
        JButton collapseButton = setupButton("Collapse", false);
        
        // info console
        final JTextArea consoleTextArea = new JTextArea();
        consoleTextArea.setEditable(false);
        consoleTextArea.setBackground(new Color(238,232,170));
        consoleTextArea.setFont(LineageDisplayConfig.CONSOLE_FONT);
        JScrollPane consolePanel = new JScrollPane(consoleTextArea);
		consolePanel.setViewportBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Node information", 
				TitledBorder.LEADING, TitledBorder.TOP, LineageDisplayConfig.GUI_FONT, Color.yellow));
		
        // setup the layout
        GroupLayout layout = new GroupLayout(licheePanel);
        licheePanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(graphPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(consolePanel, javax.swing.GroupLayout.PREFERRED_SIZE, 500, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                	.addComponent(snapshotButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                	.addGap(47, 47, 47)
                    .addComponent(removeButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGap(47, 47, 47)
                    .addComponent(collapseButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(451, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, true)
            .addGroup(layout.createSequentialGroup()
                .addComponent(graphPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE, Short.MAX_VALUE) //javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    	.addComponent(snapshotButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(removeButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(collapseButton, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(consolePanel))
                .addContainerGap())
        );
        
        // node click listeners
 		final PickedState<Integer> pickedState = treeView.getPickedVertexState();
 		
 		pickedState.addItemListener(new ItemListener() { 
 		    @Override
 		    public void itemStateChanged(ItemEvent e) {
 		    Object o = e.getItem();
 		        if (o instanceof Integer) {
 		            Integer nid = (Integer) o;
 		            if (pickedState.isPicked(nid)) {
 		                consoleTextArea.setText(config.getConsoleText(nid));
 		                consoleTextArea.setCaretPosition(0);
 		                // highlight the path
 		                PHYNode n = config.nodeInfo.get(nid);
 		                if(n.isLeaf()) {
 		                	config.subclone_path_highligth = true;
 		                	ArrayList<ArrayList<PHYNode>> lineages = new ArrayList<ArrayList<PHYNode>>();
 		                	config.spanningTree.getLineageClusters(new ArrayList<PHYNode>(), lineages, config.spanningTree.getRoot(), n.getLeafSampleId());
 		                	for(int i = 0; i < lineages.size(); i++) {
 		                		for(int j = 0; j < lineages.get(i).size(); j++) {
 		                			config.nodeHighlight[lineages.get(i).get(j).getNodeId()] = true;
 		                		}
 		                	}
 		                	config.nodeHighlight[n.getNodeId()] = true;
 		                } else {
 		                	config.subclone_sample_highligth = true;
 		                	config.nodeSelected[nid] = true;
 		                }
 		               treeView.repaint();
 		            } else {
 		            	config.subclone_path_highligth = false;
 		            	config.subclone_sample_highligth = false;
 		            	for(int i = 0; i < config.nodeHighlight.length; i++) {
 		            		config.nodeHighlight[i] = false; // reset
 		            		config.nodeSelected[i] = false;
 		            	}
 		            	treeView.repaint();
 		            }
 		        }
 		    }
 		});
		
		// remove node button
		removeButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent evt) {
				Object[] pickedNodes = new HashSet<Integer>(treeView.getPickedVertexState().getPicked()).toArray();
				if(pickedNodes.length != 1) {
					JOptionPane.showMessageDialog(frame, "One node needs to be selected.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				int selectedNode = (Integer)pickedNodes[0];
				PHYNode n = config.nodeInfo.get(selectedNode);
				if(n.isLeaf()) {
					JOptionPane.showMessageDialog(frame, "Sample nodes cannot be removed.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}		
				if(n.isRoot()) {
					JOptionPane.showMessageDialog(frame, "Cannot remove the root germline node.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				PHYNetwork constrNetwork = config.constrNetwork.removeNode(n);
				ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees(); 
				constrNetwork.evaluateLineageTrees();
				if(spanningTrees.size() > 0) {
					System.out.println("Best tree error score: " + spanningTrees.get(0).getErrorScore());
					LineageDisplayConfig new_config = new LineageDisplayConfig(constrNetwork, spanningTrees.get(0), config.sampleNames, config.beautify);
					Visualizer.showLineageTree(new_config);
				} else {
					JOptionPane.showMessageDialog(frame, "No valid lineage tree was found.", "Message", JOptionPane.PLAIN_MESSAGE);
				}
			}
		});
		
		// collapse button
		collapseButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent evt) {
				Object[] pickedNodes = new HashSet<Integer>(treeView.getPickedVertexState().getPicked()).toArray();
				if(pickedNodes.length != 2) {
					JOptionPane.showMessageDialog(frame, "Two nodes need to be selected.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				PHYNode n1 = config.nodeInfo.get((Integer)pickedNodes[0]);
				PHYNode n2 = config.nodeInfo.get((Integer)pickedNodes[1]);
				if(!n1.getSNVGroup().equals(n2.getSNVGroup())) { // can collapse only clusters
					JOptionPane.showMessageDialog(frame, "Cannot collapse nodes that are not clusters of the same group.", "Message", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				PHYNetwork constrNetwork = config.constrNetwork.collapseClusterNodes(n1, n2);
				ArrayList<PHYTree> spanningTrees = constrNetwork.getLineageTrees();  
				constrNetwork.evaluateLineageTrees();
				if(spanningTrees.size() > 0) {
					LineageDisplayConfig new_config = new LineageDisplayConfig(constrNetwork, spanningTrees.get(0), config.sampleNames, config.beautify);
					Visualizer.showLineageTree(new_config);
				} else {
					JOptionPane.showMessageDialog(frame, "No valid lineage tree was found.", "Message", JOptionPane.PLAIN_MESSAGE);
				}
			}
		});
	
		// snapshot button
		snapshotButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent evt) {
				final JButton button = (JButton) evt.getSource();
				PrintStream err = System.err;
				System.setErr(new PrintStream(new OutputStream() {
				    public void write(int b) {}
				}));
				final PDFGraphics2D bufImage = new PDFGraphics2D(0.0, 0.0, treeView.getSize().width, treeView.getSize().height);
				treeView.paint(bufImage);			
				System.setErr(err); // restore
				
				final JFileChooser fc = new JFileChooser();
				int returnVal = fc.showSaveDialog(graphPanel);
			    if (returnVal == JFileChooser.APPROVE_OPTION) {
			    	button.setEnabled(false);
			    	button.setToolTipText("Saving image...");
			    	Thread worker = new Thread() {
			    		@Override
			    		public void run() {
				            File imageFile = fc.getSelectedFile();
				            try {
				            	imageFile.createNewFile();
				            	byte[] data = bufImage.getBytes();
				            	Files.write(imageFile.toPath(), data);
				            } catch(Exception ex){
				            	System.err.println("Error while saving the image.");
				            }
				            
				            SwingUtilities.invokeLater(new Runnable() {
				                public void run() {
				                  button.setEnabled(true);
				                }
				              });
				        }
			    		
					};	
					worker.start();
			    }
			}
		});
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
}
