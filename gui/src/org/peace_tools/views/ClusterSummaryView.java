//--------------------------------------------------------------------
//
// This file is part of PEACE.
// 
// PEACE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// PEACE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
// 
// Miami University makes no representations or warranties about the
// suitability of the software, either express or implied, including
// but not limited to the implied warranties of merchantability,
// fitness for a particular purpose, or non-infringement.  Miami
// University shall not be liable for any damages suffered by licensee
// as a result of using, result of using, modifying or distributing
// this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   Dhananjai M. Rao              raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.views;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextArea;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.ProgressMonitor;
import javax.swing.SwingUtilities;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.peace_tools.core.MainFrame;
import org.peace_tools.data.ClusterFile;
import org.peace_tools.data.ClusterNode;
import org.peace_tools.data.ESTList;
import org.peace_tools.generic.Pair;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DBClassifier;
import org.peace_tools.workspace.Workspace;
import org.peace_tools.workspace.WorkspaceEvent;
import org.peace_tools.workspace.WorkspaceListener;

public class ClusterSummaryView extends JPanel 
	implements WorkspaceListener, Runnable {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor creates the various controls that
	 * are part of the summary view. In addition, it also sets up
	 * the summary graph with default view and seize.
	 * 
	 * @param frame The main frame that logically owns this view. This
	 * reference is used to launch the classifier dialog.
	 * 
	 * @param clusterFile The data model to be used to display the 
	 * clusters from a given cluster file.
	 * 
	 * @param estList The set of ESTs that contain information about 
	 * each EST in the clusters.
	 */
	public ClusterSummaryView(MainFrame frame, ClusterFile clusterFile, ESTList estList) {
		super(new BorderLayout(0, 0));
		setOpaque(false);
		// Save references to clusterFile for future use
		this.clusterFile = clusterFile;
		this.estList     = estList;
		this.mainFrame   = frame;
		this.toolbar     = new JToolBar();
		// Compute the largest cluster size.
		this.maxClusterSize = clusterFile.getRoot().getLargestClusterSize();
		// Further configure the tool bar: ** do this first **
		setupToolBar();
		// Create and configure the graph pane in the middle
		setupGraph();
		// Finally create the informational panel at the bottom
		setupInfoPanel();
		// Register as work space listener as well.
		Workspace.get().addWorkspaceListener(this);
	}
	
	/**
	 * Helper method to configure and setup various tools in the 
	 * tool bar. This method was introduced to streamline the code
	 * better and cut down the code clutter in the constructor.
	 */
	private void setupToolBar() {
		toolbar.setFloatable(false);
		// Create the x-axis slider. 
		int clusterCount = clusterFile.getRoot().getChildren().size();
		axisScale[0] = createSlider(clusterCount);
		toolbar.add(new JLabel("X Scale:"));
		toolbar.add(axisScale[0]);
		// Create the y-axis scale slider. 
		axisScale[1] = createSlider(maxClusterSize);
		toolbar.add(new JLabel(" Y Scale:"));
		toolbar.add(axisScale[1]);
		// Add spacer to make things look good
		toolbar.add(Box.createHorizontalStrut(10));
		// Create button to set/reset log scale for y-axis
		logScale = new JToggleButton(Utilities.getIcon("images/16x16/LogScale.png"));
		logScale.setToolTipText("Toggle between log & linear scale");
		logScale.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				// Resize/scale the graph.
				ClusterSummaryView.this.resizeGraph();
				ClusterSummaryView.this.repaint();
			}
		});
		toolbar.add(logScale);
		// Create button enable/disable classification information.
		classify = new JToggleButton(Utilities.getIcon("images/16x16/CheckClassifier.png"));
		classify.setToolTipText("Include EST data base classifications in graph");
		classify.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				// Repaint graph with/without classification information.
				ClusterSummaryView.this.checkAndClassify();
			}
		});
		toolbar.add(classify);
		// Add button to launch the DB classifier dialog.
		JButton button = Utilities.createButton("images/16x16/Classifier.png",
				"EST Classifiers", "classifiers", null, 
				"Launch EST classifier editing dialog to set classification rules", true);
		button.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				DBClassifierEditor editor = new DBClassifierEditor(mainFrame);
				editor.pack();
				Utilities.centerPanel(mainFrame, editor);
				editor.setVisible(true);
			}
		});
		toolbar.add(button);
		
		// Finally create the help button with suitable tool tip
		toolbar.add(Box.createHorizontalStrut(10));
		toolbar.add(Utilities.createButton("images/16x16/Help.png", 
				null, "Help", null, "Read about the cluster summary" +
						" view and how to use it", true));
		
		// Add tool bar to the main panel at the top
		add(toolbar, BorderLayout.NORTH);
	}

	/**
	 * Helper method to create a slider with suitable scale.
	 * 
	 * This method is used to create a suitable slider to set size for
	 * both x and y axis scaling. This method ensures the scale of the axis 
	 * cannot be below a value that will make the axis smaller 
	 * than 100 pixels. Note that JSlider takes only integer 
	 * values. So the scales are translated to multiplied by 100.
	 * That is rather than setting scale to 0.08 we set scale to 8.
	 * 
	 * @param size The default/normal size of the axis.
	 * 
	 * @return A JSlider with suitable scale and default value set.
	 */
	private JSlider createSlider(int size) {
		int minValue = 10000 / size; // Ensure minimum size is never < 100
		int maxValue = 500;         // Five times the full size.
		if (maxValue < minValue) {
			maxValue = minValue * 5;
		}
		// The slider to be displayed.
		JSlider slider = new JSlider(JSlider.HORIZONTAL, minValue, maxValue, Math.max(100, minValue));
		// Setup the labels (at fixed points) to be displayed in the slider.
		final int Values[] = {minValue, 100, maxValue}; 
		Hashtable<Integer, JLabel> labelTable = new Hashtable<Integer, JLabel>();
		for(int value: Values) {
			JLabel label = new JLabel("" + value / 100.0);
			// Setup font to be 7 pt to make things look good.
			Utilities.adjustFont(label, -100, 7, 0);
			labelTable.put(value, label);
		}
		// Setup the labels to be displayed in the slider
		slider.setLabelTable(labelTable);
		slider.setPaintLabels(true);
		slider.setOpaque(false);
		// Don't let the slider get too short.
		Dimension sliderSize = slider.getPreferredSize();
		sliderSize.width = 150;
		slider.setMinimumSize(sliderSize);
		// Add a suitable listener to the slider.
		slider.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				// Resize/scale the graph.
				ClusterSummaryView.this.resizeGraph();
				ClusterSummaryView.this.repaint();
			}			
		});
		// Return the slider for further use.
		return slider;
	}
	
	/**
	 * Helper method to configure and setup the informational labels
	 * at the bottom of the graph panel.  The informational panel
	 * provides information about the cluster and how it was generated.
	 * This method was introduced to streamline the code
	 * better and cut down the code clutter in the constructor.
	 */
	private void setupInfoPanel() {
		// Create the information to be displayed to the user
		JTextArea infoArea = new JTextArea(4, 40);
		infoArea.setOpaque(false);
		// Populate info area with meta data from cluster file.
		ArrayList<Pair> metadata = clusterFile.getMetadata();
		// Add each name-value (nv) pair to area info.
		for(Pair nvPair: metadata) {
			infoArea.append(nvPair.toString());
			infoArea.append("\n");
		}
		// Move back to the top of the 
		// Add the name-value pair to a scroll pane in case the 
		// information is large or the window is small.
		JScrollPane jsp = new JScrollPane(infoArea);
		jsp.setOpaque(false);
		jsp.getViewport().setOpaque(false);
		jsp.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), 
				"Clustering Parameters", TitledBorder.RIGHT, TitledBorder.TOP));
		// Add the scroll pane to the main component
		add(jsp, BorderLayout.SOUTH);
	}
	
	/**
	 * Top-level method to draw the cluster graph.
	 * 
	 * This method is the top-level method to draw the cluster graph
	 * associated with this object. This method is invoked from the 
	 * graph JPanel as the paintComponent() call back is intercepted
	 * and forwarded to this method as setup in the setupGraph() 
	 * method. This method draws the graph by suitably delegating
	 * the appropriate tasks to various helper methods.
	 * 
	 * @param g The graphics object to be used to draw the graph.
	 */
	private void drawGraph(Graphics g) {
		Graphics2D g2d = (Graphics2D) g;
		// Rotate and draw the y-axis label.
		AffineTransform oldTrans = g2d.getTransform();
		Rectangle yLabelBounds = yLabel.getBounds();
		g2d.rotate(Math.toRadians(-90), yLabelBounds.x + yLabelBounds.width / 2, yLabelBounds.y + yLabelBounds.height / 2);
		g2d.translate(0, -yLabelBounds.width / 2 + yLabelBounds.height / 2);
		Graphics subPanel = g2d.create(yLabelBounds.x, yLabelBounds.y, yLabelBounds.width, yLabelBounds.height);
		yLabel.paint(subPanel);
		g2d.setTransform(oldTrans);
		// Draw the graph area.
		g2d.setColor(Color.black);
		g2d.draw(area);
		// Now draw the core cluster information at the given
		// zoom level. Only draw the cluster information contained in
		// the clip area to save time.
		ArrayList<ClusterNode> root = clusterFile.getRoot().getChildren(); 
		Rectangle clipArea = g.getClipBounds();
		int startCluster   = toClusterIndex(Math.max(clipArea.x, area.x));
		if (startCluster == -1) {
			// The starting x is beyond the cluster range. So no need to
			// do anything further
			return;
		}
		startCluster       = Math.max(0, startCluster - 1);
		int endCluster     = toClusterIndex(Math.min(clipArea.x + clipArea.width, 
								area.x + area.width));
		endCluster         = Math.min(root.size(), endCluster + 1);
		// Setup colors for all the bars.
		Color barColor     = new Color(0, 0, 0, (axisScale[0].getValue() < 100) ? 255 : 255);
		// Loop over the clusters and draw rectangles.
		final double logMinVal = Math.log10(LogStartValue);
		final double logDenom  = Math.log10(maxClusterSize) - logMinVal;
		
  		for(int clIdx = startCluster; (clIdx < endCluster); clIdx++) {
			ClusterNode cluster = root.get(clIdx);
			int estsInCluster   = cluster.getLargestClusterSize();
			int height          = (estsInCluster * area.height - 1) /  maxClusterSize;
			if (logScale.isSelected()) {
				double tickValue =  (Math.log10(estsInCluster) - logMinVal) / logDenom;
				height = (int) (tickValue * area.height - 1);
			}
			int x1              = getStartX(clIdx); 
			int x2              = getStartX(clIdx + 1);
			Rectangle barArea   = new Rectangle(x1, area.y + area.height - height - 1,
												x2 - x1, height);
			if ((axisScale[0].getValue() > 100) && classify.isSelected()) {
				// Draw EST classification aggregate data.
				drawClassificationBars(g2d, barArea, cluster);
			}
			g2d.setColor(barColor);
			g2d.drawRect(x1, area.y + area.height - height - 1, 
					x2 - x1, height);
		}
	}

	/**
	 * Draw clustering classification information in the graph.
	 * 
	 * This method is invoked from the drawGraph() method to draw the
	 * classification information in the bar using the coloring scheme
	 * setup by the user. This method uses the EST classification 
	 * information and divides the bar area into percentages to 
	 * reflect the groups of ESTs.
	 * 
	 * @param g The graphics object to be used to draw the bars.
	 * 
	 * @param barArea The area of the bar to be filled to reflect the
	 * percentages of cluster classification information.
	 * 
	 * @param cluster The cluster whose classification information is
	 * to be displayed as percentages.
	 */
	private void drawClassificationBars(Graphics g, Rectangle barArea, ClusterNode cluster) {
		final int totalESTs = cluster.getChildren().size();
		int subBarY         = barArea.y + barArea.height; // Tracks unused space.
		final HashMap<Integer, Integer> groups = cluster.getESTClasses();
		// Draw percentage graphs for each group.
		if (groups != null) {
			final ArrayList<DBClassifier> classifiers = 
				Workspace.get().getClassifierList().getClassifiers();
			for(int i = 0; (i < classifiers.size()); i++) {
				Integer count = groups.get(new Integer(i));
				if (count == null) {
					// No entry for this DB classifier.
					continue;
				}
				// Found an entry. Compute percentage and draw bar
				double percentage   = count.doubleValue() / totalESTs;
				int    subBarHeight = (int) (barArea.height * percentage);
				if (subBarHeight < 1) {
					// No point in drawing a 
					continue;
				}
				// Update sub-bar start and end y values.
				subBarY -= subBarHeight;
				// Draw the bar with specified color.
				g.setColor(new Color(classifiers.get(i).getColor()));
				g.fillRect(barArea.x, subBarY, barArea.width, subBarHeight);
			}
		}
		// Do we fill in remaining space for unclassified ESTs with a
		// specific color or pattern?
	}
	
	/**
	 * Helper method to translate a cluster Index to an absolute x-coordinate.
	 * 
	 * This is a helper method that can be used to translate a given cluster
	 * index to an absolute x-coordinate value within the graph pane. 
	 * 
	 * @param clIndex The index of the cluster to be translated.
	 * 
	 * @return The starting x-coordinate of the cluster bar in the cluster
	 * graph.
	 */
	private int getStartX(int clIndex) {
		int maxClusters = clusterFile.getRoot().getChildren().size() + 1;
		if ((clIndex < 0) || (clIndex > maxClusters)) {
			return -1;
		}
		clIndex++;
		int tickPos = (int) (area.width * (double) clIndex / maxClusters);
		int prevTickPos = (int) (area.width * (double) (clIndex - 1) / maxClusters);
		int xPos  = (tickPos + prevTickPos)  / 2;
		xPos     += area.x;
		return xPos;
	}
	
	/**
	 * Helper method to convert an x-coordinate in the graph to a cluster index.
	 * 
	 * @param x The x-coordinate to be translated to a cluster id.
	 * 
	 * @return The index of the cluster that logically contains the given
	 * x-coordinate. If the x-coordinate is not contained by any cluster
	 * then this method returns -1.
	 */
	private int toClusterIndex(int x) {
		if ((x < area.x) || (x > area.x + area.width)) {
			return -1;
		}
		// First obtain the pixel value within bounds.
		int actualX = Math.max(area.x, x);
		actualX     = Math.min(actualX, area.x + area.width);
		actualX    -= area.x;
		// Now translate the pixel to cluster based on scale.
		int clIdx   = actualX * 100 / axisScale[0].getValue();
		clIdx = Math.min(clusterFile.getRoot().getChildren().size(), clIdx);
		return clIdx;
	}
	
	/**
	 * Helper method to configure and setup the graph panel to draw
	 * the actual cluster graph. The graph is displayed in the middle
	 * of the panel. This method was introduced to streamline the code
	 * better and cut down the code clutter in the constructor.
	 */
	private void setupGraph() {
		// Create a JPanel and intercept the paintComponent call to
		// actually draw the graph. The call is forwarded to the
		// ClusterSummaryView.drawGraph() method. If drastically different
		// graphs are to be drawn, then this component should become a
		// class hierarchy of its own.
		graph = new JPanel(null) {
			private static final long serialVersionUID = 5272843175014078369L;
			@Override
			public void paintComponent(Graphics g) {
				super.paintComponent(g);
				ClusterSummaryView.this.drawGraph(g);
			}
		};
		graph.setBackground(Color.white);
		// Let the graph use the initial default as layout
		resizeGraph();
		// Add graph to the main component
		add(new JScrollPane(graph), BorderLayout.CENTER);
	}
	
	/**
	 * Helper method to resize/scale the graph.
	 * 
	 * This method is called initially or whenever the user changes the
	 * scale of the graph. This method uses the x and y axis scale factors
	 * to suitably scale the graph size and relayout the graph.
	 */
	private void resizeGraph() {
		// Create the initial set of labels.
		setupAxisLabels();
		// Set up preferred size for drawing area
		Dimension size = area.getSize();
		size.width    += area.x;
		size.height   += area.y;
		// Account for spacing at the bottom with x-labels
		JLabel dummyLabel = new JLabel("*", JLabel.CENTER);
		Utilities.adjustFont(dummyLabel, 0, 12, 1);
		final int xLabelSpace = (dummyLabel.getPreferredSize().height * 2) + 4;
		size.height += (borderSize * 2);
		size.height += xLabelSpace;
		// Account for borders on either side.
		size.width  += (borderSize * 2);
		// Use the highest value to ensure there is sufficient space for
		// the last x-axis tick.
		final int clusterCount = clusterFile.getRoot().getChildren().size() + 1;
		dummyLabel = new JLabel("" + clusterCount);
		size.width += dummyLabel.getPreferredSize().width / 2;
		// Setup the size for the graph
		graph.setPreferredSize(size);
		graph.setMinimumSize(size);
	}

	/**
	 * Helper method to setup x and y axis labels used to draw
	 * axis information.
	 * 
	 * This helper method is invoked whenever the size of the graph
	 * is changed (by the user via tool bar buttons). This method
	 * computes the labels to be displayed in the graph and then
	 * positions them within the graph pane.
	 * 
	 * <p><b>Note:</b>  The graph (a JPanel) does not have a layout
	 * manager. Therefore the x and y locations of the labels are 
	 * completely under our control.</p>
	 */
	private void setupAxisLabels() {
		// First clear out all the old labels from the graph
		graph.removeAll();
		// Track the space available for graphs
		int clusterCount = clusterFile.getRoot().getChildren().size();
		Dimension graphSize = new Dimension((clusterCount + 2) * axisScale[0].getValue() / 100,
				maxClusterSize * axisScale[1].getValue() / 100);
		area = new Rectangle(borderSize, borderSize, 
				graphSize.width, graphSize.height);
		// Create the y-axis labels
		addYAxisLabels(area);
		// Next create the x-axis labels
		addXAxisLabels(area);
	}
	
	/**
	 * Helper method to create X-axis labels.
	 * 
	 * This is a helper method that is invoked from the setupAxisLabels()
	 * method to setup the labels on the x-axis. This method was
	 * introduced to streamline the code a bit better. 
	 * 
	 * <p><b>Note:</b>  This method must be called only after the 
	 * {@link #addYAxisLabels(Rectangle)} method has been called. This 
	 * is needed to ensure that the area of the graph has been 
	 * suitably updated.</p>
	 * 
	 * @param area The area of the graph to be used for creating the 
	 * labels.
	 */
	private void addXAxisLabels(Rectangle area) {
		final int clusterCount = clusterFile.getRoot().getChildren().size() + 1;
		// First add the main x-axis label
		JLabel axisLabel      = new JLabel("Clusters", JLabel.CENTER);
		Utilities.adjustFont(axisLabel, 0, 12, 1);
		final int labelHeight = axisLabel.getPreferredSize().height;
		axisLabel.setBounds(new Rectangle(area.x, area.y + area.height + 2 + labelHeight, 
				area.width, labelHeight));
		graph.add(axisLabel);
		// Now create all the x-axis ticks/labels
		ArrayList<JLabel> xLabels = new ArrayList<JLabel>(10); 
		createLabels(0, clusterCount, 50, area.width, true, area, xLabels);
		// Finally add y ticks
		createTicks(xLabels, true);
	}

	/**
	 * Helper method to create Y-axis labels.
	 * 
	 * This is a helper method that is invoked from the setupAxisLabels()
	 * method to setup the labels on the y-axis. This method was
	 * introduced to streamline the code a bit better. 
	 * 
	 * <p><b>Note:</b>  This method must be called before the 
	 * {@link #addXAxisLabels(Rectangle)} method. This method sets up
	 * the graph area appropriately for the 
	 * {@link #addXAxisLabels(Rectangle)} method. </p>
	 * 
	 * @param area The area of the graph to be used for creating the 
	 * labels.
	 */
	private void addYAxisLabels(Rectangle area) {
		// Add y-axis main label. 
		yLabel = new JLabel("EST per Cluster", JLabel.CENTER);
		Utilities.adjustFont(yLabel, 0, 12, 1);
		if (yLabel.getPreferredSize().getWidth() >= area.height) {
			// The long label will not fit nicely. Use a shorter version
			yLabel = new JLabel("#EST", JLabel.CENTER);
			Utilities.adjustFont(yLabel, 0, 12, 1);
		}
		// Y label will be drawn rotated. So use height as width
		int labelHeight = yLabel.getPreferredSize().height;
		yLabel.setBounds(new Rectangle(area.x, (area.y + area.height / 2) - 
				(labelHeight / 2), yLabel.getPreferredSize().width, labelHeight));
		// Update area start having added y-axis label
		area.x     += labelHeight;
		// area.width -= labelHeight;
		// Add some padding space at the top to ensure that the
		// top-most label does not spill into the borders
		JLabel dummy = new JLabel("*");
		area.y      += dummy.getPreferredSize().height / 2;
		// Create the y-axis ticks.
		ArrayList<JLabel> yLabels = new ArrayList<JLabel>(10);
		final float minValue = (logScale.isSelected() ? LogStartValue : 0);
		Dimension maxSize = createLabels(minValue, maxClusterSize, 50, 
				area.height, false, area, yLabels);
		// Y-axis labels look best when right justified rather then
		// being left-justified. So do this.
		for(JLabel label: yLabels) {
			Rectangle bounds = label.getBounds();
			bounds.width     = maxSize.width;
			label.setBounds(bounds);
		}
		// Adjust the starting point of drawable graph area
		area.x += (maxSize.width + 2);
		// area.width -= (maxSize.width + 2);
		// Finally add y ticks
		createTicks(yLabels, false);
	}

	/**
	 * Helper method to create little tick lines on x or y axis.
	 * 
	 * This is a helper method that is used create small tick lines
	 * on the x-axis or the y-axis.
	 * 
	 * @param labelList The labels to be used to determine the location
	 * of the ticks.
	 * 
	 * @param xAxis If this parameter is true, then the ticks are
	 * for the x-axis. Otherwise the ticks are to be created for the
	 * y-axis.
	 */
	private void createTicks(ArrayList<JLabel> labelList, boolean xAxis) {
		// The axis size to be used for scaling.
		final int AxisSize = (xAxis ? area.width : area.height);
		// Use the tick labels to determine start and end values
		float startValue = Float.parseFloat(labelList.get(0).getText());
		float endValue   = Float.parseFloat(labelList.get(labelList.size() - 1).getText());
		// Create the ticks using values in the labels.
		for(JLabel label : labelList) {
			float value = Float.parseFloat(label.getText());
			int tickPos = (int) (AxisSize * (double) (value - startValue) / (endValue - startValue));
			tickPos    += (xAxis ? area.x : 0);
			// Set up
			if (!xAxis && logScale.isSelected()) {
				double tickValue =  (Math.log10(value) - Math.log10(startValue)) /
					(Math.log10(endValue) - Math.log10(startValue));
				tickPos = (int) (tickValue * AxisSize);
			}
			// Create tick line and setup bounds.
			TickLine tick = new TickLine();
			if (xAxis) {
				// 1 pixel thick vertical line
				tick.setBounds(tickPos, area.y + area.height - 2, 1, 4); 
			} else {
				// Make a 1 pixel horizontal line
				tick.setBounds(area.x - 1, area.y + area.height - tickPos, 4, 1);
			}
			graph.add(tick);
		}
	}
	
	/**
	 * Helper method to create the axis labels.
	 * 
	 * This helper method is used to add labels on both the x and the
	 * y axes. This method performs some initial math to ensure that
	 * the labels are reasonably well spaced out. 
	 * 
	 * @param startValue The starting value for the axis.
	 * @param endValue The ending value for the axis.
	 * @param tickPixelSpacing The pixels between ticks.
	 * @param axisLength The net length of the axis.
	 * @param xAxis If this is true, then the labels are for the x-axis.
	 * Otherwise labels are to be created for the y-axis.
	 * @param area The net area of the graph.
	 * @param labelList An array list to which the labels created are added.
	 * This list can be null (then labels are not added to the list).
	 * @return The dimension of the largest (or longest) label created.
	 * This value is typically used to right-justify labels on the y-axis.
	 */
	private Dimension createLabels(float startValue, int endValue, final int tickPixelSpacing, 
			int axisLength, boolean xAxis, Rectangle area, 
			ArrayList<JLabel> labelList) {
		// First compute inter-label spacing based on the spacing value
		int tickCount = (int) Math.ceil((double) axisLength / tickPixelSpacing);
		// Ensure ticks are not going to be too crowded using 80% of preferred
		// size as the threshold
		if ((tickPixelSpacing * 0.8) > (axisLength / tickCount)) {
			// Cut one tick out to ensure ticks are going to fit.
			tickCount--;
		}
		// Compute actual values to skip per tick (approximately)
		final int tickValueSkip = (int) Math.ceil((double) (endValue - startValue) / tickCount);
		// Add 1 to tick count to account for last tick.
		tickCount++;
		// The following variable track max width of labels which is
		// used primarily to update y-axis labels.
		Dimension maxSize = new Dimension(0, 0);
		float prevValue   = -1;
		for(int tick = 0; (tick < tickCount); tick++) {
			float labelValue = 0;
			// Compute label value differently for y-axis in log scale.
			if (!xAxis && logScale.isSelected() && (tick > 0)) {
				// For log scale update label values into log scale.
				double tickPos= tick * tickValueSkip;
				tickPos       = (tickPos - startValue) / (endValue - startValue);
				double denom = Math.log10(endValue) - Math.log10(startValue);
				double tickValue = (tickPos * denom) + Math.log10(startValue);
				labelValue = (float) Math.pow(10, tickValue);
				// Truncate to the first digit.
				labelValue = (int) labelValue;
			} else {
				// Linear axis or the first tick.
				labelValue = startValue + (tick * tickValueSkip);
			}
			
			labelValue       = Math.min(labelValue, endValue); 
			if (tick == tickCount - 1) {
				// For last tick force the value to be the last correct one
				labelValue = endValue;
			}
			if (labelValue == prevValue) {
				// Labels can be too close if scale is big. So skip
				// redundant labels.
				continue;
			}
			prevValue        = labelValue;
			int tickPos      = (int) (axisLength * (double) (labelValue - startValue) / (endValue - startValue));
			if (!xAxis && this.logScale.isSelected()) {
				double tickValue =  (Math.log10(labelValue) - Math.log10(startValue)) /
		        	(Math.log10(endValue) - Math.log10(startValue));
				tickPos      = (int) (tickValue * axisLength);
			}
			tickPos         += (xAxis ? area.x : 0);
			JLabel tickLabel = new JLabel();
			// Round label value to integer if possible
			if (labelValue - (int) labelValue < 1e-4f) {
				tickLabel.setText("" + (int) labelValue);
			} else {
				tickLabel.setText("" + labelValue);
			}
			Utilities.adjustFont(tickLabel, -100, TickFontSize, 0);
			// Set label alignment based on axis
			tickLabel.setHorizontalAlignment(xAxis ? JLabel.CENTER : JLabel.RIGHT);
			tickLabel.setVerticalAlignment(xAxis ? JLabel.CENTER : JLabel.TOP);
			// Use the labels preferred size to setup its bounds.
			Dimension size   = tickLabel.getPreferredSize();
			if (xAxis) {
				// Setup bounds for a x-axis label such that the label is 
				// horizontally centered around the tick.
				tickLabel.setBounds(tickPos - size.width / 2, 
						area.y + area.height + 1, size.width, size.height);
			} else {
				// Setup bounds for a y-axis label such that the label is
				// vertically centered around the tick.
				tickLabel.setBounds(area.x, area.y + area.height - tickPos - size.height / 2, 
						size.width, size.height);
			}
			// Track tallest and widest label so far.
			maxSize.width  = Math.max(maxSize.width,  size.width);
			maxSize.height = Math.max(maxSize.height, size.height);
			// Add the tick to the tick list and the component
			graph.add(tickLabel);
			if (labelList != null) {
				labelList.add(tickLabel);
			}
		}
		// Return the dimension of the largest label
		return maxSize;
	}
	
	/**
	 * Helper method to check and launch classification from another thread.
	 * 
	 * This method is a utility method that is invoked from a couple of
	 * different listener call back methods to check and classifiy the
	 * data as necessary.
	 */
	private void checkAndClassify() {
		if (!classify.isSelected() || 
			(estList.isClassified() && clusterFile.isClassified())) {
			// Classifiers are not enabled. Simply post repaint
			// event and we will be done.
			repaint();
			return;
		}
		// When control drops here we need to reclassify the data
		// Classification can be a long running process. So perform
		// reclassification from a separate thread to ensure that the
		// GUI does not appear to have hung.
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				Thread t = new Thread(ClusterSummaryView.this);
				t.start();
			}
		});
	}
	
	/**
	 * Perform classification from a separate thread.
	 * 
	 * This is a helper method that is invoked from a separate thread
	 * to perform classification of EST data and subsequently cluster
	 * information.
	 */
	@Override
	public void run() {
		JLabel msg = new JLabel("Please wait while classification data is computed...",
				Utilities.getIcon("images/16x16/HourGlass.png"), JLabel.CENTER);
		ProgressMonitor pm = new ProgressMonitor(this, msg, 
				"Classifying ESTs...", 0, estList.getESTs().size() + 
				clusterFile.getRoot().getChildren().size());
		// First get the ESTs classified.
		estList.classify(pm);
		pm.setProgress(estList.getESTs().size());
		// Next get the clusters to classify their information.
		pm.setNote("Collating data for clusters...");
		clusterFile.classify(estList, pm);
		pm.setProgress(pm.getMaximum());
		// Now that the classification is done, post a repaint to
		// ensure the display is updated correctly.
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				ClusterSummaryView.this.repaint();
			}
		});
	}
	
	/**
	 * Intercept and update graph if classifiers change.
	 * 
	 * This method intercepts notifications from the work space and
	 * if the classifier list in the work space changes and the 
	 * classification option is enabled in this graph, then this method
	 * reclassifies the ESTs and redisplays the graph with classifier
	 * information.
	 */
	@Override
	public void workspaceChanged(WorkspaceEvent event) {
		if (WorkspaceEvent.EntryType.CLASSIFIER_LIST.equals(event.getEntryType())) {
			checkAndClassify();
		}
	}
	
	/**
	 * A simple component to draw vertical or horizontal ticks.
	 * 
	 * This is a simple component that merely fills its area with 
	 * black color. This component is used to draw vertical and
	 * horizontal ticks in the graph.
	 */
	private class TickLine extends JComponent {
		/**
		 * A generated serialization UID to keep compiler happy.
		 */
		private static final long serialVersionUID = 3622568645120048270L;

		/**
		 * This method overrides the default paint method to fill the
		 * area with black. This makes the ticks appear black in color.
		 * 
		 * @param g The graphics object to be used to fill the tick
		 * area.
		 */
		@Override
		public void paintComponent(Graphics g) {
			g.setColor(Color.black);
			g.fillRect(0, 0, getWidth(), getHeight());
		}
	}
	
	/**
	 * This is a constant that provides the maximum font size we
	 * would like to use for the ticks on the graph. 
	 */
	private final int TickFontSize = 10;
	
	/**
	 * The x and y axis scale controllers.
	 * 
	 * The following array holds two sliders that can be used to 
	 * control the scale of the x-axis and y-axis in the graph. The
	 * sliders are created in the setupToolBar() method. The sliders
	 * are set such that the scale of the graph cannot be set to a 
	 * value which reduces the dimension of the graph below 100x100.
	 */
	private JSlider axisScale[] = new JSlider[2];
	
	/**
	 * The main y-axis label is drawn separately so that it can be
	 * rotated and drawn vertically in the paintGraph() method.
	 * Consequently, it is maintained as a separate instance 
	 * variable (unlike other labels).
	 */
	private JLabel yLabel;
	
	/**
	 * The area in the graph pane where the actual graph is drawn. This
	 * area is set by the various method to ensure that the graph drawn
	 * does not interfere with the axis and borders etc.
	 */
	private Rectangle area;
	
	/**
	 * Empty white space around the core graph to make things look
	 * good.
	 */
	private int borderSize = 10;
	
	/**
	 * The panel that is used to draw the cluster graph. The panel
	 * is created in the setupGraph() method. The paintComponent
	 * call-back is intercepted by this component to draw the 
	 * graph.  
	 */
	private JPanel graph;
	
	/**
	 * The tool bar that contains some commonly used tools with 
	 * the summary view and for modifying the graph size / info.
	 */
	private final JToolBar toolbar;
	
	/**
	 * A toggle button which when checked (or depressed) indicates
	 * that the y-axis must be drawn on a logarithmic scale rather
	 * than on a linear scale.
	 */
	private JToggleButton logScale;
	
	/**
	 * A toggle button which when checked (or depressed) indicates
	 * that the bars must also include the classification information
	 * for each cluster displayed in the graph.
	 */
	private JToggleButton classify;
	
	/**
	 * The largest cluster in the set of clusters that we are working with.
	 * This information is set once and reused to draw the graphs.
	 */
	private final int maxClusterSize;
	
	/**
	 * The set of clusters that we are working with. The cluster
	 * data is set when this class is instantiated.
	 */
	private final ClusterFile clusterFile;

	/**
	 * The set of ESTs that the clusters correspond to. This data is
	 * actually used only for classification.
	 */
	private final ESTList estList;
	
	/**
	 * A convenience constant to set the smallest y-axis value when the
	 * data is displayed on a log scale.
	 */
	private final float LogStartValue = 0.5f;
	
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this component. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * Generated serialization UID (merely to keep the compiler happy)
	 */
	private static final long serialVersionUID = 4287833783594312723L;
}
