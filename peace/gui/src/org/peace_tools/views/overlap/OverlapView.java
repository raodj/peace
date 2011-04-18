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

package org.peace_tools.views.overlap;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.peace_tools.core.MainFrame;
import org.peace_tools.data.ClusterFile;
import org.peace_tools.data.OverlapModel;
import org.peace_tools.generic.HelpHandler;
import org.peace_tools.generic.Utilities;
import org.peace_tools.views.PropertiesTreeMaker;

/**
 * A pre-alignment overlap view of fragments organized based on clusters.
 * 
 * This is a graphical summary view of various fragments organized based
 * on the clusters in which they have been placed. This view uses the
 * pre-alignment model to obtain the information to be displayed in a 
 * convenient form. In addition, this view embellishes the data with the
 * following information:
 * 
 * <ul>
 * 
 * <li>Mechanism to sort clusters based on their sizes.</li>
 * 
 * <li>Provides a mechanism to color individual clusters with different
 * color.</li>
 * 
 *  <li>Ability to scale the set of fragments to obtain detailed or
 *  summary views of overlap view.</li>
 * 
 * </ul>
 */
public class OverlapView extends JPanel implements ChangeListener, 
	ClusterColorMapper, ListSelectionListener {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor creates the various controls (that
	 * are part of this view) that can be used to customize this view. 
	 * 
	 * @param frame The main frame that logically owns this view. This
	 * reference is used to launch the classifier dialog.
	 * 
	 * @param pam The pre-alignment in-memory model that contains the clusters
	 * and organized ESTs. This model is used by this view to render
	 * the information appropriately.
	 * 
	 * @param clusterFile The cluster file from where the clusters are being
	 * displayed. This data is used to merely populate the summary view in 
	 * this class.
	 */
	public OverlapView(MainFrame frame, OverlapModel pam, ClusterFile clusterFile) {
		super(new BorderLayout(0, 0));
		setOpaque(false);
		// Save references to clusterFile for future use
		this.pam        = pam;
		this.mainFrame  = frame;
		this.findHelper = new FindHelper(pam);
		// Create the panel to display fragments
		panel = new OverlapPanel(this.pam, this, findHelper, 10, 1.0 / MIN_COL_SCALE_FACTOR);
		// Create the cluster list using the memory model.
		clusterList = new ClusterList(this.pam, this, this);
		// Create and configure the tool bar
		this.toolbar   = new JToolBar();	
		setupToolBar();
		
		// Wrap the view panel in a scroll panel.
		JScrollPane jsp = new JScrollPane(panel);
		jsp.setViewportBorder(BorderFactory.createEmptyBorder(1, 1, 1, 1));
		// Add a horizontal ruler to display column numbers
		horizontalRuler = new Ruler(pam.getMaxCol(), pam.getMaxRow(), 
				Ruler.Orientation.HORIZONTAL);
		horizontalRuler.setScale(panel.getRowScale(), panel.getColScale());
		jsp.setColumnHeaderView(horizontalRuler);
		// Add a vertical ruler to display row numbers
		verticalRuler = new Ruler(pam.getMaxCol(), pam.getMaxRow(), 
				Ruler.Orientation.VERTICAL);
		verticalRuler.setScale(panel.getRowScale(), panel.getColScale());
		jsp.setRowHeaderView(verticalRuler);
	
		// Create and setup the highlight panel.
		highlightPanel = new HighlightPanel(pam.getESTList());
		highlightPanel.getTable().getSelectionModel().addListSelectionListener(this);
		// Set up the table model for overlap panel's highlighting feature
		panel.setESTTableModel(highlightPanel.getModel());
		
		JSplitPane mainPane =  
			PropertiesTreeMaker.createPropertiesLayout("Highlighting Panel", highlightPanel, 
					jsp, toolbar, toolbar.getComponentCount() - 2, "Highlight", 
					"images/16x16/Highlight.png", 
					"(Un)Highlight selected ESTs in the overlap view");
		
		// Next create the summary information tree...
		JTree summaryInfo = PropertiesTreeMaker.makeClusterProperties(pam.getWsEntry(), mainFrame);
		// Create a split pane with the overlap view and summary 
		JSplitPane contentPane =  
			PropertiesTreeMaker.createPropertiesLayout("Clustering Information", summaryInfo, 
					mainPane, toolbar, toolbar.getComponentCount() - 2);
		// Add the scroll pane to the main panel.
		add(contentPane, BorderLayout.CENTER);
		// Initialize colors for all clusters to grey
		clusterList.setColors(null);
	}
	
	/**
	 * Helper method to configure and setup various tools in the 
	 * tool bar. This method was introduced to streamline the code
	 * better and cut down the code clutter in the constructor.
	 */
	private void setupToolBar() {
		toolbar.setFloatable(false);
		toolbar.add(Box.createHorizontalStrut(5));
		toolbar.add(new JLabel("Clusters: ", Utilities.getIcon("images/16x16/Cluster.png"), JLabel.LEFT));
		toolbar.add(clusterList);
		toolbar.add(Box.createHorizontalStrut(5));
		// Create the x-axis slider. Note that there is a bit of kludge on
		// the range to permit fractional values on the slider (a slider
		// by default only permits integral values. However, here if the
		// value is less than 4, then it is interpreted as a fractional scale.
		double colScale = panel.getColScale();
		colScale = (colScale < 1) ? (colScale * MIN_COL_SCALE_FACTOR) : 
			(colScale + MIN_COL_SCALE_FACTOR);
		colScale = Math.min(1.0, colScale);
		axisScale[0] = new JSlider(JSlider.HORIZONTAL, 1, 20, (int) colScale);
		axisScale[0].setOpaque(false);
		toolbar.add(new JLabel("Column:"));
		toolbar.add(axisScale[0]);
		// Create the y-axis scale slider. 
		axisScale[1] = new JSlider(JSlider.HORIZONTAL, 1, 18, (int) panel.getRowScale());
		axisScale[1].setOpaque(false);
		toolbar.add(new JLabel(" Height:"));
		toolbar.add(axisScale[1]);
		// Setup suitable listeners to the sliders
		axisScale[0].addChangeListener(this);
		axisScale[1].addChangeListener(this);
		// ------------------------------------------------------------------------
				
		// Finally create the help button with suitable tool tip
		toolbar.add(Box.createHorizontalStrut(10));
		JButton button = Utilities.createToolButton("images/16x16/Help.png", 
						null, "Help", null, "Read about the overlap (pre-alignment) summary" +
						" view and how to use it", true);
		button.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				HelpHandler.showHelp(mainFrame, "http://www.peace-tools.org/downloads/manual.pdf#page=35");
			}
		});
		toolbar.add(button);		
		toolbar.add(Box.createHorizontalStrut(5));
		// Add tool bar to the main panel at the top
		add(toolbar, BorderLayout.NORTH);
	}
	
	@Override
	public Color getColor(Integer clusterID) {
		return clusterColorMap.get(clusterID);
	}

	@Override
	public void setColor(Integer clusterID, Color color) {
		clusterColorMap.put(clusterID, color);
	}
	
	/**
	 * Helper method to compute and update scale for the actual display.
	 * 
	 * This method is is part of the ChangeListener interface. This method
	 * is invoked whenever the user modifies the sliders for "Width"
	 * and "Height" in this panel. This method computes a new scaling
	 * factor for rows and columns and updates the scaling factor in the
	 * actual fragment display panel.
	 * 
	 * @param e The change event associated with this method. This
	 * method is currently not used.
	 * 
	 */
	public void stateChanged(ChangeEvent e) {
		// Compute and set appropriate scale for the fragments rendering
		// panel in this view class.
		double colScale = axisScale[0].getValue();
		// Normalize column scale
		colScale = (colScale > MIN_COL_SCALE_FACTOR) ? 
				(colScale - MIN_COL_SCALE_FACTOR) : (colScale / MIN_COL_SCALE_FACTOR);
		// Setup the scale in fragment display panel
		panel.setScale(axisScale[1].getValue(), colScale);
		// Setup the scale in the vertical and horizontal rulers as well.
		horizontalRuler.setScale(axisScale[1].getValue(), colScale);
		verticalRuler.setScale(axisScale[1].getValue(), colScale);
		// Repaint to reflect any changes in bounds
		repaint();
	}
	
	/**
	 * This method is called whenever the user selects a different row in the
	 * highlight table. This method uses the EST index from the currently selected
	 * row in the {@link HighlightPanel} and sets the currently selected entry to
	 * that value.
	 * 
	 * @param lse The list event associated with this call. Currently this method
	 * is ignored and the currently selected row is directly obtained from the table.
	 */
	@Override
	public void valueChanged(ListSelectionEvent lse) {
		JTable estTable = highlightPanel.getTable();
		if (estTable.getSelectedRow() == -1) {
			return;
		}
		// First obtain EST at the specified index
		int estIndex = (Integer) estTable.getValueAt(estTable.getSelectedRow(), 1);
		// Have the corresponding fragment found.
		if (findHelper.find(estIndex) != null) {
			// Have the overlap panel scroll to this est
			this.panel.scrollToCurrent();
			// Ensure that the panel gets rfreshed
			this.panel.repaint();
		}
	}
	
	/**
	 * The row and column size scale controllers.
	 * 
	 * The following array holds two sliders that can be used to 
	 * control the scale of the row and column in the graph. The
	 * sliders are created in the setupToolBar() method. The sliders
	 * are set such that the scale for the row cannot be decreased below
	 * 1 pixel and the column size cannot drop below 1/4 pixel (four
	 * nucleotides per pixel).
	 */
	private JSlider axisScale[] = new JSlider[2];	
		
	/**
	 * The tool bar that contains some commonly used tools with 
	 * the overlap view and for modifying the scale used to
	 * draw the overlap view.
	 */
	private final JToolBar toolbar;
	
	/**
	 * The list that displays the list of clusters currently present in
	 * the memory model and permits the user to change colors. This object
	 * is created in the constructor and is added to the tool bar as the
	 * first component.
	 */
	private final ClusterList clusterList;
	
	/**
	 * The PreAlignmentModel that provides the necessary information
	 * for rendering the core information required for this view.
	 * This instance variable is initialized by the constructor and
	 * it is never changed during the life time of this object.
	 */
	private final OverlapModel pam;

	/**
	 * The helper class that tracks the currently active EST fragment
	 * and assists in search operations. This class is used to determine
	 * the currently active entry so that it can be suitably highlighted
	 * in the overlap view.
	 */
	private final FindHelper findHelper;
	
	/**
	 * The overlap panel that actual renders the fragments for us.
	 * This instance variable holds an immutable reference to the
	 * overlap panel that actually renders the fragments.
	 */
	private final OverlapPanel panel;
	
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this component. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * The horizontal ruler that is displayed at the top of the overview
	 * panel. This ruler is added to the scroll pane that contains the
	 * {@link #panel}. The horizontal ruler displays a ruler with 
	 * column numbers indicating the relative position of each nucleotide.
	 */
	private final Ruler horizontalRuler;

	/**
	 * The vertical ruler that is displayed to the left of the overview
	 * panel. This ruler is added to the scroll pane that contains the
	 * {@link #panel}. The vertical ruler displays a ruler with 
	 * row numbers indicating the relative row number where a given
	 * fragment was accomodated.
	 */
	private final Ruler verticalRuler;
	
	/**
	 * Factor to determine minimum column scale value.
	 * 
	 * This constant provides a convenient mechanism to quickly change
	 * the minimum scaling factor for columns in the fragment display
	 * panel. The minimum scaling factor is computed as:
	 * 1.0 / MIN_COL_SCALE_FACTOR. This value is used in various
	 * spots in the code to configure and work with sliders that requires
	 * integral values to operate (it cannot deal with fractional values
	 * like 0.25)
	 */
	private final double MIN_COL_SCALE_FACTOR = 10.0;
	
	/**
	 * A hash map to quickly look-up color for a given cluster ID.
	 * 
	 * This hash map is used when rendering fragments. The hash map is
	 * used to look up the colors set for a given cluster, given its
	 * cluster ID. This hash map is used to implement methods in the
	 * ClusterColorMapper interface.
	 */
	private final HashMap<Integer, Color> clusterColorMap =
		new HashMap<Integer, Color>();
	
	/**
	 * The panel (that is displayed to the left of this overlap view) that
	 * permits the user to conveniently select the set of ESTs to be 
	 * highlighted in the overlap view.
	 */
	private final HighlightPanel highlightPanel;
	
	/**
	 * Generated serialization UID (merely to keep the compiler happy)
	 */
	private static final long serialVersionUID = 4287833783594312723L;
}
