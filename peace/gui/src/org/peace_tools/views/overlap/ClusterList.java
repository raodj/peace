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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.peace_tools.data.ClusterNode;
import org.peace_tools.data.OverlapModel;
import org.peace_tools.generic.Utilities;

/**
 * A combo-box like custom component to display list of clusters from
 * a tool bar.
 *
 * <p>This component provides a custom GUI component that looks similar
 * to a combo-box but uses a JList in its popup to display a set of 
 * choices that permit multiple-selections. This permits the user to
 * select multiple clusters and apply colors (and possibly other
 * transformations) to the selected set. This component consists of
 * the following set of independent sub-components:</p>
 * 
 * <ul>
 * <li>A text-box which displays the current number of selected 
 * components.</li>
 * <li>A small pull-down button that is used to display the full
 * list of clusters for operations.</li>
 * <li>Several tool bar buttons to change colors of selected cluster
 * entries.</li> 
 * </ul>
 * 
 * <p><b>Note</b>This class is meant to be instantiated and used
 * only by the OverlapView class. Consequently it has been made 
 * package private.</p>
 */
class ClusterList extends JPanel implements ActionListener {
	/**
	 * The default constructor. 
	 * 
	 * @param overlapModel The overlap model from where the data regarding
	 * the set of clusters is retrieved.
	 * 
	 * @param colorMapper The helper object to be used to determine the
	 * color for rendering fragments belonging to various clusters.
	 * 
	 * @param owner The owner to be repainted whenever colors for clusters
	 * are changed.
	 */
	public ClusterList(OverlapModel overlapModel, ClusterColorMapper colorMapper, 
			JComponent owner) {
		// Save references to clusterFile for future use
		this.overlapModel = overlapModel;
		this.colorMapper  = colorMapper;
		this.owner        = owner;
		// The border color for static controls.
		borderColor = getBackground().darker();
		// Setup a default background color.
		setBackground(Color.white);
		// Create adjacent drop down button.
		dropDownButton = Utilities.createButton("images/16x16/MoreDetails.png", 
				null, "popup", this, "Show/hide cluster list", true);
		dropDownButton.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createLineBorder(borderColor), 
				BorderFactory.createEmptyBorder(2, 0, 2, 0)));
		dropDownButton.setDefaultCapable(false);
		dropDownButton.setFocusable(false);
		// Create a label to show selected cluster count.
		selectionStats = new JLabel("0 Selected");
		Dimension prefSize = selectionStats.getPreferredSize();
		prefSize.width = prefSize.width * 3 / 2;
		selectionStats.setPreferredSize(prefSize);
		// Create a tool bar button to set colors for selected or all clusters
		changeColor = Utilities.createButton("images/16x16/Colors.png", 
				null, "changeColor", this, "Change color for selected or all clusters", true);
		changeColor.setBorder(BorderFactory.createLineBorder(borderColor));
		Utilities.makeToolBarButton(changeColor, true);
		changeColor.setDefaultCapable(false);
		// Create button for randomly assigning colors.
		randomColors = Utilities.createToolButton("images/16x16/ColorsGrey.png", 
				null, "randomColor", this, "Randomly assign colors for selected or all clusters", true);
		randomColors.setBorder(BorderFactory.createLineBorder(borderColor));
		Utilities.makeToolBarButton(randomColors, true);
		randomColors.setDefaultCapable(false);
		// Setup a Box layout to layout all the components horizontally
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
		add(dropDownButton); add(Box.createHorizontalStrut(3));
		add(selectionStats); add(Box.createHorizontalStrut(3));
		add(Box.createHorizontalGlue());
		add(changeColor);    add(Box.createHorizontalStrut(3));
		add(randomColors);   add(Box.createHorizontalStrut(3));
		// Setup a nice etched border around the whole pane.
		// setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEtchedBorder(), 
		//		BorderFactory.createEmptyBorder(1, 1, 1, 1)));
		// Finally create the popup cluster list.
		createClusterListAndPopup();
	}
	
	@Override
	public void actionPerformed(ActionEvent event) {
		final String cmd = event.getActionCommand();
		if ("popup".equals(cmd)) {
			// Show or hide the popup menu depending on its visibility
			if (!popupList.isVisible()) {
				// Set suitable size for the pop-up menu depending on 
				// current dimension of this component.
				Dimension prefSize = clusterList.getPreferredSize();
				prefSize.width     = this.getWidth();
				if (clusterList.getModel().getSize() > 8) {
					// Account for vertical scroll bar in width to make layout nicer
					prefSize.width -= clusterListJSP.getVerticalScrollBar().getPreferredSize().getWidth();
				}
				clusterList.setPreferredSize(prefSize);
				// Show the pop-up at the appropriate relative location
				popupList.show(this, 0, this.getHeight());
			} else {
				popupList.setVisible(false);
			}
		} else if ("changeColor".equals(cmd)) {
			// Let user choose a new color.
			Color color = JColorChooser.showDialog(this, "Select cluster color", null);
			if (color != null) {
				// Setup color using utility method.
				ArrayList<Color> colorList = new ArrayList<Color>(1);
				colorList.add(color);
				setColors(colorList);
				// Get the owner to repaint its children
				owner.repaint();
			}
		} else if ("randomColor".equals(cmd)) {
			// Assign random colors to selected clusters
			Color stdColors[] = {Color.blue, Color.cyan, Color.darkGray, Color.gray, Color.green, 
					Color.lightGray, Color.magenta, Color.orange, Color.pink, Color.red, 
					Color.yellow, new Color(0xc0, 0xc0, 0xff)};
			// Use helper method to assign colors.
			setColors(Arrays.asList(stdColors));
			// Get the owner to repaint its children
			owner.repaint();
			// Check to display warning about potential confusion in colors if needed.
			if (stdColors.length < overlapModel.getSize()) {
				JOptionPane.showMessageDialog(this, LOW_COLOR_MSG, 
						"Insufficient color resolution", JOptionPane.WARNING_MESSAGE);
			}
		}
	}
		
	@Override
	public void paintComponent(Graphics g) {
		// Get the bounds using the reference button size.
		Rectangle bounds = this.dropDownButton.getBounds();
		// Paint the current background color first.
		g.setColor(getBackground());
		g.fillRect(bounds.x, bounds.y, getWidth() - 1, bounds.height - 1);
		// Paint an internal type border to make things look good.
		g.setColor(borderColor);
		g.drawRect(bounds.x, bounds.y, getWidth() - 1, bounds.height - 1);
	}
	
	/**
	 * Sets color for selected or all clusters to given color.
	 * 
	 * This method is a common helper method that is used to change 
	 * (or set) colors for either all clusters or just the selected
	 * clusters (as the case may be) to a given set of colors. The 
	 * set of colors to be used are passed in as the parameter. This
	 * method reuses colors in the list if the number of selected
	 * clusters are greater than the number of colors.
	 * 
	 * @param colorList The list of colors to be used for assigning
	 * colors to various clusters. If this list is null, then a default
	 * color of grey is assigned to each cluster.
	 */
	protected void setColors(List<Color> colorList) {
		// Setup a default list to make logic easier below.
		if (colorList == null) {
			colorList = new ArrayList<Color>();
		}
		// Setup a default color to make logic easier below.
		if (colorList.isEmpty()) {
			colorList.add(new Color(0xe0, 0xe0, 0xe0));
		}
		// Set colors for all clusters or just selected clusters
		// as the case may be.
		int[] selectedList = clusterList.getSelectedIndices();
		if ((selectedList == null) || (selectedList.length == 0)) {
			// Set colors for all entries. For this take a short cut
			// by temporarily selecting all entries.
			clusterList.setSelectionInterval(0, overlapModel.getSize() - 1);
			selectedList = clusterList.getSelectedIndices();
			clusterList.clearSelection();
		}
		assert ( selectedList != null );
		int colorIdx = 0;
		for(int i = 0; (i < selectedList.length); i++) {
			ClusterNode clsNode = overlapModel.getElementAt(selectedList[i]);
			colorMapper.setColor(clsNode.getClusterId(), colorList.get(colorIdx));
			colorIdx = (colorIdx + 1) % colorList.size();
		}
	}
	
	/**
	 * Creates the cluster list used to display clusters.
	 * 
	 * This method is invoked from the constructor to create the cluster
	 * list and the {@link #popupList} that is used to display the cluster
	 * list whenever the user clicks on the {@link #dropDownButton}. This
	 * method is invoked only once from the constructor. It was primarily
	 * introduced to streamline the constructor. 
	 */
	private void createClusterListAndPopup() {
		// Create the cluster list using given memory model.
		clusterList = new JList(overlapModel);
		clusterList.setBorder(null);
		// Setup a custom cell renderer so that the colors
		// used for clusters are efficiently rendered
		clusterList.setCellRenderer(new ClusterListCellRenderer(colorMapper));
		// Setup selection mode to permit multiple selections
		clusterList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		// Add a list selection listener to update the selected entry count
		// in the text field
		clusterList.addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent event) {
				if (!event.getValueIsAdjusting()) {
					// The selection has settled.
					int selCount = 0;
					if (clusterList.getSelectedIndex() != -1) {
						selCount = clusterList.getSelectedIndices().length;
					}
					ClusterList.this.selectionStats.setText("" + selCount + " selected");
				}
			}
		});
		// Create the pop up for displaying the cluster list.
		clusterListJSP = new JScrollPane(clusterList);
		clusterListJSP.setBorder(null);
		popupList = new JPopupMenu();
		popupList.add(clusterListJSP);
	}

	/**
	 * A simple non-editable text field to display cluster selection
	 * statistics. This field is meant to provide the user with some
	 * summary information for convenience.
	 */
	private JLabel selectionStats;
	
	/**
	 * The drop down button (with a small triangle) to display cluster
	 * pop up list. This button is created indirectly from the constructor
	 * and either hides or displays the pop up that contains the cluster
	 * list. 
	 */
	private JButton dropDownButton;
	
	/**
	 * A tool bar style button that permits the user to change the colors
	 * associated with the selected clusters. This button is created in
	 * the constructor.
	 */
	private final JButton changeColor;
	
	/**
	 * A tool bar style button that permits the user to assign random
	 * colors to each cluster. The colors are chosen at random from a given
	 * palate.  This button is created in the constructor.
	 */
	private final JButton randomColors;
	
	/**
	 * The pop up that contains the cluster list to be displayed to the user.
	 * The usage of pop up menu in this context is a bit unconventional because
	 * typically menus would have only menu items. However, in this case, it has
	 * a JList. However, the pop up menu provides the convenience to display a
	 * pop up without having to handle the overheads of focus management and
	 * displaying the popup.
	 */
	private JPopupMenu popupList;
	
	/**
	 * The list that displays the list of clusters currently present in
	 * the memory model. This list is created in the constructor and is
	 * further customized in the {@link #adaptClusterListAndButton(OverlapPanel)}
	 * method. 
	 */
	private JList clusterList;
	
	/**
	 * The scroll pane that contains the {@link #clusterList} so that 
	 * the list can scroll when the there are a large number of clusters.
	 */
	private JScrollPane clusterListJSP;
	
	/**
	 * The OverlapModel that provides the necessary information
	 * for rendering the core information required by this panel.
	 * This instance variable is initialized by the constructor and
	 * is never changed during the life time of this object.
	 */
	private final OverlapModel overlapModel;
	
	/**
	 * The final top-level component that logically owns this cluster
	 * list. This component is notified to refresh the display whenver
	 * cluster colors are changed or updated.
	 */
	private final JComponent owner;
	
	/**
	 * Helper to look-up color for a given cluster ID.
	 * 
	 * This interface is used when rendering fragments. The color mapper
	 * is used to look up the colors set for a given cluster, given its
	 * cluster ID. The cluster colors are set and updated indirectly by
	 * the OverlapView class.
	 */
	private final ClusterColorMapper colorMapper;

	/**
	 * The color of the border to be used for this component to provide a 
	 * simple line border that makes it look like a static control. This
	 * value is set in the constructor before setting the background of 
	 * this control to white.
	 */
	private final Color borderColor;
	
	/**
	 * A fixed message that is displayed to the user to indicate that some of
	 * the clusters have duplicate colors. This constant is used in the
	 * {@link #actionPerformed(ActionEvent)} method.
	 */
	private static final String LOW_COLOR_MSG = "<html>" +
		"The number of default colors used to randomly color the clusters<br>" +
		"were fewer than the number of clusters. Therefore, two or more<br>" +
		"clusters have the same color.  You may choose clusters from the <br>" +
		"list and selectively assign different colors to them to suitably<br>" +
		"distinguish them from each other.<br>" +
		"</html>";
	
	/**
	 * Generated serialization UID (merely to keep the compiler happy)
	 */
	private static final long serialVersionUID = 5321838842105746304L;
}
