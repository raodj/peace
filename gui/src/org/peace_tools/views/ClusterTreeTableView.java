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
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.views;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.Icon;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.table.TableColumnModel;

import org.netbeans.swing.etable.ETableColumn;
import org.netbeans.swing.outline.DefaultOutlineModel;
import org.netbeans.swing.outline.Outline;
import org.netbeans.swing.outline.OutlineModel;
import org.netbeans.swing.outline.RenderDataProvider;
import org.peace_tools.core.MainFrame;
import org.peace_tools.data.ClusterNode;
import org.peace_tools.data.ClusterTreeTableModel;
import org.peace_tools.data.ESTList;
import org.peace_tools.generic.HelpHandler;
import org.peace_tools.generic.Pair;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;

/**
 * The tree table model that provides a combination of a tree and table
 * to display cluster information in a convenient form.
 *  
 * <p>The tree-table or Outline model is a generic non-standard GUI component
 * developed by Sun Microsystems as part of the NetBeans IDE. This jar has been
 * obtain from the NetBeans package as a part. The tree-table is a combination 
 * both a JTree and a JTable. The first column in the tree-table is a tree that
 * provides the user with a convenient interface to access and control view of
 * hierarchical information. The remaining columns in the tree-table display 
 * detailed information about each entry in the tree table. Additional details
 * on the tree table is available at: 
 * http://bits.netbeans.org/dev/javadoc/index.html</p> 
 */
public class ClusterTreeTableView extends JPanel implements ActionListener {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor sets up the job list table and 
	 * configures the tree-table to the default configuration.
	 * 
	 * @param model The data model to be used to display the cluster
	 * and EST information in the tree-table.
	 * 
	 * @param frame The main frame that logically owns this view. This
	 * reference is used to launch the classifier dialog.
	 */
	public ClusterTreeTableView(ClusterTreeTableModel model, MainFrame frame) {
		super(new BorderLayout(0, 0));
		// Save reference to the memory model
		this.model     = model;
		this.mainFrame = frame;
		// First create a default outline model
		OutlineModel outModel = 
			DefaultOutlineModel.createOutlineModel(model, model);
		// Create the tree-table
		clusterTable = new Outline(outModel);
		// Set some table properties
		clusterTable.setBorder(null);
		// clusterTable.setShowHorizontalLines(true);
		clusterTable.setRenderDataProvider(new DataSetRenderer()); 
		clusterTable.setFillsViewportHeight(true);
		clusterTable.setDragEnabled(false);
		clusterTable.setDropTarget(null);
		clusterTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		clusterTable.setGridColor(new Color(0xe0, 0xe0, 0xe0));
		clusterTable.setRootVisible(false);
		// Setup some column properties.
		TableColumnModel tcm = clusterTable.getColumnModel();
		for(int colId = 0; (colId < tcm.getColumnCount()); colId++) {
			ETableColumn column = (ETableColumn) tcm.getColumn(colId);
			column.setPreferredWidth(200);
			if (colId > 0) {
				// Set a default cell renderer with mono-spaced font.
				column.setCellRenderer(new CustomTableCellRenderer());
			}
		}
		// Place the tree-table in a scroll pane to it can scroll
		JScrollPane scroller = new JScrollPane(clusterTable);
		scroller.setBorder(null);
		scroller.getViewport().setBackground(clusterTable.getBackground());
		// Create summary information.
		Component summaryPanel = createSummaryInfo();
		// Place the summary information and the tree-table into a panel.
		JPanel subPanel = new  JPanel(new BorderLayout(0, 0));
		subPanel.add(summaryPanel, BorderLayout.NORTH);
		subPanel.add(scroller, BorderLayout.CENTER);
		// Add sub panel to ourselves at the center to it can take up all
		// the space we have.
		add(subPanel, BorderLayout.CENTER);
		// Create the toolbar at the top of the job list
		createToolbar();
	}
	
	/**
	 * Helper method to configure and setup various tools in the
	 * tool bar. This method was introduced to streamline the code
	 * better and cut down the code clutter in the constructor.
	 */
	private void createToolbar() {
		toolbar = new JToolBar();
		toolbar.setFloatable(false);
		// Add button to launch the DB classifier dialog.
		toolbar.add(Utilities.createToolButton("images/16x16/ClusterSummary.png",
				"Summary Graph", "graph", this, 
				"Displays summary information about clusters as a graph", true));
		// Careate and add combo box to select the type of sorting to be used.
		final String SortOptions[] = {
			"None", "Smallest clusters first", "Smallest clusters last"
		};
		sortOrder = new JComboBox(SortOptions);
		sortOrder.setEditable(false);
		sortOrder.setSelectedIndex(0);
		sortOrder.addActionListener(this);
		sortOrder.setActionCommand("resort");
		sortOrder.setMaximumSize(sortOrder.getPreferredSize());
		// Add combo box to the tool bar.
		toolbar.add(Box.createHorizontalStrut(10));
		toolbar.add(new JLabel("Sort clusters: "));
		toolbar.add(sortOrder);
		// Add button to save selected ESTs.
		toolbar.add(Box.createHorizontalStrut(10));
		toolbar.add(Utilities.createToolButton("images/16x16/SaveEST.png",
				"", "save", this, 
				"Saves ESTs in currently selected clusters to a given file", true));
		// Add button to show or hide summary information
		toolbar.add(Box.createHorizontalStrut(10));
		JToggleButton summaryButton = 
			new JToggleButton("Summary ...", 
					Utilities.getIcon("images/16x16/Information.png"), true);
		summaryButton.setToolTipText("Show/Hide summary information about this cluster file");
		summaryButton.addActionListener(this);
		summaryButton.setActionCommand("summary");
		// Make button a tool bar button.
		Utilities.makeToolBarButton(summaryButton, true);
		toolbar.add(summaryButton);
		// Create the help button to provide some help.
		toolbar.add(Box.createHorizontalStrut(10));
		toolbar.add(Utilities.createToolButton("images/16x16/Help.png", 
				null, "help", this, 
				"Read about the cluster display and various controls", true));
		// Add tool bar to the north
		add(toolbar, BorderLayout.NORTH);
	}

	/**
	 * Method to compute and setup summary information.
	 * 
	 * This method is invoked from the constructor to create and populate the 
	 * summary information tab. The summary information is computed once when
	 * the view is created. After that the summary information (which is 
	 * reasonably small) is held in a text area. The user has the option
	 * to hide the summary information to gain some more screen real estate.
	 * 
	 * @return The scroll pane that contains the summary information.
	 */
	private JScrollPane createSummaryInfo() {
		// Obtain the two distinct summary information.
		double[] fastaStats   = model.getESTList().computeStatistics();
		double[] clusterStats = model.getRoot().computeStatistics();
		assert ( fastaStats.length == clusterStats.length );
		// Build the summary information into a string buffer first.
		StringBuilder sb = new StringBuilder();
		for(int idx = 0; (idx < SUMMARY_STATS_STRING.length); idx++) {
			// Account for unformatted header (3 lines) and trailer (last line)
			if ((idx > 2) && (idx < SUMMARY_STATS_STRING.length - 1)) {
				sb.append(String.format(SUMMARY_STATS_STRING[idx], fastaStats[idx - 3], 
						clusterStats[idx - 3]));
			} else {
				// Unformatted static string information
				sb.append(SUMMARY_STATS_STRING[idx]);
			}
		}
		// Next add general meta data from the cluster file as well
		sb.append("\nCLUSTERING JOB SUMMARY INFORMATION\n");
		ArrayList<Pair> metadata = model.getClusterFile().getMetadata();
		for(Pair info: metadata) {
			sb.append(String.format("%s: %s\n", info.getName(), info.getValue()));
		}
		// Now store all the summary information into a text area.
		JTextArea summaryInfo = new JTextArea(7, 50);
		summaryInfo.setBorder(null);
		summaryInfo.setEditable(false);
		// summaryInfo.setBackground(new Color(0xe0, 0xe0, 0xe0));
		// Make the font mono spaced.
		int fontSize = summaryInfo.getFont().getSize();
		summaryInfo.setFont(new Font(Font.MONOSPACED, Font.PLAIN, fontSize));
		// Store information to be displayed.
		summaryInfo.append(sb.toString());
		summaryInfo.setCaretPosition(0);
		// Wrap the text area into a scroll pane to enable scrolling.
		summaryPane = new JScrollPane(summaryInfo);
		return summaryPane;
	}
	
	@Override
	public void actionPerformed(ActionEvent event) {
		if ("graph".equals(event.getActionCommand())) {
			mainFrame.getViewFactory().createSummaryView(model.getWsEntry());
		} else if ("save".equals(event.getActionCommand())) {
			saveSelectedESTs();
		} else if ("resort".equals(event.getActionCommand())) {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			model.sort(sortOrder.getSelectedIndex());
			setCursor(Cursor.getDefaultCursor());
			repaint();
		} else if ("help".equals(event.getActionCommand())) {
			HelpHandler.showHelp(mainFrame, "http://www.peace-tools.org/downloads/manual.pdf#page=34");
		} else if ("summary".equals(event.getActionCommand())) {
			summaryPane.setVisible(!summaryPane.isVisible());
			revalidate();
			repaint();
		}
	}
	
	/**
	 * Helper method to save the selected ESTs.
	 * 
	 * This method is invoked from the actionPerformed() method whenever
	 * the user clicks on the tool bar button to save selected ESTs. This
	 * method performs the following tasks:
	 * 
	 * <ol>
	 * 
	 * <li>If no entries are selected, then this method displays a warning
	 * message and exits immediately with no further action.</li>
	 * 
	 * <li>Otherwise, it prompts the user to select a file. If the user
	 * cancels the file selection operation then this method exits without
	 * any side effects.</li>
	 * 
	 * <li>Once the user selects a valid file, it writes the ESTs 
	 * corresponding to the selected entries to the file. Any errors
	 * that occur are suitably reported to the user.</li>
	 * 
	 * </ol>
	 */
	private void saveSelectedESTs() {
		int[] selectedRows = clusterTable.getSelectedRows();
		if ((selectedRows == null) || (selectedRows.length == 0)) {
			JOptionPane.showMessageDialog(this, SELECT_EST_MSG,
					"Selection is empty", JOptionPane.INFORMATION_MESSAGE);
			return;
		}
		// Now that the selection is valid, let the user choose a valid
		// file to save the data to.
        JFileChooser jfc = new JFileChooser();
        jfc.setDialogTitle("Choose target FASTA file");
        jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        jfc.setDialogType(JFileChooser.SAVE_DIALOG);
        if (jfc.showDialog(this, "Save (FASTA)") != JFileChooser.APPROVE_OPTION) {
        	// The user canceled the file selection. bail out
        	return;
        }
        // Write the selected ESTs to the file.
        String fileName = jfc.getSelectedFile().getAbsolutePath();
        try {
			FileOutputStream fos = new FileOutputStream(fileName);
			PrintStream      ps  = new PrintStream(fos);
			for(int i = 0; (i < selectedRows.length); i++) {
				Object value = clusterTable.getValueAt(selectedRows[i], 0);
				if (value instanceof ClusterNode) {
					// The nodes recursively write ESTs if needed.
					ClusterNode node = (ClusterNode) value;
					node.write(model.getESTList(), ps);
				}
			}
			// Close the stream.
			ps.close();
			// report successful file save to the user.
			UserLog.log(UserLog.LogLevel.NOTICE, "ClusterTreeTableView",
					"Selected ESTs have been saved to FASTA file " + fileName);
		} catch (Exception e) {
			JPanel msg = Utilities.collapsedMessage("Unable to export ESTs to specified FASTA file.",
					Utilities.toString(e));
			
			JOptionPane.showMessageDialog(this, msg, "Save error", JOptionPane.ERROR_MESSAGE);
		}
	}
	
	/**
	 * A tree cell renderer that essentially provides a better 
	 * representative set of Icons to make the overall display 
	 * look a bit prettier.
	 */
	private class DataSetRenderer implements RenderDataProvider {
		/**
		 * The default constructor. Currently it does nothing
		 * and is here merely to adhere to coding conventions.
		 */
		public DataSetRenderer() {
			// Nothing to be done for now.
		}

		/**
		 * Provide background for a given object. Provides a default
		 * implementation for the RenderDataProvider API.
		 * 
		 * @return This method always returns null for default
		 * (or no) background color. 
		 */
		@Override
		public java.awt.Color getBackground(Object o) {
			return null;
		}

		/**
		 * Provides a name for the node to be displayed on the tree. 
		 * 
		 * Provides a default implementation for the RenderDataProvider API.
		 * 
		 * @return This method returns a suitable string corresponding
		 * to the given object.  
		 */
		@Override
		public String getDisplayName(Object o) {
			ClusterNode node = (ClusterNode) o;
			if (node.isESTNode()) {
				final int estID    = node.getESTId();
				final ESTList ests = model.getESTList();
				String name = "EST #" + node.getESTId() + " (" +
				ests.getESTs().get(estID).getInfo() + ")";
				return "<html>" + name + "</html>";
			} else {
				return "<html><i>" + o.toString() + "</i></html>";
			}
		}

		/**
		 * Provide foreground for a given object. Provides a default
		 * implementation for the RenderDataProvider API.
		 * 
		 * @return This method always returns null for default
		 * (or no) foreground color. 
		 */
		@Override
		public java.awt.Color getForeground(Object o) {
			return null;
		}

		/**
		 * Provide an icon for a given object. Provides a default
		 * implementation for the RenderDataProvider API.
		 * 
		 * @return This method returns a suitable icon for a the
		 * given object assuming that the object is either an EST
		 * node or a Cluster. 
		 */
		public javax.swing.Icon getIcon(Object o) {
			ClusterNode node = (ClusterNode) o;
			return TreeIcons[node.isESTNode() ? 1 : 0];
		}

		/**
		 * Provide an tool tip for the node. Provides a default
		 * implementation for the RenderDataProvider API.
		 * 
		 * @return This method currently returns null for all 
		 * objects.
		 */
		public String getTooltipText(Object o) {
			return null;
		}

		/**
		 * Indicate all labels are in HTML. Provides a default
		 * implementation for the RenderDataProvider API.
		 * 
		 * @return This method currently returns true to indicate
		 * that the labels are always in HTML.
		 */
		public boolean isHtmlDisplayName(Object o) {
			return true;
		}
	}

	/**
	 * The array of icons that are displayed in a the cluster/EST
	 * tree adjacent to entries. The icons are meant to provide
	 * quick visual cues as to the contents of an tree and possibly
	 * make the display look customized and nice.
	 */
	private static final Icon TreeIcons[] = { 
			Utilities.getIcon("images/16x16/ClusterSmall.png"),
			Utilities.getIcon("images/16x16/ESTSmall.png")
	};
	
	/**
	 * The toolbar that contains some commonly used tools with 
	 * the jobs.
	 */
	private JToolBar toolbar;

	/**
	 * The actual Tree-table that provides a graphical view of the
	 * list of clusters in a given cluster file.
	 */
	private Outline clusterTable;

	/**
	 * The table-tree model that is being used to display information
	 * in this view. This value is set in the constructor and is
	 * never changed.
	 */
	private final ClusterTreeTableModel model;
	
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this component. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;

	/**
	 * Various options on how the clusters in the cluster tree table view
	 * are to be sorted. The options are populated when this combo box is
	 * created in the createToolbar() method.
	 */
	private JComboBox sortOrder;
	
	/**
	 * The scroll pane that contains a text area that displays summary
	 * statistics regarding ESTs and the clusters displayed in this 
	 * view. The summary pane is created by the constructor and can 
	 * be shown/hidden by the user by clicking on a toggle button.
	 */
	private JScrollPane summaryPane;
	
	/**
	 * A simple message that is displayed to the user. 
	 */
	private String SELECT_EST_MSG = "<html>" +
			"You must first select one or more ESTs or Clusters to be<br>" +
			"exported to another FASTA file. You can select ranges using the<br>" +
			"control and shift keys while selecting ESTs.</html>";
	
	/**
	 * Various format strings that are formatted (to fill-in values) to
	 * display summary statistics to the user. These strings are used in
	 * the createSummaryInfo(). These strings are pulled out of the method
	 * to keep the code cleaner in the method. 
	 */
	private final String SUMMARY_STATS_STRING[] = {
		"---------------------------------+---------------------------------\n",
		" FASTA SUMMARY STATISTICS        |       CLUSTER SUMMARY STATISTICS\n",
		"---------------------------------+---------------------------------\n",
		" Fragment count      : %-6.0f    |     Cluster count     : %-6.0f\n",
		" Min fragement length: %-6.0f    |     Min cluster size  : %-6.0f\n",
		" Max fragement length: %-6.0f    |     Max cluster size  : %-6.0f\n",
		" Avg fragement length: %-6.2f    |     Avg cluster size  : %-6.2f\n",
		" Standard deviation  : %-6.2f    |     Standard deviation: %-6.2f\n",
		"---------------------------------+---------------------------------\n"
	};
	
	/**
	 * A generated serial version ID for serialization (more
	 * realistically to keep the compiler happy). 
	 */
	private static final long serialVersionUID = 80617431851108817L;
}
