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
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.FileOutputStream;
import java.io.PrintStream;

import javax.swing.Box;
import javax.swing.DefaultCellEditor;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import org.peace_tools.core.MainFrame;
import org.peace_tools.core.ViewFactory;
import org.peace_tools.data.EST;
import org.peace_tools.data.ESTTableModel;
import org.peace_tools.generic.FindDialog;
import org.peace_tools.generic.FindEvent;
import org.peace_tools.generic.FindListener;
import org.peace_tools.generic.HelpHandler;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DataFileStats;

/**
 * The table model that provides a graphical display of fragments in
 * a FASTA file.
 *  
 * This class provides a tabular view of the fragments stored in a 
 * given FASTA file using a JTable. The JTable uses the data
 * exposed via the ESTTableView class. The view is created on demand
 * via the ViewFactory.  
 */
public class ESTTableView extends JPanel implements ActionListener, ChangeListener, FindListener {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor sets up the table to display FASTA 
	 * entries and configures the table to the default configuration.
	 * 
	 * @param model The data model to be used to display the FASTA entries
	 * and associated information in the table.
	 * 
	 * @param frame The main frame that logically owns this view.
	 * 
	 */
	public ESTTableView(ESTTableModel model, MainFrame frame) {
		super(new BorderLayout(0, 0));
		// Save reference to the memory model
		this.model     = model;
		this.mainFrame = frame;
		// Create a custom cell editor (to make things look pretty) for
		// the second column that shows a 'revert' button 
		JCheckBox revertCB = new JCheckBox();
		revertCB.setIcon(Utilities.getIcon("images/16x16/Blank.png"));
		revertCB.setSelectedIcon(Utilities.getIcon("images/16x16/Revert.png"));
		final DefaultCellEditor secondColEditor = new DefaultCellEditor(revertCB);
		// Create the table with custom cell renderer
		estTable = new JTable(model) {
			private static final long serialVersionUID = 6507239570257555103L;
			@Override
			public TableCellRenderer getCellRenderer(int row, int col) {
				if (col == 0) {
					return RevertBtnRenderer;
				} else if (col > 3) {
					return ESTRenderer;
				} else {
					return super.getCellRenderer(row, col);
				}
			}
			@Override
			public TableCellEditor getCellEditor(int row, int col) {
				if (col == 0) {
					return secondColEditor;
				}
				return super.getCellEditor(row, col);
			}
		};
		
		// Set some table properties
		estTable.setBorder(null);
		// clusterTable.setShowHorizontalLines(true);
		estTable.setFillsViewportHeight(true);
		estTable.setDragEnabled(false);
		estTable.setDropTarget(null);
		estTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		estTable.setGridColor(new Color(0xe0, 0xe0, 0xe0));
		estTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		// Setup column sizes
		setColumnWidths();
		// Prevent reordering of columns
		estTable.getTableHeader().setReorderingAllowed(false);
		
		// Place the table in a scroll pane to it can scroll
		JScrollPane scroller = new JScrollPane(estTable);
		scroller.setBorder(null);
		scroller.getViewport().setBackground(estTable.getBackground());
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
	 * This method implements the find interface to use the find dialog. This 
	 * method intercepts the events fired by the find dialog whenever the
	 * user clicks the 'find' button (in that dialog). This method 
	 * performs a search (with the help of the model) and appropriately
	 * updates the currently selected row.
	 * 
	 * @param event The find event that contains information about the search
	 * requested by the user.
	 */
	@Override
	public boolean find(FindEvent event) {
		// If no row is selected then the selected row is -1 but our search will
		// commence at row 0 by default.
		final int currRow = Math.max(0, estTable.getSelectedRow());
		// Get the model to do the actual search.
		int matchingRow = model.find(event, currRow);
		// Check and update the currently selected row.
		if (matchingRow != -1) {
			// Scroll to the matching row
		    Rectangle rect = estTable.getCellRect(matchingRow, 0, true);
		    estTable.scrollRectToVisible(rect);
		    // Set the row to be selected.
			estTable.getSelectionModel().setSelectionInterval(matchingRow, matchingRow);
			return true;
		}
		// No match found.
		return false;
	}
	
	/**
	 * Utility method to set the column widths.
	 * 
	 * This is a utility method that can be used set up default column sizes
	 * for the various columns displayed in this table.
	 */
	private void setColumnWidths() {
		// Setup some column properties.
		final TableColumnModel tcm = estTable.getColumnModel();
		final int ColCount   = tcm.getColumnCount();
		final int ColSizes[] = {20, 20, 75, 150, (ColCount > 3) ? 250 : 800};
		// Setup column sizes for each column.
		for(int colId = 0; (colId < tcm.getColumnCount()); colId++) {
			tcm.getColumn(colId).setPreferredWidth((colId > 3) ? ColSizes[4] : ColSizes[colId]);
		}		
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
				"Text View", "text", this, 
				"Displays raw text version of FASTA file", 
				(model.getESTList().getName() != null))); // Enable button only if file name
		// Careate and add combo box to select the type of sorting to be used.
		final String SortOptions[] = {
			"None", "Shortest fragments first", "Shortest fragments last"
		};
		sortOrder = new JComboBox(SortOptions);
		sortOrder.setEditable(false);
		sortOrder.setSelectedIndex(0);
		sortOrder.addActionListener(this);
		sortOrder.setActionCommand("resort");
		sortOrder.setMaximumSize(sortOrder.getPreferredSize());
		// Add combo box to the tool bar.
		toolbar.add(Box.createHorizontalStrut(10));
		toolbar.add(new JLabel("Sort entries: "));
		toolbar.add(sortOrder);
		// Add button to save selected ESTs.
		toolbar.add(Box.createHorizontalStrut(10));
		JButton findButton = Utilities.createToolButton("images/16x16/Find.png",
				"", "find", this, 
				"Search for information in the FASTA file", true);
		findButton.setMnemonic(KeyEvent.VK_F);
		toolbar.add(findButton);
		toolbar.add(Box.createHorizontalStrut(5));
		toolbar.add(Utilities.createToolButton("images/16x16/SaveEST.png",
				"", "save", this, 
				"Saves selected entries to another FASTA file", true));
		
		// Add button to enable/disable column-wise display and set bases per column.
		columization = new JCheckBox("Multiple columns with", false);
		columization.setActionCommand("columization");
		columization.addActionListener(this);
		columization.setOpaque(false);
		toolbar.add(columization);
		SpinnerNumberModel model = new SpinnerNumberModel(100, 20, 500, 20); 
		// Create spinner for bases/column and add to tool bar
		basesPerCol = new JSpinner(model);
		toolbar.add(basesPerCol);
		basesPerCol.addChangeListener(this);
		basesPerCol.setEnabled(false);
		// Fix the maximum width of this field to make thing prettier
		Dimension maxSize = basesPerCol.getPreferredSize();
		maxSize.width     = maxSize.width * 3 / 2;
		basesPerCol.setPreferredSize(maxSize);
		basesPerCol.setMaximumSize(maxSize);
		// Put wrapping message.
		toolbar.add(new JLabel(" nt/column"));
		
		toolbar.add(Box.createHorizontalStrut(10));
		// Create the help button to provide some help.
		toolbar.add(Box.createHorizontalStrut(10));
		toolbar.add(Utilities.createToolButton("images/16x16/Help.png", 
				null, "help", this, 
				"Read about the the FASTA display and various controls", true));
		// Create a horizontal expandable strut
		toolbar.add(Box.createHorizontalGlue());
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
	private JLabel createSummaryInfo() {
		// Obtain the two distinct summary information.
		DataFileStats fastaStats   = model.getESTList().computeStatistics();
		// Build the summary information into a string 
		String summary = String.format(" #Entries: %d, Shortest fragment: %d nt, " +
				"Longest fragment: %d nt, Avg. fragment length: %.2f nt (SD: %.2f nt)", 
				fastaStats.getCount(), fastaStats.getMinLength(), 
				fastaStats.getMaxLength(), fastaStats.getAvgLength(), 
				fastaStats.getLengthSD());
		// Now store all the summary information into a JLabel.
		JLabel infoLabel = new JLabel(summary);
		infoLabel.setIcon(Utilities.getIcon("images/16x16/Information.png"));
		infoLabel.setBorder(new EmptyBorder(2, 10, 2, 10));
		infoLabel.setBackground(Color.white);
		return infoLabel;
	}
	
	@Override
	public void actionPerformed(ActionEvent event) {
		if ("text".equals(event.getActionCommand())) {
			String fileName = model.getESTList().getName();
			mainFrame.getViewFactory().createView(fileName, fileName, null, ViewFactory.ViewType.TEXT_VEIW, false, true, null);
		} else if ("save".equals(event.getActionCommand())) {
			saveSelectedESTs();
		} else if ("resort".equals(event.getActionCommand())) {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			model.sort(sortOrder.getSelectedIndex());
			setCursor(Cursor.getDefaultCursor());
			repaint();
		} else if ("columization".equals(event.getActionCommand())) {
			basesPerCol.setEnabled(columization.isSelected());
			stateChanged(null);
		} else if ("help".equals(event.getActionCommand())) {
			HelpHandler.showHelp(mainFrame, "http://www.peace-tools.org/downloads/manual.pdf#page=34");
		} else if ("find".equals(event.getActionCommand())) {
			FindDialog fd = FindDialog.getDialog();
			fd.setFindListener(this);
			fd.setVisible(true);
		}
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		// The user changed the bases per column value. Update model.
		int colSize = -1;
		if (this.columization.isSelected()) {
			Integer value = (Integer) basesPerCol.getValue();
			colSize = value;
		}
		model.setBasesPerCol(colSize);
		// Update preferred column sizes.
		setColumnWidths();
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
		int[] selectedRows = model.getSelectedRows();
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
			FileOutputStream fos      = new FileOutputStream(fileName);
			PrintStream      ps       = new PrintStream(fos);
			for(int i = 0; (i < selectedRows.length); i++) {
				EST est = model.getESTAt(selectedRows[i]);
				est.write(ps);
			}
			// Close the stream.
			ps.close();
			// report successful file save to the user.
			UserLog.log(UserLog.LogLevel.NOTICE, "ClusterTreeTableView",
					"Selected ESTs have been saved to FASTA file " + fileName);
		} catch (Exception e) {
			// Format and display exception
			JPanel msg = Utilities.collapsedMessage("Unable to export ESTs to specified FASTA file.",
					Utilities.toString(e));
			JOptionPane.showMessageDialog(this, msg, "Save error", JOptionPane.ERROR_MESSAGE);
		}
	}

	/**
	 * The icon that is visible and looks similar to the revert button image used
	 * by the corresponding cell editor. This image is loaded once and reused
	 */
	private static final ImageIcon RevertIcon = Utilities.getIcon("images/16x16/Revert.png");

	/**
	 * The icon that is invisible and merely takes up the same space as the Revert
	 * icon (used when no reverting is needed). This image is loaded once and reused
	 */
	private static final ImageIcon BlankIcon = Utilities.getIcon("images/16x16/Blank.png");

	
	private static class RevertCellRenderer extends JLabel
    implements TableCellRenderer {
		/**
		 * The general serial version UID to enable serialization of this class as
		 * per Java requirements.
		 */
		private static final long serialVersionUID = 5654108539980884223L;
		
		
		/** Default constructor.
		 * 
		 * The default constructor merely sets the standard information in
		 * the label.
		 */
		public RevertCellRenderer() {
			super(BlankIcon);
			setOpaque(true);
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			Boolean data = (Boolean) value;
			setIcon(data ? RevertIcon : BlankIcon);
			setBackground(isSelected ? table.getSelectionBackground() : table.getBackground());
			return this;
		}
	}
	
	/**
	 * The toolbar that contains some commonly used tools with 
	 * the jobs.
	 */
	private JToolBar toolbar;

	/**
	 * The actual Tree-table that provides a graphical view of the
	 * list of ESTs in a given data file.
	 */
	private JTable estTable;

	/**
	 * The table model that is being used to display information
	 * in this view. This value is set in the constructor and is
	 * never changed.
	 */
	private final ESTTableModel model;
	
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
	 * A checkbox that can be used by the user to enable or disable 
	 * column-wise display of nucleotide sequences. 
	 */
	private JCheckBox columization;
	
	/**
	 * A spinner that is displayed in the toolbar. This spinner can be used
	 * by the user to set the number of bases/nucleotides to be displayed
	 * in a single column in the table.
	 */
	private JSpinner basesPerCol;
	
	/**
	 * A simple/custom cell renderer that displays data using a monospaced
	 * font. This cell renderer is used to render all the sequence data so
	 * that entries in multiple rows appear to align consistently on the
	 * same boundary.
	*/
	private static final CustomTableCellRenderer ESTRenderer = new CustomTableCellRenderer();
	
	/**
	 * A simple/custom cell renderer that displays a blank or revert button
	 * when a specific EST entry has been modified (or is unmodified). This renderer
	 * is used to display the data in the second column of the table display.
	 */
	private static final RevertCellRenderer RevertBtnRenderer = new RevertCellRenderer();

	/**
	 * A simple message that is displayed to the user. 
	 */
	private String SELECT_EST_MSG = "<html>" +
			"You must first select one or more entries to be exported to another<br>" +
			"FASTA file. You can select ranges using the control and shift keys<br>" +
			"while selecting entries.</html>";
	
	/**
	 * A generated serialization UID (included just to keep the compiler happy).
	 */
	private static final long serialVersionUID = 6547127105109078286L;
}
