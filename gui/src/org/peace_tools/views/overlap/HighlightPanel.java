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
import java.awt.Rectangle;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.TableColumnModel;

import org.peace_tools.data.ESTList;
import org.peace_tools.data.ESTTableModel;
import org.peace_tools.generic.CustomBorder;
import org.peace_tools.generic.Utilities;


/**
 * A panel that lists ESTs to be highlighted.
 * 
 * This panel (can be shown and hidden by the user) permits the user
 * to quickly select ESTs to be highlighted and highlight them. This
 * panel utilizes an existing table model for displaying ESTs to be
 * highlighted. 
 */
class HighlightPanel extends JPanel implements ChangeListener {
	HighlightPanel(ESTList estList) {
		super(new BorderLayout(0, 5));
		// Create and add the GUI components to a sub-panel
		JPanel subPanel = new JPanel(new BorderLayout(0, 5));
		subPanel.add(createGotoPanel(estList.getESTs().size()), BorderLayout.NORTH);
		subPanel.add(createTable(estList), BorderLayout.CENTER);
		// The top-panel contains a title and sub-panel that contains
		// various other components.
		JLabel heading = new JLabel("<html><b>Highlight Panel</b></html>",
				Utilities.getIcon("images/16x16/Information.png"), JLabel.LEFT);
		heading.setBorder(BorderFactory.createCompoundBorder(
				new CustomBorder("ssds"), 
				BorderFactory.createEmptyBorder(2, 5, 2, 5)));
		add(heading, BorderLayout.NORTH);
		add(subPanel, BorderLayout.CENTER);
	}

	/**
	 * Obtain the table model used by the EST table in this panel.
	 * 
	 * @return The table model used by the EST table in this panel.
	 */
	public ESTTableModel getModel() {
		return this.model;
	}
	
	/**
	 * Obtain the GUI table component used in this panel.
	 * 
	 * @return The GUI table used in this panel to list the ESTs.
	 */	
	public JTable getTable() {
		return this.estTable;
	}
	
	/**
	 * Method to handle user changes to the EST index spinner. This method 
	 * is invoked whenever the user changes the entry in the
	 * {@link #estIndexEntry} spinner. This method changes the currently
	 * selected row in the table to the value in the spinner (assuming
	 * that the value in the spinner is valid).
	 * 
	 * @param e The event associated with the user action. This event is
	 * currently ignored.
	 */
	@Override
	public void stateChanged(ChangeEvent e) {
		int row = ((Number) estIndexEntry.getValue()).intValue();
		if ((row >= 0) && (row < model.getRowCount())) {
			// Row is valid. First scroll the table to the appropriate row
		    Rectangle rect = estTable.getCellRect(row, 0, true);
		    estTable.scrollRectToVisible(rect);
			this.estTable.getSelectionModel().setSelectionInterval(row, row);
		}
	}
	
	/**
	 * Helper method to create the panel that has a brief message followed by a 
	 * spinner that enables the user to rapidly navigate to a specific EST entry 
	 * given its index.
	 * 
	 * @param estCount The total number of ESTs being displayed in the model.
	 * 
	 * @return This method returns a panel containing a quick help message along
	 * with a spinner.
	 */
	private JPanel createGotoPanel(final int estCount) {
		// Create the spinner to quickly locate entries
		estIndexEntry = new JSpinner(new SpinnerNumberModel(0, 0, estCount - 1, 1));
		estIndexEntry.addChangeListener(this);
        Utilities.adjustDimension(estIndexEntry, 10, 4); // Adjust size to look right
		// Place the spinner with suitable labels into a panel.
		JPanel spinPanel = Utilities.createLabeledComponents("Enter EST index to select:", 
				"(Helps navigate entries below)", 0, false, estIndexEntry);
		JLabel msg = new JLabel(INFO_MSG);
		msg.setBorder(BorderFactory.createTitledBorder(""));
		// Now place all the components in a suitable panel.
		JPanel topPanel = new JPanel(new BorderLayout(5, 5));
		topPanel.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 5));
		topPanel.add(msg, BorderLayout.NORTH);
		topPanel.add(spinPanel, BorderLayout.SOUTH);
		return topPanel;
	}
	
	/**
	 * Helper method to create the table that lists all the ESTs. This 
	 * method is called only once from the constructor. This method was introduced
	 * to streamline the constructor and keep the code clutter to a minimum.
	 * 
	 * @return This method returns a scroll pane containing the EST table.
	 */
	private JComponent createTable(ESTList estList) {
		// Create the table model to be used for the EST list.
		model = new ESTTableModel(estList, true);
		// Create the table with custom model
		estTable = new JTable(model);
		// Set some table properties
		estTable.setBorder(null);
		estTable.setShowHorizontalLines(true);
		estTable.setShowVerticalLines(false);
		estTable.setFillsViewportHeight(true);
		estTable.setDragEnabled(false);
		estTable.setDropTarget(null);
		estTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		estTable.setGridColor(new Color(0xe0, 0xe0, 0xe0));
		estTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		estTable.getSelectionModel().setSelectionInterval(0, 0);
		// Prevent reordering of columns and rows
		estTable.getTableHeader().setReorderingAllowed(false);
		estTable.setAutoCreateRowSorter(false);
		// Setup column sizes
		final TableColumnModel tcm = estTable.getColumnModel();
		tcm.getColumn(0).setPreferredWidth(20);
		tcm.getColumn(1).setPreferredWidth(75);
		tcm.getColumn(2).setPreferredWidth(125);
		
		// Place the table in a scroll pane to it can scroll
		JScrollPane scroller = new JScrollPane(estTable);
		scroller.getViewport().setBackground(estTable.getBackground());
		// Return the scroller for use in the main panel.
		return scroller;
	}
	
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
	private ESTTableModel model;

	/**
	 * This spinner permits the user to quickly navigate to a given EST
	 * entry using the index of the EST entry. This spinner is created in
	 * the {@link #createGotoPanel(int)} method.
	 */
	private JSpinner estIndexEntry;
	
    /**
     * A simple informational message that is displayed to the user. 
     */
	private static final String INFO_MSG = "<html><font size=\"-2\"><i>"     +
		"Select the ESTs to be hilited "         +
		"by checking the box next to an " +
		"EST entry in the table below. " +
		"Highlighted ESTs flash. "  +
		"You may use the spinner below "+
		"to locate ESTs based on their " +
		"logical index number."             +
		"</i></font></html>";
	
	/**
	 * Generated serialization GUID just to keep compiler happy. 
	 */
	private static final long serialVersionUID = -6499741676563416904L;

}