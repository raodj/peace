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
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import javax.swing.AbstractCellEditor;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.DefaultCellEditor;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.ListSelectionModel;
import javax.swing.border.Border;
import javax.swing.event.CellEditorListener;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableColumnModel;

import org.peace_tools.core.MainFrame;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.ClassifierList;
import org.peace_tools.workspace.DBClassifier;
import org.peace_tools.workspace.Workspace;

/**
 * This class provides a comprehensive editor for modifying the 
 * DBClassifierEditor. 
 */
public class DBClassifierEditor extends JDialog
	implements ActionListener, ListSelectionListener, CellEditorListener {
	/**
	 * The constructor.
	 * 
	 * The constructor initializes the base class and then creates the
	 * various controls for editing the current list of database 
	 * classifiers associated with this work space.
	 * 
	 * @param parent The parent window to which this editor dialog logically
	 * belongs.
	 */
	protected DBClassifierEditor(MainFrame parent) {
		super(parent, "Data base Classifier Editor", true); // creates a modal dialog.
		this.mainFrame = parent;
		// Set the preferred layout to be border layout.
		setLayout(new BorderLayout(0, 0));
		// Setup the overall preferred size
		setPreferredSize(new Dimension(500, 300));
		// Create and setup the tool bar to the north of this component
		toolbar = new JToolBar();
		createTools();
		add(toolbar, BorderLayout.NORTH);
		// Create the actual table.
		createTable();
		// Finally create the OK, Apply, and Cancel buttons.
		JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 10, 7));
		buttonPanel.add(Utilities.createButton(null, "  OK  ", 
				"ok", this, "Apply all changes and close this dialog", true));
		buttonPanel.add(Utilities.createButton(null, " Apply ", 
				"apply", this, "Apply all changes and continue editing", true));
		buttonPanel.add(Utilities.createButton(null, " Cancel ", 
				"cancel", this, "Discard all changes and close this dialog", true));
		// Add button panel to main pane.
		add(buttonPanel, BorderLayout.SOUTH);
		// Setup closing listeners for the dialog.
		setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		this.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				DBClassifierEditor.this.checkAndClose();
			}
		});
	}

	/**
	 * Helper method to configure and setup various tools in the
	 * tool bar. This method was introduced to streamline the code
	 * better and cut down the code clutter in the constructor.
	 */
	private void createTools() {
		toolbar.setFloatable(false);
		// Create tool to add and delete classifier entries
		toolbar.add(Utilities.createToolButton("images/16x16/AddClassifier.png", 
				"", "add", this, "Add a new database classifier entry", true));
		toolbar.add(Utilities.createToolButton("images/16x16/DeleteClassifier.png", 
				"", "delete", this, "Delete currently selected database classifier entry", false));
		toolbar.add(Box.createHorizontalStrut(10));

		// Create tools to move entries up and down in the list.
		toolbar.add(Utilities.createToolButton("images/16x16/UpArrow.png", 
				"", "up", this, "Move current entry higher up in the list", false));
		toolbar.add(Utilities.createToolButton("images/16x16/DownArrow.png", 
				"", "down", this, "Move current entry lower down in the list", false));
		toolbar.add(Box.createHorizontalStrut(10));

		// Create tools for checking entries and for help.
		toolbar.add(Utilities.createToolButton("images/16x16/CheckClassifier.png", 
				"", "check", this, "Check the strings this classifier matches", false));
		toolbar.add(Utilities.createToolButton("images/16x16/RevertClassifier.png", 
				"", "revert", this, "Revert data base classifiers back to defaults", true));
		toolbar.add(Utilities.createToolButton("images/16x16/Help.png", 
				"", "help", this, "Help on DB classifiers and this editor", true));
	}

	/**
	 * Helper method to configure and setup various columns in the
	 * DB classifier list table. This method was introduced to 
	 * streamline the code better and cut down the code clutter 
	 * in the constructor.
	 */
	private void createTable() {
		// Create and populate the table with the current set of 
		// entries.
		final String[] ColTitles =  {"Description", 
				"Regular Expression", "Color", "Enabled"};
		tableData = new DefaultTableModel(ColTitles, 0) {
			private static final long serialVersionUID = -1523552352906125456L;
			/**
			 * Override default class type of column three to be boolean.
			 * 
			 * The data type class for the third column needs to be returned
			 * as a boolean to ensure that the JTable displays the data
			 * as convenient check boxes rather than as plain strings.
			 * 
			 * @return The data type class to be used for a given column.
			 */
			@Override
			public Class<?> getColumnClass(int columnIndex) {
				return (columnIndex == 3) ? Boolean.class : 
					super.getColumnClass(columnIndex);
			}
		};
		// Populate the data values.
		setTableData();
		
		table = new JTable(tableData);
		// Set some table properties
		table.setBorder(null);
		table.setShowHorizontalLines(true);
		table.setFillsViewportHeight(true);
		table.setDragEnabled(false);
		table.setDropTarget(null);
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		table.getSelectionModel().addListSelectionListener(this);
		// Ensure table rows are not too small as our icons are 16x16
		// 19 is visual magic
		table.setRowHeight(Math.max(19, table.getRowHeight()));
		// Setup some column properties to make table look good.
		final int ColSizes[]    = {220, 140, 75, 60};
		TableColumnModel tcm = table.getColumnModel();
		for(int i = 0; (i < ColSizes.length); i++) {
			tcm.getColumn(i).setPreferredWidth(ColSizes[i]);
		}
		// Setup listener to validate reg exp after user modifies it.
		tcm.getColumn(1).setCellEditor(new DefaultCellEditor(new JTextField()));
		tcm.getColumn(1).getCellEditor().addCellEditorListener(this);
		// Create the color renderer for the table.
		ColorColumnRenderer = new ColorRenderer(table);
		ColorColumnEditor   = new ColorEditor();
		tcm.getColumn(2).setCellEditor(ColorColumnEditor);
		tcm.getColumn(2).setCellRenderer(ColorColumnRenderer);
		// Place the log in a scroll pane so that servers can be scrolled
		JScrollPane scroller = new JScrollPane(table);
		scroller.setBorder(null);
		scroller.getViewport().setBackground(table.getBackground());
		add(scroller, BorderLayout.CENTER);
	}
	
	/**
	 * Helper method to copy data from workspace into table model.
	 * 
	 * This is a helper method that is used to copy the data from
	 * the in-memory work space data structures into the table model
	 * to be displayed in the JTable.
	 */
	private void setTableData() {
		// Populate the data values.
		ArrayList<DBClassifier> classifiers = Workspace.get().getClassifierList().getClassifiers();
		for(DBClassifier dbClas: classifiers) {
			Object[] rowData = {dbClas.getDescription(), dbClas.getRegExp(),
					new Integer(dbClas.getColor()), 
					dbClas.isEnabled() ? Boolean.TRUE : Boolean.FALSE};
			tableData.addRow(rowData);
		}
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String ValidCmds = "add    delete up     down   check  revert " +
		"help   ok     apply  cancel";
		int cmdCode = ValidCmds.indexOf(e.getActionCommand());
		if (cmdCode == -1) {
			// Invalid command!
			return;
		}
		cmdCode /= 7; // normalize to make life easier
		final int selRow = table.getSelectedRow();
		switch(cmdCode) {
		case 0: // add
			tableData.addRow(new Object[]{"Genbank", "[a-zA-Z]+", 
					new Integer(0), Boolean.TRUE});
			// Select new row if another row is not selected.
			if (table.getSelectedRow() == -1) {
				int newRow = tableData.getRowCount() - 1;
				table.addRowSelectionInterval(newRow, newRow);
			}
			break;
		case 1: //delete
			tableData.removeRow(selRow);
			break;
		case 2: // up
			tableData.moveRow(selRow, selRow, selRow - 1);
			break;
		case 3: // down
			tableData.moveRow(selRow, selRow, selRow + 1);
			break;
		case 4: // test regular expression
			testRegExp();
			break;
		case 5: // revert
			setTableData();
			table.setRowSelectionInterval(0, 0);
			break;
		case 6: //help
			mainFrame.showHelp("http://www.peace-tools.org/");
			break;
		case 7: // ok
		case 8: // apply
			applyClassifiers(cmdCode == 7);
			break;
		case 9: // cancel
			checkAndClose();
			break;
		default: // command not yet implemented
		}
	}

	/**
	 * Method to check the working status of a DB classifier. 
	 * This method is invoked whenever the user clicks on the tool bar
	 * button to validate a classifier entry. This method first prompts
	 * the user to enter the data string to be matched with the regular
	 * expression. It then attempts to match the data string and provides
	 * the result back to the user.
	 */
	private void testRegExp() {
		int row = table.getSelectedRow();
		String regExp  = tableData.getValueAt(row, 1).toString();
		String message = String.format(CHECK_MSG, regExp);
		String data = JOptionPane.showInputDialog(this, message, 
				"Validating Regular Expression", JOptionPane.OK_CANCEL_OPTION);
		// Validate the regExp
		try {
			boolean result = data.matches(regExp);
			String resultMsg = String.format(CHECK_RESULT_MSG, 
					(result ? "" : "<b>not</b> "), regExp, data);
			JOptionPane.showMessageDialog(this, resultMsg, 
					"Results", (result ? JOptionPane.INFORMATION_MESSAGE : 
								JOptionPane.WARNING_MESSAGE));
		} catch (PatternSyntaxException pse) {
			// The regular expression is invalid. Report error and
			// force user to edit the entry again.
            JPanel msg = Utilities.collapsedMessage(REG_EX_ERR_MSG, Utilities.toString(pse));
            JOptionPane.showMessageDialog(this, msg, 
            	"Invalid Regular Expression", JOptionPane.ERROR_MESSAGE);
            // Force re-editing of this entry.
            table.editCellAt(table.getSelectedRow(), 1);
            // Don't proceed further
            return;
		}
	}
	
	/**
	 * Helper method to check and close the dialog.
	 * This method is a helper method that is called from a couple of
	 * different places to cross check with the user before abandoning
	 * edits and exiting.
	 */
	private void checkAndClose() {
		int choice = 
			JOptionPane.showConfirmDialog(this, CANCEL_MSG, 
					"Are you sure?", JOptionPane.YES_NO_OPTION);
		if (choice == JOptionPane.YES_OPTION) {
			setVisible(false);
			dispose();
		}
	}
	
	/**
	 * Helper method to validate regular expressions and apply them
	 * to the work space.
	 * 
	 * @param closeDialog If this flag is true, then this method also
	 * hides and disposes the dialog.
	 */
	private void applyClassifiers(boolean closeDialog) {
		// Validate all the entries in the classifiers.
		for(int i = 0; (i < tableData.getRowCount()); i++) {
			String regExp = tableData.getValueAt(i, 1).toString();
			// Validate the regExp
			try {
				// Compile regex to ensure it is valid.
				Pattern.compile(regExp);
			} catch (PatternSyntaxException pse) {
				// The regular expression is invalid. Report error and
				// force user to edit the entry again.
	            JPanel msg = Utilities.collapsedMessage(REG_EX_ERR_MSG, Utilities.toString(pse));
	            JOptionPane.showMessageDialog(this, msg, 
	            	"Invalid Regular Expression", JOptionPane.ERROR_MESSAGE);
	            // Force re-editing of this entry.
	            table.editCellAt(table.getSelectedRow(), 1);
	            // Don't proceed further
	            return;
			}
		}
		// All the regular expressions are valid. Save them to the 
		// work space. But first create a list of classifiers.
		ArrayList<DBClassifier> classifiers = new ArrayList<DBClassifier>();
		for(int i = 0; (i < tableData.getRowCount()); i++) {
			int color = ((Integer) tableData.getValueAt(i, 2)).intValue();
			boolean enabled = ((Boolean) tableData.getValueAt(i, 3)).booleanValue();
			DBClassifier cls = new DBClassifier(tableData.getValueAt(i, 1).toString(),
					tableData.getValueAt(i, 0).toString(), color, enabled); 
			classifiers.add(cls);
		}
		// Set classifiers in work space.
		ClassifierList clsList = Workspace.get().getClassifierList();
		clsList.set(classifiers);
		// Close the dialog as needed.
		if (closeDialog) {
			setVisible(false);
			dispose();
		}
	}
	
	/**
	 * Method intercepts table row selection notifications to enable/disable
	 * tool bar buttons.
	 * 
	 */
	@Override
	public void valueChanged(ListSelectionEvent event) {
		// Enable/disable various toolbar buttons.
		int index = table.getSelectedRow();
		// Delete button (enable only if selection is valid).
		toolbar.getComponent(1).setEnabled((index >=0) && (index < tableData.getRowCount()));
		// Up button (enable only if selection is not the first one)
		toolbar.getComponent(3).setEnabled((index > 0) && (index < tableData.getRowCount()));
		// Down button (enable only if selection is not the last one)
		toolbar.getComponent(4).setEnabled((index >= 0) && (index < tableData.getRowCount() - 1));
		// Validate button (enable only if selection is valid)
		toolbar.getComponent(6).setEnabled((index >=0) && (index < tableData.getRowCount()));
	}
	
	@Override
	public void editingCanceled(ChangeEvent event) {
		editingStopped(event);
	}

	@Override
	public void editingStopped(ChangeEvent event) {
		// Validate the regular expression entered by the user. 
		// Recollect that we set a default cell editor for this column.
		// Therefore the type casts below are safe.
		DefaultCellEditor dce = (DefaultCellEditor) event.getSource();
		String regExp = (String) dce.getCellEditorValue().toString();
		// Validate the regExp
		try {
			// Compile regex to ensure it is valid.
			Pattern.compile(regExp);
		} catch (PatternSyntaxException pse) {
			// The regular expression is invalid. Report error and
			// force user to edit the entry again.
            JPanel msg = Utilities.collapsedMessage(REG_EX_ERR_MSG, Utilities.toString(pse));
            JOptionPane.showMessageDialog(this, msg, 
            	"Invalid Regular Expression", JOptionPane.ERROR_MESSAGE);
            // Force re-editing of this entry.
            table.editCellAt(table.getSelectedRow(), 1);
		}
	}
	
	/**
	 * A simple renderer that uses a JLabel to render the color for a
	 * DB classifier.
	 * 
	 * This is an internal class that is used by this editor/view
	 * to suitably display the color associated with a DB classifier
	 * entry as a icon in this list. This renderer is set for the 
	 * third column.
	 */
	private class ColorRenderer extends DefaultTableCellRenderer {
		/**
		 * The constructor to setup initial configuration for rendering colors.
		 * 
		 * The constructor creates a few simple borders that are reused to 
		 * render colors in the column.
		 * 
		 * @param table The table which which this renderer is associated.
		 * The table is used to obtain suitable colors to ensure that the
		 * selections don't clash.
		 */
		ColorRenderer(JTable table) {
			setOpaque(true);
			selectedBorder = BorderFactory.createMatteBorder(2, 5, 2, 5,
					table.getSelectionBackground());
			unselectedBorder = BorderFactory.createMatteBorder(2, 5, 2, 5,
					table.getBackground());
		}
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			// Set the background color using value.
			Integer colorCode = (Integer) value;
			setBackground(new Color(colorCode.intValue()));
			setBorder(isSelected ? selectedBorder : unselectedBorder);
			return this;
		}

		/**
		 * The border to be used to draw color boxes when the row is
		 * not selected in the table. This border has a simple matte
		 * border drawn in the tables' selection background color.
		 */
		final Border unselectedBorder;

		/**
		 * The border to be used to draw color boxes when the row is
		 * selected in the table. This border has a simple matte
		 * border drawn in the tables' default background color.
		 */
		final Border selectedBorder;

		/**
		 * The general serial version UID to enable serialization of this class as
		 * per Java requirements.
		 */
		private static final long serialVersionUID = 5654108539980884223L;
	}


	/**
	 * Custom cell editor to display color chooser dialog to change colors.
	 * 
	 * This class provides a custom implementation for the color column 
	 * displayed by this classifier editor. This class permits the user
	 * to change the color via a color chooser dialog. This class provides
	 * a slightly non-traditional implementation for a cell editor in that
	 * it popups up another dialog for editing cell contents (rather than
	 * permitting the user to modify the contents directly inline in the
	 * table's cell).
	 */
	private class ColorEditor extends AbstractCellEditor implements TableCellEditor {
		/**
		 * Custom implementation for obtaining cell editor component.
		 * 
		 * This method implements the method in AbstractCellEditor interface.
		 * This method is supposed to return the editor component to be displayed
		 * inline in the table's cell to edit the color. However, this method
		 * implements this feature is a different way. Rather than returning an
		 * a component, this method pops up a color chooser dialog that permits
		 * a user to choose the color. Once the color is chosen this method updates
		 * the value in the column.
		 * 
		 * @param table The JTable that is asking the editor to edit.
		 * @param value The value of the cell to be edited. This method assumes that
		 * this value is an Integer object that represents the RGB color values.
		 * @param isSelected This value is true if the cell is to be rendered with highlighting. 
		 * This parameter is currently unused.
		 * @param row The row of the cell being edited.
		 * @param column The column of the cell being edited.
		 * @return This method always returns null, pretending like there is
		 * no editor to be displayed. 
		 */
		@Override
		public Component getTableCellEditorComponent(JTable table,
				Object value, boolean isSelected, int row, int column) {
			// Launch the color chooser dialog.
			color = new Color((Integer) value);
			Color changedColor = JColorChooser.showDialog(table, "Choose color for classifier", color);
			if (changedColor != null) {
				// The user actually changed the color. Save info.
				color = changedColor;
				table.setValueAt(new Integer(color.getRGB()), row, column);
			}
			// Return null pretending there is no editor.
			return null;
		}

		/**
		 * Implementation for method in CellEditor interface to return the
		 * color value set by the user.
		 * 
		 * This method returns the color value set by the user for the 
		 * given cell.
		 * 
		 * @return The color value set by the user for this cell.
		 */
		@Override
		public Object getCellEditorValue() {
			return color;
		}

		/**
		 * Override default method in parent class to correctly handle
		 * mouse button clicks.
		 *
		 * The default base class always returns true and causes the 
		 * dialog to popup whenever the user clicks on a cell. This method
		 * overrides the default implementation to permit the cell to be
		 * editable only after a double click.
		 * 
		 * @return This method returns true only if the user double clicks on
		 * the column.
		 */
		@Override
		public boolean isCellEditable(EventObject event) {
			if (event instanceof MouseEvent) {
				MouseEvent me = (MouseEvent) event;
				return (me.getClickCount() >= 2);
			}
			return false;
		}

		/**
		 * The color value that was last set by the user. This instance 
		 * variable is used to pass information between two API methods
		 * getTableCellEditorComponent() and  getCellEditorValue()
		 */
		private Color color;
		
		/**
		 * The generated serialization UID to keep compiler happy.
		 */
		private static final long serialVersionUID = 6443184579706912472L;
	}

	/**
	 * An instance of the color renderer that is used to render the
	 * third column in the table with a suitable colored label.
	 */
	private ColorRenderer ColorColumnRenderer;

	/**
	 * An instance of the color renderer that is used to render the
	 * third column in the table with a suitable colored label.
	 */
	private ColorEditor ColorColumnEditor;

	/**
	 * The JTable with the list of DB classifiers that are being 
	 * currently edited by the user.
	 */
	private JTable table;

	/**
	 * This instance variable contains the mutable data model that 
	 * contains the actual classifier data displayed in the table.
	 */
	private DefaultTableModel tableData;

	/**
	 * The tool bar that contains the various buttons/tools with
	 * with the classifier list can be edited and modified.
	 */
	private final JToolBar toolbar;

	/**
	 * Convenient reference to the main frame class that contains this
	 * dialog box. This value is set in the constructor
	 * and is never changed. The main frame reference is primarily used
	 * to handle help requests.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * A simple error message to be displayed to the user whenever
	 * a invalid regular expression is entered.
	 */
	private String REG_EX_ERR_MSG = 
		  "<html>Invalid regular expression syntax encountered. You must change<br>" +
		  "it to a valid regular expression before applying the changes.</html>";

	/**
	 * A conventional message to be displayed to the user prior to
	 * canceling out of this dialog.
	 */
	private String CANCEL_MSG = 
		  "<html>Any unapplied edits you have made will be lost.<br>" +
		  "Are you sure you want to cancel out of this dialog?</html>";

	/**
	 * A String.format() style message to be displayed to the 
	 * user while prompting to enter data string for validating a given
	 * regular expression. The actual regular expression is substituted
	 * by the testRegExp() method.
	 */
	private String CHECK_MSG = 
		  "<html>Enter the data string to be used to verify if the selected<br>" +
		  "regular expression mathes the supplied data string. The data string<br>" +
		  "is a sample EST data base FASTA identifier for testing.<br><br>" +
		  "The regular expression is: %s</html>";

	/**
	 * A String.format() style message to be displayed to the 
	 * user after validating a regular expression against a data string
	 * entered by the user. The actual regular expression and data string
	 * is substituted by the testRegExp() method.
	 */
	private String CHECK_RESULT_MSG = 
		  "<html><i>The supplied data string does %smatch the selected classifier entry.</i><br><br>" +
		  "The regular expression is: %s<br>" +
		  "The data string is: %s</html>";
	
	/**
	 * Generated serialization UID to keep compiler happy.
	 */
	private static final long serialVersionUID = -6093428974800271396L;
}
