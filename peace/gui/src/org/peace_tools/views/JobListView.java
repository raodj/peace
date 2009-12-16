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
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import org.peace_tools.data.JobListTableModel;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;

/**
 * This class provides a tabular view of the list of jobs that are
 * currently configured in this work space. This table uses the 
 *  JobListTableModel class that provides the Job data from the
 *  work space in a form that is easily displayed in a table.
 */
public class JobListView extends JPanel implements ActionListener {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor sets up the job list table and 
	 * configures the table to the default configuration. 
	 */
	public JobListView() {
		super(new BorderLayout(0, 0));
		// First create the model and the table.
		model = new JobListTableModel();
		jobTable = new JTable(model) {
			private static final long serialVersionUID = 1430052270991586572L;
			@Override
			public TableCellRenderer getCellRenderer(int row, int column) {
				if (column == 0) {
					return JobNameRenderer; 
				}
				if (column == 1) {
					return StatusRenderer;
				}
				return super.getCellRenderer(row, column);
			}
		};
		// Ensure rows are not too small as our icons are 16x16
		// 19 is visual magic
		jobTable.setRowHeight(Math.max(19, jobTable.getRowHeight()));
        // Setup some column properties.
        TableColumnModel tcm = jobTable.getColumnModel();
		tcm.getColumn(0).setPreferredWidth(50);
		tcm.getColumn(1).setPreferredWidth(100);
		tcm.getColumn(2).setPreferredWidth(50);
		tcm.getColumn(3).setPreferredWidth(50);
        // Set some table properties
		jobTable.setBorder(null);
        jobTable.setShowHorizontalLines(true);
        jobTable.setFillsViewportHeight(true);
        jobTable.setDragEnabled(false);
        jobTable.setDropTarget(null);
        jobTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        // Place the log in a scroll pane so that jobs can be scrolled
        JScrollPane scroller = new JScrollPane(jobTable);
        scroller.setBorder(null);
        scroller.getViewport().setBackground(jobTable.getBackground());
        add(scroller, BorderLayout.CENTER);
        // Create the toolbar at the top of the job list
        toolbar = new JToolBar();
        toolbar.setFloatable(false);
        toolbar.add(Utilities.createButton("images/16x16/Delete.png", 
        		null, "DeleteJob", null, 
        		"Delete the currently selected job entry", false));
        toolbar.add(Utilities.createButton("images/16x16/DeleteOldJobs.png", 
        		null, "DeleteOldJobs", null, 
        		"Delete all old job entries", true));
        // Add tool bar to the north
        add(toolbar, BorderLayout.NORTH);
        // Finally create pop up menu for various job options
        createPopupMenu();
        // Add a mouse handler to trigger popups
        addMouseAdapter(jobTable);
	}
	
	/**
	 * This is a helper method to create the pop-up menu.
	 * 
	 * This is a helper method that was introduced to streamline the
	 * code in the constructor. This method creates a popup menu
	 * that provides options for monitoring and controlling 
	 * jobs and outputs. 
	 */
	private void createPopupMenu() {
		popupMenu = new JPopupMenu();
		popupMenu.add(Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Start job monitor (daemon thread)", "startMonitor", this, 
				"images/16x16/StartJobMonitor.png", null, true, false));
		popupMenu.add(Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Stop job monitor (daemon thread)", "stopMonitor", this, 
				"images/16x16/StopJobMonitor.png", null, true, false));
		popupMenu.addSeparator();
		//--------------------------------------------------
		popupMenu.add(Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Show outputs from the job", "view", this, 
				"images/16x16/DefaultView.png", null, true, false));
		popupMenu.add(Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Show all jobs on server", "jobListAll", this, 
				"images/16x16/ServerInfo.png", null, true, false));
		popupMenu.add(Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Show my jobs on server", "jobList", this, 
				null, null, true, false));
		popupMenu.addSeparator();
		//--------------------------------------------------
		popupMenu.add(Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Remove entry from Workspace", "delete", this, 
				"images/16x16/Delete.png", null, true, false));

	}
	
	/**
	 * A refactored helper method to add a mouse adapter. This method
	 * adds a mouse adapter to intercept certain mouse events occurring
	 * on the data set tree to trigger various operations. The mouse
	 * adapter simply delegates the actual operations to other methods 
	 * in this class.
	 * 
	 * @param list The list object to which the mouse adapter is to be
	 * added.
	 */
	private void addMouseAdapter(JComponent list) {
		list.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				if (e.getClickCount() == 2) {
					handleDoubleClick(e);
				}
			}
			@Override
			public void mousePressed(MouseEvent e) {
				if (e.isPopupTrigger()) {
					handlePopup(e);
				}
			}
			@Override
			public void mouseReleased(MouseEvent e) {
				if (e.isPopupTrigger()) {
					handlePopup(e);
				}
			}
		});
	}
	
	/**
	 * Helper method to left mouse click on a table item.
	 * 
	 * This method is invoked whenever the user clicks on the left
	 * mouse button on a item in the job table. This method checks
	 * to see if the entry is valid and if so, pops up a menu with
	 * valid operations for the selected job entry.
	 * 
	 * @param me The mouse event associated with the mouse click.
	 */
	public void handlePopup(MouseEvent me) {
		int row = jobTable.rowAtPoint(me.getPoint());
		// Get the file at the given row and column.
		Job job = model.getJob(row);
		if (job == null) {
			return;
		}
        // Select the table entry
		jobTable.setRowSelectionInterval(row, row);
        // Now enable/disable popup menu items based on the item
        // selected in the list.
        // popupMenu.updateItems(entry, false);
        // Show pop-up menu.
        popupMenu.show(jobTable, me.getX(), me.getY());
	}
	
	/**
	 * Helper method to handle double click of the mouse on a list item.
	 * 
	 * This method is invoked whenever the user double clicks on a
	 * row in the  set tree. This method checks to see if the
	 * entry is valid and if so, opens a view of the specified
	 * data file. 
	 * 
	 * @param me The mouse event associated with the double click.
	 */
	private void handleDoubleClick(MouseEvent me) {
		assert (me.getClickCount() == 2);
		// Obtain the row that was selected 
		int row = jobTable.rowAtPoint(me.getPoint());
		// Get the file at the given row and column.
		Object entry = model.getJob(row);
		if (entry == null) {
			return;
		}
	}

	/**
	 * A simple renderer that uses a JLabel to render JobName and
	 * icon indicating job status.
	 * 
	 * This is an internal class that is used by the job list view
	 * to suitably display the name of the job along with a status
	 * icon. This renderer is set for the first column in the Job list
	 * view.
	 */
	private class JobNameStatusRenderer extends DefaultTableCellRenderer {
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			// Let the user class set up the renderer with defaults
			super.getTableCellRendererComponent(table, value, isSelected, hasFocus, 
					row, column);
			// Get status for this job from the 3rd column.
			JobBase.JobStatusType status = (JobBase.JobStatusType) model.getValueAt(row, 1);
			setIcon(JobStatusIcons[status.ordinal()]);
			return this;
		}
		/**
		 * The general serial version UID to enable serialization of this class as
		 * per Java requirements.
		 */
		private static final long serialVersionUID = 5654108539980884223L;
	}

	private class JobStatusRenderer implements TableCellRenderer {
		/**
		 * The label with a given status string that is used when the job
		 * is in a non-running status.
		 */
		private DefaultTableCellRenderer defaultRenderer;
		
		/**
		 * The progress bar is used when the job is in running status to
		 * indicate the amount of progress made thus far.
		 */
		private JProgressBar progressBar;
		
		/** Default constructor.
		 * 
		 * The default constructor merely sets the default font to be used
		 * by the label to be non-bold.
		 */
		public JobStatusRenderer() {
			// Create defaults.
			defaultRenderer = new DefaultTableCellRenderer();
			progressBar     = new JProgressBar();
			progressBar.setBorder(new CompoundBorder(
					new EmptyBorder(2, 5, 2, 5), new EtchedBorder()));
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			JobBase.JobStatusType status = (JobBase.JobStatusType) value;
			if (!JobBase.JobStatusType.RUNNING.equals(status)) {
				// In non running status, just use the default
				// cell renderer.
				return defaultRenderer.getTableCellRendererComponent(table, value, 
						isSelected, hasFocus, row, column);
			}
			// When control drops here that means status is running. In
			// this mode return the progress bar suitably updated.
			int progress[] = {-1, -1}; // #ESTs, #Total ESTs
			Job job = model.getJob(row);
			if (job != null) {
				progress = job.getProgressInfo();
			}
			// Set progress bar configuration based on progress made.
			progressBar.setBackground(isSelected ? table.getSelectionBackground() : table.getBackground());
			if ((progress[0] != -1) && (progress[1] != -1)) {
				progressBar.setMaximum(Math.max(1, progress[1]));
				progressBar.setValue(progress[0]);
			}
			// Use progress bar for rendering this column information
			return progressBar;
		}
		/**
		 * The general serial version UID to enable serialization of this class as
		 * per Java requirements.
		 */
		private static final long serialVersionUID = 5654108539980884223L;
	}
	
	/**
	 * The list of icons that are used by job name status table cell 
	 * renderer to provide visual cues about the current status of 
	 * a job. The entries are listed in the same order as they are
	 * enumerated in the JobBase.JobListType enumeration. The order
	 * is important to ensure that the icons match up with the status
	 * information correctly.
	 */
	private static final Icon JobStatusIcons[] = {
		Utilities.getIcon("images/16x16/JobStarting.png"),
		Utilities.getIcon("images/16x16/JobQueued.png"),
		Utilities.getIcon("images/16x16/JobRunning.png"),
		Utilities.getIcon("images/16x16/JobFinishing.png"),
		Utilities.getIcon("images/16x16/JobSuccess.png"),
		Utilities.getIcon("images/16x16/JobError.png")
	};
	
	/**
	 * The model that we are using to render the information in the
	 * tablular view of jobs.
	 */
	private final JobListTableModel model;
	
	/**
	 * The toolbar that contains some commonly used tools with 
	 * the jobs.
	 */
	private JToolBar toolbar;
	
	/**
	 * The actual JTable that provides a graphical view of the
	 * list of jobs currently configured on this workspace.
	 */
	private JTable jobTable;

	/**
	 * The pop up menu that is displayed when the user left-clicks
	 * on an item in the job list view. 
	 */
	private JPopupMenu popupMenu;
	
	/**
	 * An instance of the cell renderer that is used to render the
	 * first column in the table with a suitable status icon.
	 */
	private final JobNameStatusRenderer JobNameRenderer = 
		new JobNameStatusRenderer();

	/**
	 * An instance of the cell renderer that is used to render the
	 * second column in the table to display status or progress
	 * information.
	 */
	private final JobStatusRenderer StatusRenderer = 
		new JobStatusRenderer();

	/**
	 * A generated serial version ID for serialization (more
	 * realistically to keep the compiler happy). 
	 */
	private static final long serialVersionUID = 80617431851108817L;

	@Override
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub
		
	}
}
