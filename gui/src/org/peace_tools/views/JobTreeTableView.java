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
import java.awt.Dimension;
import java.awt.Point;
import java.awt.datatransfer.Transferable;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.activation.ActivationDataFlavor;
import javax.activation.DataHandler;
import javax.swing.Box;
import javax.swing.DropMode;
import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTable.DropLocation;
import javax.swing.JToolBar;
import javax.swing.ListSelectionModel;
import javax.swing.TransferHandler;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.swing.tree.TreePath;

import org.netbeans.swing.outline.DefaultOutlineModel;
import org.netbeans.swing.outline.Outline;
import org.netbeans.swing.outline.OutlineModel;
import org.netbeans.swing.outline.RenderDataProvider;
import org.peace_tools.core.AbstractMenuHelper;
import org.peace_tools.core.MainFrame;
import org.peace_tools.data.JobsTreeTableModel;
import org.peace_tools.generic.Log.LogLevel;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.JobList.JobOrListWrapper;

/**
 * This class provides a tabular view of the list of jobs that are
 * currently configured in this work space. This table uses the 
 *  JobListTableModel class that provides the Job data from the
 *  work space in a form that is easily displayed in a table.
 */
public class JobTreeTableView extends JPanel {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor sets up the job list table and 
	 * configures the table to the default configuration.
	 * 
	 * @param mainFrame The main frame that logically owns this job list
	 * view. The main frame is primarily used as the job listener which
	 * receives notifications on job completion. 
	 */
	public JobTreeTableView(MainFrame mainFrame, JobList root) {
		super(new BorderLayout(0, 0));
		// Set reference to main frame
		this.mainFrame = mainFrame;
		// First create the model and the table.
		model = new JobsTreeTableModel(root);
		// Create the actual GUI component for tree-table view
		OutlineModel outModel = 
			DefaultOutlineModel.createOutlineModel(model, model);
		jobTreeTable = new Outline(outModel) {
			private static final long serialVersionUID = -1248711551439232568L;
			@Override
			public TableCellRenderer getCellRenderer(int row, int column) {
				if (column == 1) {
					return StatusRenderer;
				}
				return super.getCellRenderer(row, column);
			}
		};
		// Ensure rows are not too small as our icons are 16x16
		// 19 is visual magic
		jobTreeTable.setRowHeight(Math.max(19, jobTreeTable.getRowHeight()));
        // Setup some column properties.
        TableColumnModel tcm = jobTreeTable.getColumnModel();
        tcm.getColumn(0).setHeaderValue("Testing");
		tcm.getColumn(0).setPreferredWidth(100);
		tcm.getColumn(1).setPreferredWidth(100);
		tcm.getColumn(2).setPreferredWidth(75);
		tcm.getColumn(3).setPreferredWidth(100);
		tcm.getColumn(4).setPreferredWidth(20);
        // Set some table properties
		jobTreeTable.setRenderDataProvider(new NameStatusRenderer(jobTreeTable));
		jobTreeTable.setBorder(null);
		jobTreeTable.setShowHorizontalLines(true);
		jobTreeTable.setGridColor(new Color(0xe0, 0xe0, 0xe0));
		jobTreeTable.setRootVisible(false);
		jobTreeTable.setMinimumSize(new Dimension(50, 50));
		// Setup DnD operations on the job tree
		jobTreeTable.setDragEnabled(true);
		jobTreeTable.setTransferHandler(new JobTransferHandler(jobTreeTable));
		jobTreeTable.setDropMode(DropMode.ON);
		jobTreeTable.setColumnSelectionAllowed(false);
        // Place the tree-table in a scroll pane so that jobs can be scrolled
        JScrollPane scroller = new JScrollPane(jobTreeTable);
        scroller.setBorder(null);
        scroller.getViewport().setBackground(jobTreeTable.getBackground());
        add(scroller, BorderLayout.CENTER);
        // Create the tool-bar at the top of the job list
        toolbar = new JToolBar();
        toolbar.setFloatable(false);
        AbstractMenuHelper helper = mainFrame.getMenuHelper(AbstractMenuHelper.HelperType.JOB_MENU);
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.START_JOB_MONITOR, false));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.STOP_JOB_MONITOR, false));
        toolbar.add(Box.createHorizontalStrut(5));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.SHOW_JOB_DETAILS, false));
        toolbar.add(Box.createHorizontalStrut(5));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.REMOVE_JOB, false));
        toolbar.add(Box.createHorizontalStrut(5));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.CREATE_SUBLIST, false));        
        // Add tool bar to the north
        add(toolbar, BorderLayout.NORTH);
        // Finally create pop up menu for various job options
        createPopupMenu();
        // Add a mouse handler to trigger pop-up menu
        addMouseAdapter(jobTreeTable);
		// Set the selection model for this table.
        jobTreeTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        jobTreeTable.setRowSelectionAllowed(true);
        jobTreeTable.getSelectionModel().addListSelectionListener(helper.getListSelectionListener(jobTreeTable));
		// Also setup the helper as a model listener as well.
		model.addTableModelListener(helper);
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
		// Add various menu options to the main menu. 
		AbstractMenuHelper helper = mainFrame.getMenuHelper(AbstractMenuHelper.HelperType.JOB_MENU);
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.COMPUTE_MST, false));
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.COMPUTE_CLUSTERS, false));
		popupMenu.addSeparator();
		// Context sensitive job operation menu items.
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.START_JOB_MONITOR, false));
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.STOP_JOB_MONITOR, false));
		popupMenu.addSeparator();
		
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.SHOW_JOB_DETAILS, false));
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.ABORT_JOB, false));
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.REMOVE_JOB, false));
		popupMenu.addSeparator();
		
		// Create menu options to add/delete job sublists
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.CREATE_SUBLIST, false));
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.DELETE_SUBLIST, false));
		popupMenu.addSeparator();
		
		// Finally create the server entries in the menu.
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.SHOW_JOBS_ON_SERVER, false));
		popupMenu.add(helper.getMenuItem(AbstractMenuHelper.ActionType.SHOW_MY_JOBS_ON_SERVER, false));
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
		final Point mousePos = me.getPoint();
		int row = jobTreeTable.rowAtPoint(mousePos);
		jobTreeTable.clearSelection();
		if (jobTreeTable.getClosestPathForLocation(mousePos.x, mousePos.y) != null) {
	        // Select the table entry
			jobTreeTable.setRowSelectionInterval(row, row);
		}
        // JobMenuHelper enables/disables pop-up menu items based 
		// on the item selected.
        popupMenu.show(jobTreeTable, me.getX(), me.getY());
	}

	/**
	 * Helper method to handle double click of the mouse on a list item.
	 * 
	 * This method is invoked whenever the user double clicks on a
	 * row in the job list. This method checks to see if the
	 * entry is valid and if so, opens a tab with information about 
	 * the job running on the server.
	 * 
	 * @param me The mouse event associated with the double click.
	 */
	private void handleDoubleClick(MouseEvent me) {
		assert (me.getClickCount() == 2);
		// Obtain the row that was selected 
		int row = jobTreeTable.rowAtPoint(me.getPoint());
		// Get the job at the given row
		JobOrListWrapper wrapper = 
			(JobOrListWrapper) jobTreeTable.getOutlineModel().getValueAt(row, -1);
		if (wrapper == null) {
			return;
		}
	}
	
	/**
	 * A simple tree cell renderer that that essentially provides a better 
	 * representative set of Icons to make the overall display 
	 * look a bit prettier.
	 * 
	 * This is an internal class that is used by this view
	 * to suitably display the job ID along with a status
	 * icon. This renderer is set for the first column in the Job tree
	 * table view.
	 */
	private class NameStatusRenderer implements RenderDataProvider {
		/**
		 * The tree table model that is using this renderer. This reference is
		 * used to detect if a given object is a drop target to decide if it
		 * get's a special background color to highlight that it is a candidate
		 * drop target.
		 */
		private final Outline treeTable;
		
		/**
		 * The default constructor. Currently it does nothing
		 * and is here merely to adhere to coding conventions.
		 * 
		 * @param treeTable The Outline object that is logically using this
		 * renderer.
		 */
		public NameStatusRenderer(Outline treeTable) {
			this.treeTable = treeTable;
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
			DropLocation dropLoc = treeTable.getDropLocation();
			if ((dropLoc != null) && 
				(treeTable.getValueAt(dropLoc.getRow(), dropLoc.getColumn()) == o)) {
				return treeTable.getSelectionBackground();
			}
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
			JobOrListWrapper wrapper = (JobOrListWrapper) o;
			if (wrapper.isSubList()) {
				JobList jl = wrapper.getSubList();
				return "<html><i>" + jl.getName() + "</i></html>";
			}
			// This is a job entry
			final Job job = wrapper.getJob();
			return "<html>" + job.getJobID() + "</html>";
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
			// Get status for this object depending on whether it is a job or a sublist
			JobOrListWrapper wrapper = (JobOrListWrapper) o;
			if (wrapper.isSubList()) {
				return JobStatusIcons[8];
			}
						
			JobBase.JobStatusType status = wrapper.getJob().getStatus();
			return JobStatusIcons[status.ordinal()];
		}

		/**
		 * Provide an tool tip for the node. Provides a default
		 * implementation for the RenderDataProvider API.
		 * 
		 * @return This method currently returns null for all 
		 * objects.
		 */
		public String getTooltipText(Object o) {
			JobOrListWrapper wrapper = (JobOrListWrapper) o;
			if (wrapper.isSubList()) {
				return "<html>" + wrapper.getSubList().getDescription() + "</html>";
			} else {
				return "<html>" + wrapper.getJob().getDescription() + "</html>";
			}
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
			if (!(value instanceof JobBase.JobStatusType)) {
				return defaultRenderer.getTableCellRendererComponent(table, value, 
						isSelected, hasFocus, row, column);
			}
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
			JobOrListWrapper wrapper = (JobOrListWrapper) table.getValueAt(row, -1);
			if ((wrapper != null) && !wrapper.isSubList()) {
				progress = wrapper.getJob().getProgressInfo();
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
		Utilities.getIcon("images/16x16/JobWaiting.png"),
		Utilities.getIcon("images/16x16/JobQueued.png"),
		Utilities.getIcon("images/16x16/JobRunning.png"),
		Utilities.getIcon("images/16x16/JobFinishing.png"),
		Utilities.getIcon("images/16x16/JobSuccess.png"),
		Utilities.getIcon("images/16x16/JobError.png"),
		Utilities.getIcon("images/16x16/JobError.png"),
		Utilities.getIcon("images/16x16/JobSet.png")
	};
	
	/**
	 * An inner class that facilitates drag-n-drop (DnD) of jobs and
	 * job lists.
	 * 
	 * This class extends Java's DnD infrastructure class and implements
	 * selected method that are needed to accomplish DnD to help
	 * users organize hierarchical JobList and Job entries.
	 * 
	 */
	private class JobTransferHandler extends TransferHandler {
		/**
		 * A serialization GUID (just to keep compiler happy) 
		 */
		private static final long serialVersionUID = -3845109707292047030L;
		
		/**
		 * The top-level tree-table GUI component to which this 
		 * transfer handler has been associated. This value is set
		 * in the constructor.
		 */
		private final Outline treeTable;
		
		/**
		 * A special DnD data flavor with associated mime-type to be
		 * used to uniquely detect the DnD objects created by this
		 * class. 
		 */
		private final ActivationDataFlavor JobTreeTableFlavor =
			new ActivationDataFlavor(TreePath.class, "application/jobTreeTable", "JobTreeTableEntry");
		
		/**
		 * The constructor.
		 * 
		 * The constructor merely sets up the references to the GUI
		 * components for future references.
		 * 
		 * @param treeTable The tree-table GUI components with which
		 * this transfer handler is logically associated with. This
		 * parameter cannot be null.
		 */
		public JobTransferHandler(Outline treeTable) {
			this.treeTable = treeTable;
		}
		
		/**
		 * Override base class method to add a node at the appropriate
		 * location in the job list hierarchy.
		 * 
		 * This method is invoked by the underlying Java components
		 * when the user drags an entry (may it be a Job or JobList) 
		 * and drops it on the appropriate location. This method
		 * uses the drop location to determine the position in
		 * the hierarchical job list to insert the entry.
		 * 
		 * @param dndInfo The Java DnD infrastructure object that
		 * contains the information about the object being dragged
		 * and the location where it is to be dropped.
		 */
		@Override
		public boolean importData(TransferSupport dndInfo) {
			// Extract the job node to be inserted from the drag-n-drop object
			final TreePath srcPath = getPath(dndInfo.getTransferable());
			// Obtain the final drop location where the hierarchy is to be inserted
			final DropLocation loc = dndInfo.getDropLocation();
			// Validate to ensure things look good for further processing
			if ((srcPath == null) || (loc == null) || 
				!(loc instanceof JTable.DropLocation)) {
				// The drop location is not on a JTable like we expect.
				// That does not work for this target.
				return false;
			}
			// Obtain the closest path to the drop location to determine
			// the job list into which the hierarchy is to be inserted.
			JTable.DropLocation tableLoc = (JTable.DropLocation) loc;
			final TreePath targetPath = treeTable.getLayoutCache().getPathForRow(tableLoc.getRow());
			if (targetPath == null) {
				// Can't figure out target path. Nothing can be done.
				return false;
			}
			// Now locate the appropriate target job list and target position
			// (within the target job list) where the element is to be inserted.
			int targetPos = -1;
			JobList targetList = null;
			final JobOrListWrapper targetLeaf = (JobOrListWrapper) 
				targetPath.getLastPathComponent();
			if (!targetLeaf.isSubList()) {
				// The target is a job. in this case we need to find the 
				// immediately enclosing list node.
				JobOrListWrapper jobListNode = (JobOrListWrapper) 
					targetPath.getPathComponent(targetPath.getPathCount() - 2);
				assert(jobListNode.isSubList());
				targetPos  = jobListNode.getSubList().getJobs().indexOf(targetLeaf);
				targetList = jobListNode.getSubList();
			} else {
				targetList = targetLeaf.getSubList();
				targetPos  = targetList.getJobs().size();
			}
			// Insert the target node at the target position. The insertions
			// will fire suitable events and updating the GUI.
			JobOrListWrapper srcWrapper = (JobOrListWrapper) srcPath.getLastPathComponent();
			if (srcWrapper.isSubList()) {
				targetList.insert(srcWrapper.getSubList(), targetPos);
			} else {
				targetList.insert(srcWrapper.getJob(), targetPos);
			}
			// Set the newly inserted node as the selected item.
			treeTable.changeSelection(targetLeaf.isSubList() ? (tableLoc.getRow() + targetPos + 1) : tableLoc.getRow(), 
					0, false, false);
			// If the node was previously expanded ensure the new node remains expanded.
			if (treeTable.isExpanded(srcPath)) {
				treeTable.expandPath(targetPath);
			}
			return true;
		}

		/**
		 * Helper method that is used to obtain the tree path from the
		 * current transferable data. 
		 * 
		 * @param dndData The drag-n-drop object from which the transfer
		 * hierarchy and data is to be obtained.
		 * 
		 * @return The tree path that indicates the entry being
		 * dragged.
		 */
		private TreePath getPath(Transferable dndData) {
			TreePath hierarchy = null;
			try {
				Object tmpTreePath = dndData.getTransferData(JobTreeTableFlavor);
				if ((tmpTreePath != null) && (tmpTreePath instanceof TreePath)) {
					// Obtain the job hierarchy being moved. 
					hierarchy = (TreePath) tmpTreePath;
					// Ensure we have at least two levels as the root cannot be really moved.
					final int pathSize = hierarchy.getPathCount();
					if (pathSize < 2) {
						UserLog.log(LogLevel.WARNING, "JobTreeTableView", 
							"The path being moved has less than 2 levels. This is not expected!");
					}
				}
			} catch (Exception e) {
				ProgrammerLog.log(e);
			}
			return hierarchy;	
		}

		/**
		 * Implementation for Java's DnD method to indicate if currently
		 * cursor location is a valid location to drop.
		 * 
		 * This method is invoked by the underlying DnD methods to
		 * determine if the current mouse position is a valid
		 * location to drop the current item being dragged.
		 * This method returns true only when:
		 * 
		 * <ol>
		 * <li>The current operation is to move the entry (and not
		 * to copy the entry).</li>
		 * 
		 * <li>If the object being dragged is of type JobTreeTableFlavor.</li>
		 * 
		 * <li>If the drop location corresponds to the zero-th or Job
		 * column.</li>
		 * </ol>
		 * 
		 * @return This method returns true if the current mouse position
		 * is a valid drop location for the current entry being dragged.
		 */
		@Override
		public boolean canImport(TransferSupport dndInfo) {
			if (!dndInfo.isDrop() || (dndInfo.getDropAction() != MOVE)) {
				// Currently we only support move operations.
				return false;
			}
			dndInfo.setShowDropLocation(true);
			// Check if the data is same as the flavor this class creates.
			if (!dndInfo.isDataFlavorSupported(JobTreeTableFlavor)) {
				return false;
			}
			if (getPath(dndInfo.getTransferable()) == null) {
				// The object being dragged is not the kind we are looking for.
				return false;
			}
			// The job tree table currently accepts drops only on the primary node column
			DropLocation loc = dndInfo.getDropLocation();
			if (!(loc instanceof JTable.DropLocation)) {
				// The drop location is not on a JTable like we expect.
				// That does not work for this target.
				return false;
			}
			JTable.DropLocation tableLoc = (JTable.DropLocation) loc;
			// Return true only if the object being dropped on column zero
			return (tableLoc.getColumn() == jobTreeTable.convertColumnIndexToView(0));
		}

		/**
		 * Return the DnD actions supported by this method.
		 * 
		 * This method overrides the default implementation in the base
		 * class to report that only Move operations are supported by
		 * this transfer handler.
		 * 
		 * @return This method returns only {@link TransferHandler#MOVE}
		 * as the acceptable operation. 
		 */
		@Override
		public int getSourceActions(JComponent c) {
			return MOVE;
		}

		/**
		 * Creates the DnD object containing the information to be transfered.
		 * 
		 * This method is invoked by the underlying DnD classes when
		 * the user attempts to start DnD operations. This method uses
		 * the currently selected entry in the GUI component as the
		 * component being dragged. It creates a suitable 
		 * DataHandler component with the JobTreeTableFlavor to be
		 * used for DnD operations.
		 * 
		 * @param c The component on which the DnD operation is being
		 * performed. Currently this parameter is not really used
		 * by this method.
		 * 
		 * @return The transferable object that contains information
		 * about the node in the job hierarchy being dragged.
		 */
		@Override
		protected Transferable createTransferable(JComponent c) {
			final int selRow    = treeTable.getSelectedRow();
			final TreePath path = treeTable.getOutlineModel().getLayout().getPathForRow(selRow);
			if (path != null) {
				// System.out.printf("The path %s %s expanded\n", path.toString(),
				//		(treeTable.isExpanded(path) ? " is " : " is NOT "));
				return new DataHandler(path, JobTreeTableFlavor.getMimeType());
			}
			// Maybe the base class can handle it...
			return super.createTransferable(c);
		}

		/**
		 * Method to remove the node that has been successfully dragged.
		 * 
		 * This method is invoked by the underlying Java DnD infrastructure
		 * after the user has successfully dropped a node on the final
		 * drop location. This method is invoked to remove the node
		 * from its original position in the hierarchy.
		 * 
		 * @param source The source GUI component that was involved in 
		 * the DnD operation. Currently, this parameter is not used and the
		 * {@link #treeTable} is used.
		 * 
		 * @param data The component that encapsulates information about the 
		 * component being dragged.
		 */
		@Override
		protected void exportDone(JComponent source, Transferable data,
				int action) {
			// Extract the tree path component being dragged from the data.
			final TreePath hierarchy = getPath(data);
			if ((action != TransferHandler.MOVE) || (hierarchy == null)) {
				return;
			}
			// Now remove the node being moved from our model.
			final JobOrListWrapper wrapper = (JobOrListWrapper) 
				hierarchy.getPathComponent(hierarchy.getPathCount() - 2);
			if (wrapper.isSubList()) {
				// Remove the node. The underlying model will
				// automatically fire events to update the GUI. 
				wrapper.getSubList().remove((JobOrListWrapper) hierarchy.getLastPathComponent());
			}
		}
	}
	
	/**
	 * The model that we are using to render the information in the
	 * tree-table view of jobs.
	 */
	private final JobsTreeTableModel model;
	
	/**
	 * The tool-bar that contains some commonly used tools with 
	 * the jobs.
	 */
	private JToolBar toolbar;
	
	/**
	 * The actual Tree-table that provides a graphical view of the
	 * hierarchical structure of job lists and jobs they contain.
	 */
	private Outline jobTreeTable;

	/**
	 * The pop up menu that is displayed when the user left-clicks
	 * on an job entry in this view. 
	 */
	private JPopupMenu popupMenu;

	/**
	 * An instance of the cell renderer that is used to render the
	 * second column in the table to display status or progress
	 * information.
	 */
	private final JobStatusRenderer StatusRenderer = 
		new JobStatusRenderer();

	/**
	 * The main frame that logically owns this job list view. The 
	 * main frame is primarily used as the job listener which
	 * receives notifications on job completion.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * A generated serial version ID for serialization (more
	 * realistically to keep the compiler happy). 
	 */
	private static final long serialVersionUID = 80617431851108817L;
}
