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

package org.peace_tools.data;

import java.util.ArrayList;

import javax.swing.table.AbstractTableModel;

import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;
import org.peace_tools.workspace.WorkspaceEvent;
import org.peace_tools.workspace.WorkspaceListener;
import org.peace_tools.workspace.JobList.JobOrListWrapper;

/**
 * A bridge class between Job entries in a workspace and a JTable.
 * 
 * This class serves as a bridge between the in-memory representation of
 * jobs in a Workspace (represented by the set of classes in the 
 * <code>workspace</code> package). This class enables reusing the jobs
 * hierarchy maintained by the Workspace object to display it in a
 * JTable.
 * 
 * <p><b>Note:</b>  This table model currently provides the following information
 * for each job: JobID, Server, Monitor, CPUs, Status. The monitor column
 * indicates the status of the background job monitoring thread.</p>
 */
public class JobListTableModel extends AbstractTableModel 
implements WorkspaceListener {
	/**
	 * The top-level wrapper that contains the top-level job list that 
	 * this table model is meant to operate on. This value is set when 
	 * the view associated with this model is created.
	 */
	private final JobOrListWrapper root;

	/**
	 * This array maintains a flattened list of job entries.
	 * This list is dynamically built (on demand) whenever
	 * this model is used in a table and quick responses
	 * to individual row entries are required. Whenever the
	 * tree-model changes this array list is reset to null,
	 * forcing list to be rebuilt.
	 */
	private ArrayList<JobOrListWrapper> jobEntries = null; 

	/**
	 * The only constructor for this class.
	 * 
	 * The constructor merely initializes the instance variables
	 * to their default initial values.
	 * 
	 * @param root The top-level job list from where this model must
	 * obtain the necessary information.
	 */
	public JobListTableModel(final JobList root) {
		this.root = root.new JobOrListWrapper(root);
		Workspace.get().addWorkspaceListener(this);
	}

	/**
	 * Obtain the class that describes the data type of a given 
	 * column.
	 * 
	 * This method overrides the API method in RowModel.
	 * 
	 * @param column The zero-based index of the column whose
	 * Class type is to be returned.
	 */
	@Override
	public Class<?> getColumnClass(int column) {
		if (column == 1) {
			return JobBase.JobStatusType.class;
		}
		return String.class;
	}

	/**
	 * This method overrides the interface method in RowModel to 
	 * return the number of columns in the TreeTable.
	 * 
	 * This method overrides the API method in AbstractTableModel (and
	 * also satisfies the RowModel API).
	 *  
	 * @return The number of columns to be displayed in this table. There
	 * are always 5 columns to be displayed.
	 */
	@Override
	public int getColumnCount() {
		return 5;
	}

	/**
	 * Override the default names that are set for the columns 
	 * displayed in the job table.
	 * 
	 * @param col The zero-based column index whose title is to 
	 * be returned.
	 * 
	 * @return The title associated with the column.
	 */
	@Override
	public String getColumnName(int col) {
		String titles[] = {"Job ID", "Status", "Monitor", "Server", "CPUs"};
		return titles[col];
	}

	/**
	 * Interface method to determine if an entry in the model is 
	 * mutable.
	 * 
	 * @return This method always returns false to indicate that
	 * the data cannot be modified.
	 */
	@Override
	public boolean isCellEditable(int row, int column) {
		return false;
	}

	/**
	 * Helper method to recursively populate all the job entries in a
	 * given job list.
	 * 
	 * This method is used to add the job entries contained within 
	 * a given job list. This method recursively descends into a job
	 * list and adds job entries to the job list.
	 * 
	 * @param jobList The job list whose entries are to be recursively
	 * processed.
	 * 
	 * @param targetList The array list to which job wrapper's are to
	 * be added.
	 */
	private void addJobEntries(JobOrListWrapper jlw, 
			ArrayList<JobOrListWrapper> targetList) {
		if (jlw.isSubList()) {
			for(JobOrListWrapper wrapper: jlw.getSubList().getJobs()) {
				addJobEntries(wrapper, targetList);
			} 
		} else {
			targetList.add(jlw);
		}
	}

	/**
	 * Helper method to populate the {@link #jobEntries} list with jobs
	 * for rapid access to entries to display in a tabular view.
	 * 
	 * This method was primarily introduced to streamline building
	 * the {@link #jobEntries}. This essential aspect of this method
	 * is that is synchronized and consequently multi-threading safe.
	 * Note that
	 */
	private synchronized void buildJobEntries() {
		if (jobEntries == null) {
			ArrayList<JobOrListWrapper> jobEntryList = 
				new ArrayList<JobOrListWrapper>();
			addJobEntries(root, jobEntryList);
			jobEntries = jobEntryList;
		}
	}

	/**
	 * This method implements the abstract method defined
	 * in {@link AbstractTableModel#getRowCount()}.
	 * 
	 * This method ensures that the {@link #jobEntries} list is
	 * populated by invoking the {@link #buildJobEntries()}.
	 * 
	 * @return This method return the number of job entries in the
	 * model.
	 */
	@Override
	public int getRowCount() {
		buildJobEntries();
		return jobEntries.size();
	}

	/**
	 * Obtain the value to be displayed at a given row and column.
	 * 
	 * <p><b>NOTE</b>: There is a difference in the values returned
	 * by this method and the {@link #getValueAt(int, int)} method
	 * for different columns.</p> 
	 * 
	 * This method implements the abstract method defined in 
	 * {@link AbstractTableModel#getValueAt(int, int)} method.
	 * 
	 * This method ensures that the {@link #jobEntries} list is
	 * populated by invoking the {@link #buildJobEntries()}.
	 *
	 * @param row The row in the list of jobs for which a value
	 * is to be returned.
	 * 
	 * @param column The column in the jobs table for 
	 * @return The object to be displayed at a given row and column.
	 * The value depends on the column being displayed.
	 */
	@Override
	public Object getValueAt(int row, int column) {
		// Ensure our job information is populated.
		buildJobEntries();
		if ((row < 0) || (row > jobEntries.size())) {
			// Invalid row index. Nothing can be returned.
			return null;
		}
		// Obtain the job object whose data is to be returned
		final Job job = jobEntries.get(row).getJob();
		assert(job != null);
		switch (column) {
		case 0: // The job ID
			return job.getJobID();
		case 1: // return the current status.
			return job.getStatus();
		case 2: // return monitor status.
			if (job.getMonitor() != null) {
				return "Running";
			} else if (!job.isDone() && !job.isWaiting()) {
				return "Needed";
			}
			return "Not Needed";
		case 3: // return the server on which job is running.
			final Workspace ws = Workspace.get();
			final Server srvr  = ws.getServerList().getServer(job.getServerID());
			return (srvr != null) ? srvr.getName() : "<n/a>";
		case 4: // return number of CPUs 
			return "" + (job.getCPUsPerNode() * job.getNodes());
		}
		// A invalid column!
		return null;
	}

	@Override
	public void workspaceChanged(WorkspaceEvent event) {
		// Figure out the index of the entry that has changed 
		// This works fine for inserts and updates only.
		JobList jl = Workspace.get().getJobList();
		int firstRow = jl.getJobs().indexOf(event.getSource());
		int lastRow  = firstRow + 1;
		if (firstRow == -1) {
			firstRow = 0;
			lastRow  = jl.getJobs().size();
		}
		// Clear out our cached entries so that they are rebuilt
		jobEntries = null;
		// Translate workspace event to a model event.
		switch (event.getOperation().ordinal()) {
		case 0: // INSERT
			this.fireTableRowsInserted(firstRow, lastRow);
			break;
		case 1: // UPDATE
			this.fireTableRowsUpdated(firstRow, lastRow);
			break;
		case 2: // DELETE
			this.fireTableRowsDeleted(firstRow, lastRow);
			break;			
		}
	}
	
	/**
	 * Convenience method to obtain the job at a given row.
	 * 
	 * This is a convenience method that is used by different parts of
	 * the GUI to obtain the Job information from a given row.
	 *  
	 * @param row The row in the model whose associated job object
	 * is to be returned.
	 * 
	 * @return The job object associated with the row, assuming
	 * the row is valid. If not, this method returns null.
	 */
	public Job getJob(final int row) {
		buildJobEntries();
		if ((row >= 0) && (row < jobEntries.size())) {
			return jobEntries.get(row).getJob();
		}
		return null;
	}

	/**
	 * Convenience method to obtain the entry at a given row.
	 * 
	 * This is a convenience method that is used by different parts of
	 * the GUI to obtain the entry (which can contain a Job or a
	 * JobList) information from a given row.
	 *  
	 * @param row The row in the model whose associated entry
	 * is to be returned.
	 * 
	 * @return The entry associated with the row, assuming
	 * the row is valid. If not, this method returns null.
	 */
	public JobOrListWrapper getEntry(final int row) {
		buildJobEntries();
		if ((row >= 0) && (row < jobEntries.size())) {
			return jobEntries.get(row);
		}
		return null;
	}
	
	/**
	 * A generated serial version ID for serializing the information in
	 * this class (if needed).
	 */
	private static final long serialVersionUID = -5174610250739852226L;
}
