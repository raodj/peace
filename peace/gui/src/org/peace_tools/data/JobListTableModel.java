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

import javax.swing.table.AbstractTableModel;

import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;
import org.peace_tools.workspace.WorkspaceEvent;
import org.peace_tools.workspace.WorkspaceListener;

/**
 * A bridge class between Job entries in a workspace and a JTable.
 * 
 * This class serves as a bridge between the in-memory representation of
 * jobs in a Workspace (represented by the set of classes in the 
 * <code>workspace</code> package). This class enables reusing the data
 * set hierarchy maintained by the Workspace object to display it in a
 * JTable. In addition, this class also acts to monitor and update data
 * set views.
 * 
 * @note This table model currently provides the following information
 * for each job: JobID, Server, CPUs, Status.
 */
public class JobListTableModel extends AbstractTableModel implements WorkspaceListener {
	/**
	 * The default constructor.
	 * 
	 * The constructor registers this JobListTableModel with the workspace
	 * to receive notifications on Job status updates.
	 */
	public JobListTableModel() {
		Workspace.get().addWorkspaceListener(this);
	}
	
	/**
	 * Method to obtain the columns that are to be displayed by in a
	 * Job table.
	 * 
	 * @return This method currently always returns 4.
	 */
	@Override
	public int getColumnCount() {
		return 4;
	}

	/**
	 * Method to return the number of rows to be displayed in the 
	 * Job table.
	 * 
	 * @param This method returns the number of jobs currently in
	 * this workspace.
	 */
	@Override
	public int getRowCount() {
		Workspace ws = Workspace.get();
		return ws.getJobList().getJobs().size();
	}

	/**
	 * This is a convenience method to obtain the job data.
	 * 
	 * @param row The row id whose Job entry is to be returned.
	 * @return The job entry corresponding to a given row. If the
	 * entry is unavailable, this method returns null.
	 */
	public Job getJob(int row) {
		Workspace ws = Workspace.get();
		if (ws.getJobList().getJobs().size() < row) {
			return null;
		}
		// Obtain the job object whose data is to be returned
		return ws.getJobList().getJobs().get(row);
	}
	
	/**
	 * Obtain the value to be displayed at a given row and column.
	 * 
	 * @return The object to be displayed at a given row and column.
	 * The value depends on the column being displayed.
	 */
	@Override
	public Object getValueAt(int row, int column) {
		Workspace ws = Workspace.get();
		if (ws.getJobList().getJobs().size() < row) {
			return null;
		}
		// Obtain the job object whose data is to be returned
		Job job = ws.getJobList().getJobs().get(row);
		Server srvr = ws.getServerList().getServer(job.getServerID());
		switch (column) {
		case 0: // The job ID
			return job.getJobID();
		case 1: // return the current status.
			return job.getStatus();
		case 2: // return the server on which job is running.
			return (srvr != null) ? srvr.getName() : "<n/a>";
		case 3: // return number of CPUs 
			return "" + (job.getCPUsPerNode() * job.getNodes());
		}
		// A invalid column!
		return null;
	}
	
	/**
	 * Obtain the class that describes the data type of a given 
	 * column.
	 */
	@Override
	public Class<?> getColumnClass(int column) { 
		if (column == 1) {
			return JobBase.JobStatusType.class;
		}
		return String.class;
	}
	
	/**
	 * Interface method to determine if an entry in the JTable is 
	 * editable.
	 * 
	 * @return This method always returns false to indicate that
	 * the data in the Job table is not editable.
	 */
	@Override
	public boolean isCellEditable(int row, int column) {
		return false;
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
		String titles[] = {"Job ID", "Status", "Server", "CPUs"};
		return titles[col];
	}
	
	/**
	 * A generated serial version ID for serializing the information in
	 * this class (if needed).
	 */
	private static final long serialVersionUID = -5174610250739852226L;
}
