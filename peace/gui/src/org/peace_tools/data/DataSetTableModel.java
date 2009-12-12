package org.peace_tools.data;

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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;

import javax.swing.table.AbstractTableModel;

import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Workspace;
import org.peace_tools.workspace.WorkspaceEvent;
import org.peace_tools.workspace.WorkspaceListener;

public class DataSetTableModel extends AbstractTableModel implements WorkspaceListener {
	/**
	 * The default constructor.
	 * 
	 * The constructor registers this DataSetTableModel with the work
	 * space to receive notifications on Server status updates.
	 */
	public DataSetTableModel() {
		Workspace.get().addWorkspaceListener(this);
	}
	
	/**
	 * Method to obtain the columns that are to be displayed in a
	 * 
	 * 
	 * @return This method currently always returns 4.
	 */
	@Override
	public int getColumnCount() {
		return 4;
	}

	/**
	 * Method to return the number of rows to be displayed in the 
	 * table. This method counts the number of files in the various
	 * data sets and returns that number.
	 * 
	 * @param This method returns the number of servers currently in
	 * this work space.
	 */
	@Override
	public int getRowCount() {
		Workspace ws = Workspace.get();
		ArrayList<DataSet> dataSets = ws.getDataSets();
		int fileCount = 0;
		for(DataSet ds: dataSets) {
			fileCount += ds.getClusterList().size();
			fileCount += ds.getMSTList().size();
			fileCount++; // Account for EST FASTA file itself
		}
		return fileCount;
	}

	/**
	 * This is a helper method to obtain a file name corresponding
	 * to a given row  in the work space.
	 * 
	 * @param row The row in the work space for which the file name
	 * is to be returned.
	 * 
	 * @return The file name corresponding to the given row.
	 */
	private String getFileName(int row) {
		Workspace ws = Workspace.get();
		ArrayList<DataSet> dataSets = ws.getDataSets();
		for(DataSet ds: dataSets) {
			if (row == 0) {
				// Return the EST FASTA file path.
				return ds.getPath(); 
			}
			row--; // account for EST file skipped.
			if (ds.getMSTList().size() > row) {
				return ds.getMSTList().get(row).getPath();
			}
			row -= ds.getMSTList().size(); // track mst files skipped
			if (ds.getClusterList().size() > row) {
				return ds.getClusterList().get(row).getPath();
			}
			row -= ds.getClusterList().size(); // track cluster files skipped
		}
		// invalid row!
		return null;
	}
	
	/**
	 * This is a helper method to obtain a file name corresponding
	 * to a given row  in the work space.
	 * 
	 * @param row The row in the work space for which the file name
	 * is to be returned.
	 * 
	 * @return The file name corresponding to the given row.
	 */
	public Object getEntry(int row) {
		Workspace ws = Workspace.get();
		ArrayList<DataSet> dataSets = ws.getDataSets();
		for(DataSet ds: dataSets) {
			if (row == 0) {
				// Return the EST FASTA file path.
				return ds; 
			}
			row--; // account for EST file skipped.
			if (ds.getMSTList().size() > row) {
				return ds.getMSTList().get(row);
			}
			row -= ds.getMSTList().size(); // track mst files skipped
			if (ds.getClusterList().size() > row) {
				return ds.getClusterList().get(row);
			}
			row -= ds.getClusterList().size(); // track cluster files skipped
		}
		// invalid row!
		return null;
	}
	
	/**
	 * Obtain the value to be displayed at a given row and column.
	 * 
	 * @return The object to be displayed at a given row and column.
	 * The value depends on the column being displayed.
	 */
	@Override
	public Object getValueAt(int row, int column) {
		// First obtain the full path for the file at the
		// given row in the work space.
		String fileName = getFileName(row);
		if (fileName == null) {
			return null;
		}
		try {
			// Create a file object to obtain information.
			File file = new File(fileName);
			switch (column) {
			case 0: return file;
			case 1: return file.getCanonicalPath().toString();
			case 2: return new Long(file.length()).toString();
			case 3: 
				Calendar timestamp = Calendar.getInstance();
				timestamp.setTimeInMillis(file.lastModified());
				return timestamp.getTime().toString();
			}
		} catch (IOException ioe) {
			ProgrammerLog.log(ioe);
		}
		// A invalid column!
		return null;
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
	
	/**
	 * Obtain the class that describes the data type of a given 
	 * column.
	 * 
	 * @param int column The zero-based index of the column of the
	 * table for which class is required.
	 */
	@Override
	public Class<?> getColumnClass(int column) { 
		if (column == 0) {
			return File.class;
		}
		return String.class;
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
		String titles[] = {"File name", "Path", 
				"Size (bytes)", "Date Modified"};
		return titles[col];
	}
	
	/**
	 * Translates work space change notifications to table changes.
	 * 
	 * This method implements the only method in the WorkspaceListener
	 * interface to intercept events notifying changes in the work
	 * space. This method in-turn fires table changed events.
	 * 
	 * @note This method simply fires a table data changed event
	 * indicating that all the cells in the model have changed.
	 */
	@Override
	public void workspaceChanged(WorkspaceEvent event) {
		this.fireTableDataChanged();
	}
	
	/**
	 * A generated serial version ID for serializing the information in
	 * this class (if needed).
	 */
	private static final long serialVersionUID = -5174610250739852226L;
}
