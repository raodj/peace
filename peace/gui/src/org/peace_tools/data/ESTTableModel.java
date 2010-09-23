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
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

import javax.swing.table.AbstractTableModel;

import org.peace_tools.generic.FindEvent;

/**
 * A bridge class between FASTA entries in a FASTA file and a JTable.
 * 
 * This class serves as a bridge between the in-memory representation of
 * fragments in a FASTA file to the actual display in JTable. This class 
 * enables reusing/sharing fragment data and displaying it different views
 * with minimal memory footprint. The terminology used in this class is
 * from the context of the Model-View-Controller (MVC) object-oriented
 * design pattern.
 * 
 * <p><b>Note:</b>  The number of columns displayed by this table model can be
 * varied based on the column width set via a call to setColumnSize() method.
 * By default this model provides only the following 3 columns:
 * 
 *  <ol>
 *  
 *  <li>The first column indicates if the entry has been modified. The column
 *  appears as a custom icon/button when the user can revert the choice.</li>
 *  
 *  <li>The second column is a boolean value that indicates if the row has been
 *  selected for export. This column appears as a check box.</li>
 *  
 *  <li>The third column contains the index of the EST being displayed.</li>
 *  
 *  <li>The fourth column contains a modifiable FASTA header for the sequence.</li>
 *  
 *  <li>The fifth (and remainder) column contains the modifiable nucleotide sequence.</li>
 *  
 *  </ol>
 *   
 */
public class ESTTableModel extends AbstractTableModel {
	/**
	 * The default constructor.
	 * 
	 * The constructor computes the core information required for quick 
	 * operation of this model class and sets up various instance variables
	 * to appropriate values.
	 * 
	 * @param estList The list of fragments from a FASTA file to be used and
	 * suitably exposed by this model class. This parameter cannot be null.
	 * 
	 * @param miniView In miniView Only three columns are reported by the model.
	 * In addition, the model entries cannot be edited. This mode is used by the
	 * OverlapView to permit users to conveniently select ESTs to be highlighted
	 * in that view.
	 */
	public ESTTableModel(ESTList estList, boolean miniView) {
		assert ( estList != null );
		this.estList     = estList;
		this.basesPerCol = -1;
		this.miniView    = miniView;
		// Compute the maximum length of ESTs to determine column count.
		int maxLen = 0;
		for(EST est: estList.getESTs()) {
			maxLen = Math.max(maxLen, est.getSequence().length());
		}
		// Save max len for future use. It will change when EST sequences
		// are edited by the user.
		this.maxESTLen = maxLen;
		// Setup the initial selection list to default 'false' values
		this.selectionList = new ArrayList<Boolean>(estList.getESTs().size());
		for(int i = 0; (i < estList.getESTs().size()); i++) {
			selectionList.add(Boolean.FALSE);
		}
		// Create an empty hash map without any modified entries in it
		this.modifiedESTs = new HashMap<Integer, EST>();
	}

	/**
	 * Method to return the number of rows to be displayed in the 
	 * table.
	 * 
	 * @return This method returns the number of entries in the FASTA
	 * file to be displayed in a table.
	 */
	@Override
	public int getRowCount() {
		return estList.getESTs().size();
	}
 	
	/**
	 * This method overrides the interface method in TableModel to 
	 * return the number of columns in the Table. This value is
	 * determined using a combination of maxLen and basesPerCol values.
	 * 
	 * @return The number of columns to be displayed in this table. There
	 * are always a minimum of 3 columns. In miniView the number of columns
	 * is fixed to three and the basesPerCol value is ignored.
	 */
	public int getColumnCount() {
		if (this.miniView) {
			// In miniView we always report only 3 columns
			return 3;
		}
		int colCount = 4;
		if (basesPerCol > 0) {
			colCount += (int) Math.ceil(maxESTLen / (double) basesPerCol);
		} else {
			colCount++;
		}
		return colCount;
	}

	/**
	 * This method overrides the interface method in RowModel to 
	 * return the title for the columns. The column titles are dynamically
	 * determined depending on the number of columns displayed in the table.
	 * 
	 * @return The title for the columns. The title for the first two columns
	 * that contains the fragment index and fasta header information is fixed.
	 */	
	public String getColumnName(int column) {
		if (miniView) {
			// In miniView we skip over the first column with revert button
			// as the model entries cannot be edited.
			column++;
		}
		if ((column == 0) || (column == 1)) {
			return "";
		} else if (column == 2) {
			return "Index";
		} else if (column == 3) {
			return "FASTA Header";
		} else if (basesPerCol == -1) {
			return "Full Fragment";
		} else {
			// Compute title based on the bases per column value.
			column -= 2;
			return String.format("Frag: %d - %d", column * basesPerCol, 
				(column + 1) * basesPerCol);
		}
	}

	/**
	 * This method returns the subset of base pairs to be displayed in a
	 * given column for a given row. This method overrides the 
	 * interface method in TableModel to return the data to be displayed
	 * in a given row. The return value is computed as follows:
	 * 
	 * <ol>
	 * 
	 * <li>For column 0, this method returns a Boolean value indicating if
	 * the entry has been modified.</li>
	 * 
	 * <li>For column 1, this method returns a Boolean value indicating if
	 * the entry has been selected.</li>
	 * 
	 * <li>For column 2, this method returns the logical index of the 
	 * entry in the FASTA file.</li>
	 * 
	 * <li>For column 3, this method always returns the FASTA identifier
	 * associated with the entry in the given row.</li>
	 * 
	 * <li>For other columns, if basesPerEST == -1, then this method returns
	 * the complete fasta sequence to be displayed in the second column.</li>
	 * 
	 * <li>If the node is an EST node and column * basesPerCol is less than
	 * EST length, then the bases logically associated with the column are 
	 * returned.</li>
	 * 
	 * <li>Otherwise this method returns an empty string.</li>
	 * 
	 * </ol>
	 * 
	 * @param row The row for which the data is desired.
	 * 
	 * @param col The column for which the data is desired.
	 * 
	 * @return The data to be displayed for a given node (or row) and a 
	 * given column. 
	 */
	@Override
	public Object getValueAt(int row, int col) {
		if (miniView) {
			// In miniView we skip over the first column with revert button
			// as the model entries cannot be edited.
			col++;
		}
		if ((row < 0) || (row > estList.getESTs().size()) ||
			(col < 0) || (col > this.getColumnCount())) {
			return "";
		}
		// Obtain the entry for the given row while heeding sorting
		if (sortedIndexs != null) {
			row = sortedIndexs[row];
		}
		// Get original entry.
		EST est = estList.getESTs().get(row);
		// Check and get modified entry for this EST (if any)
		EST newEST = this.modifiedESTs.get(est.getID());
		// Use the modified EST we have a modified entry for it.
		est = (newEST != null) ? newEST : est;
		// Return the appropriate object based on the column
		if (col == 0) {
			// Return boolean flag indicating if this EST is modified
			return new Boolean(newEST != null);			
		} else if (col == 1) {
			// Return boolean flag indicating if this EST is selected
			return this.selectionList.get(row);
		} else if (col == 2) {
			return est.getID();
		} else if (col == 3) {
			return est.getInfo();
		} else if (basesPerCol == -1) {
			// The full fragment is displayed in one column.
			return est.getSequence();
		}
		// Return substring of the sequence based on basesPerCol and
		// the col for which data is requested.
		col -= 2;
		final String seq = est.getSequence();
		final int start  = (col * basesPerCol);
		final int stop   = Math.min(start + basesPerCol, seq.length());
		String retVal    = "";
		if (start < seq.length()) {
			retVal = seq.substring(start, stop);
		}
		return retVal;
	}

	/**
	 * Obtain the EST entry at a given row.
	 * 
	 * @param row The logical row from where the EST is to be retrieved.
	 * Note that this method pays heed to the sorted order. So the row is
	 * the logical row in the sorted list of entries and not the absolute
	 * index in the EST list.
	 * 
	 * @return The EST at the given logical row.
	 */
	public EST getESTAt(int row) {
		// Obtain the entry for the given row while heeding sorting
		if (sortedIndexs != null) {
			row = sortedIndexs[row];
		}
		// Get original entry value first to look up EST index.
		EST est = estList.getESTs().get(row);
		// Check and pick modified entry (if any)
		EST newEST = this.modifiedESTs.get(est.getID());
		// Return with preference to the newly modified entry
		return (newEST != null) ? newEST : est;
	}
	
	/**
	 * Method to change the order in which fragments are ordered.
	 * 
	 * This method resets the order in which fragments are presented
	 * by this class to the "view" that is displaying the FASTA
	 * fragments.
	 * 
	 * @param order The order in which the clusters are to be sorted. The 
	 * following values are valid for this parameter: 0: no sorting,
	 * 1: shorter fragments first, 2: shorter fragments last.
	 * 
	 */
	public void sort(final int order) {
		if (order == 0) {
			// Reset sorting.
			this.sortedIndexs = null;
			// Fire notification
			fireTableStructureChanged();
			return;
		}
		// The comparator object to be used for sorting below.
		Comparator<Integer> comp = new Comparator<Integer>() {
			@Override
			public int compare(Integer param1, Integer param2) {
				final EST frag1 = estList.getESTs().get(param1);
				final EST frag2 = estList.getESTs().get(param2);
				// Compute difference.
				final int diff  = frag1.getSequence().length() - 
					frag2.getSequence().length();
				// Interpret difference based on order.
				return (order == 1) ? diff : -diff;
			}
		};
		// Create a default list of values.
		final int ESTCount = estList.getESTs().size();
		sortedIndexs = new Integer[ESTCount];
		for(int i = 0; (i < ESTCount); i++) {
			sortedIndexs[i] = new Integer(i);
		}
		// Sort the list based on the sort criteria
		Arrays.<Integer>sort(sortedIndexs, comp);
		// Fire notification
		fireTableStructureChanged();
	}
	
	/**
	 * Method to change the number of bases to be displayed in each column.
	 * 
	 * This method must be used to change the number of bases per column to 
	 * be presented my this model class.
	 * 
	 * @param basesPerCol If this parameter is -1, then all sequences are 
	 * presented as a single column. Otherwise the sequences are presented
	 * as multiple columns, with each column having no more than the 
	 * specified number of neucleotides.
	 */
	public void setBasesPerCol(int basesPerCol) {
		this.basesPerCol = basesPerCol;
		fireTableStructureChanged();
	}
	
	/**
	 * Obtain the class that describes the data type of a given 
	 * column.
	 * 
	 * @param column The zero-based index of the column whose
	 * Class type is to be returned.
	 * 
	 * @return This method always returns Boolean.class for the 
	 * first two columns and for all other columns it returns 
	 * String.class as the class of the data as this model 
	 * exposes all the information as Strings.
	 */
	@Override
	public Class<?> getColumnClass(int column) {
		if (miniView) {
			// In miniView we skip over the first column with revert button
			// as the model entries cannot be edited.
			column++;
		}
		return (column < 2) ? Boolean.class : String.class;
	}
	
	/**
	 * Interface method to determine if an entry in the JTable is 
	 * editable.
	 * 
	 * @return This method return true as long as the user is not
	 * trying to edit the EST ID (or index) values to indicate that
	 * the data in the FASTA file cannot be edited.
	 */
	@Override
	public boolean isCellEditable(int row, int column) {
		return (miniView ? (column == 0) : (column != 2));
	}
	
	/**
	 * Obtain the list of ESTs associated with this model.
	 * 
	 * @return The EST list set for use by this model.
	 */
	public ESTList getESTList() { return estList; }
	
	/**
	 * This method overrides the empty base class implementation to maintain
	 * a set of modified ESTs. The modified ESTs are maintained in a separate
	 * hash map to permit the user to revert the modifications on a per-EST
	 * basis. This method appropriately interprets the changes depending
	 * on the column to which the changes have been made.
	 * 
	 * @param value The new value to be assigned to the specified entry.
	 * @param row The row in the table that has been modified.
	 * @param column The column whose value has been modified.
	 */
	@Override
	public void setValueAt(Object value, int row, int column) {
		if (miniView) {
			// In miniView mode entries can be merely selected.
			if (column == 0) {
				column++;
			} else {
				// In miniView mode other columns cannot be modified.
				return;
			}
		}
		// Obtain the entry for the given row while heeding sorting
		if (sortedIndexs != null) {
			row = sortedIndexs[row];
		}
		// Obtain the actual ID (.ie. index) of the EST being modified
		final int estIndex = estList.getESTs().get(row).getID();
		// The operations depend on the column which has been modified.
		if (column == 0) {
			// The user has managed to click the revert button. In this
			// case we simply clear out the entry in the modified ESTs 
			// hash map.
			modifiedESTs.remove(estIndex);
		} else if (column == 1) {
			// The user has toggled selection. 
			selectionList.set(row, !selectionList.get(row));
		} else if (column != 2) {
			// When control drops here that means the user has modified the
			// FASTA header or the sequence. We need to create and/or manage a
			// modified entry as needed.
			EST currEST = modifiedESTs.get(estIndex);
			if (currEST == null) {
				// Use the original value as the entry
				currEST = estList.getESTs().get(row);
			}
			// Appropriately obtain the necessary updates
			String fastaHeader = (column == 3) ? value.toString() : currEST.getInfo();
			String sequence    = currEST.getSequence();
			if (column > 3) {
				// Sequence handling is a bit more tricky as sequences can be displayed
				// as multi-column value
				if (this.basesPerCol == -1) {
					// Sequence is being displayed as a single column
					sequence = value.toString();
				} else {
					// Splice in the modified value at the appropriate point in the
					// sequence. 
					final int col = column - 2;
					final String seq = currEST.getSequence();
					final int start  = (col * basesPerCol);
					final int stop   = Math.min(start + basesPerCol, seq.length());
					if (start < seq.length()) {
						sequence = seq.substring(0, start) + value.toString() + seq.substring(stop + 1);
					}
				}
			}
			// Now update the EST entry in the modifiedESTs hash map
			EST newEST = new EST(estIndex, fastaHeader, sequence);
			modifiedESTs.put(estIndex, newEST);
			// Update the maximum EST length if it has changed.
			if (maxESTLen < sequence.length()) {
				// Yup. Max length has changed. The whole table may have to rerendered.
				maxESTLen = sequence.length();
				fireTableStructureChanged();
			}
		}
		// In addition the 'revert' status for this entry has changed.
		// Fire an event to the listeners to ensure the entry is properly
		// refreshed in the GUI view.
		fireTableRowsUpdated(row, row);
	}
	
	/**
	 * Method to determine the set of entries that have been selected by the user.
	 * 
	 * This method can be used to obtain the list of selected entries in the model.
	 * This method is synonymous to JTable.getSelectedRows() except that it directly
	 * deals with the selection information stored in the model rather than with 
	 * the table selection model.
	 * 
	 * @return A native array containing the list of rows that have been selected.
	 * This method always returns a non-null array. However, if no entries are 
	 * selected then the return value is an array with zero entries.
	 */
	public int[] getSelectedRows() {
		// First determine the selected rows in a array list
		ArrayList<Integer> selectedRows = new ArrayList<Integer>(estList.getESTs().size());
		for(int row = 0; (row < estList.getESTs().size()); row++) {
			// Handle mapping when list is logically sorted.
			final int actualRow = (sortedIndexs != null) ? sortedIndexs[row] : row;
			if (selectionList.get(actualRow)) {
				selectedRows.add(row);
			}
		}
		// Now convert the array list to a native array for interface compliance
		int[] selections = new int[selectedRows.size()];
		for(int entry = 0; (entry < selections.length); entry++) {
			selections[entry] = selectedRows.get(entry);
		}
		// Return the native array containing selected rows
		return selections;
	}
	
	/**
	 * Helper method to search through the list of ESTs for a given substring
	 * or a regular expression. This helper method is used by the view to
	 * perform "Find" operations in this data set.
	 * 
	 * @param fe The find event from which the data for the search is extracted
	 * and used.
	 * 
	 * @param startRow The starting row for the search. This is typically the
	 * currently selected row in the table.
	 * 
	 * @return This method returns the index of the row where a match for the
	 * given find request was encountered. If a valid match was not found, then
	 * this method returns -1.
	 */
	public int find(final FindEvent fe, final int startRow) {
		final int direction    = fe.isForwardSearch() ? 1 : -1; // row increment.
		int  currRow           = startRow + direction;
		final int MaxRow       = estList.getESTs().size();
		
		// Setup search string depending on type of search requested by user
		String searchStr = null;
		if (!fe.isRegExSearch()) {
			searchStr = fe.getSearchString();
			if (!fe.isCaseSensitive()) {
				searchStr = searchStr.toLowerCase();
			}
		}
		// Search row-by-row
		while ((currRow != startRow) && (currRow >= 0) && (currRow < MaxRow)) {
			final EST currEST = getESTAt(currRow);
			boolean match = (searchStr != null) ? currEST.contains(searchStr, fe.isCaseSensitive()) : 
				currEST.contains(fe.getRegExPattern());
			if (match) {
				return currRow;
			}
			// Now match found in this row. Onto the next row.
			currRow += direction;
			// Handle the wrap around cases.
			if (fe.isWrapAround()) {
				if (currRow < 0) {
					currRow = MaxRow - 1; 
				} else if (currRow >= MaxRow) {
					currRow = 0;
				}
			}
		}
		return -1;
	}
	
	/**
	 * Reference to the list of fragments that is actually exposed by this
	 * table model. This value is set in the constructor and is never changed
	 * during the life time of this class.
	 */
	private final ESTList estList;
	
	/**
	 * The maximum length of an EST sequence to be adapted by this model.
	 * This value along with the basesPerCol determines the total number of 
	 * columns that are logically represented by this model. 
	 */
	private int maxESTLen;
	
	/**
	 * The number of base pairs per column. This value along with the
	 * maxESTLen determines the total number of columns that are logically
	 * represented by this model. 
	 */
	private int basesPerCol;
	
	/**
	 * This array contains a sorted list of indexes of fragments if a sorting
	 * scheme has been applied to the data set. If this array is not null then
	 * this table model returns entries using the order specified in this
	 * array.  The entries in this array are indexes into the original set of
	 * fragments contained in the estList instance variable.
	 *  clusters in the root. This list changes depending on the order of 
	 *  sorting requested by the user.
	 */
	private Integer[] sortedIndexs = null;
	
	/**
	 * An array list that is used to maintain information on whether the given
	 * row is currently selected by the user. Initially this array is set to
	 * all false values. The values are displayed as a check box that the
	 * user can check/clear to select an entry to be saved.
	 */
	private ArrayList<Boolean> selectionList;

	/**
	 * A hash map that is used to maintain information on ESTs that have 
	 * been modified. These modified ESTs are maintained in this hash map.
	 * Whenever sequence data is requested, this hash map is used to report
	 * the updated information. The entries are also used to display a "revert"
	 * button in the second column in the table.
	 */
	private HashMap<Integer, EST> modifiedESTs;

	
	/**
	 * Flag to indicate if this model is acting in 'mini view' mode. In miniView 
	 * (this value is true) only three columns are reported by the model. In 
	 * addition, the model entries cannot be edited. This mode is used by the
	 * OverlapView to permit users to conveniently select ESTs to be highlighted
	 * in that view.
	 */
	private final boolean miniView;
	
	/**
	 * A generated serialization UID included just to keep the compiler happy.
	 */
	private static final long serialVersionUID = -3412623208152260577L;
}
