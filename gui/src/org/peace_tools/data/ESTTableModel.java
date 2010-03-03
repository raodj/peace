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

import java.util.Arrays;
import java.util.Comparator;

import javax.swing.table.AbstractTableModel;

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
 * By default this model provides only 2 columns.
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
	 */
	public ESTTableModel(ESTList estList) {
		assert ( estList != null );
		this.estList     = estList;
		this.basesPerCol = -1;
		// Compute the maximum length of ESTs to determine column count.
		int maxLen = 0;
		for(EST est: estList.getESTs()) {
			maxLen = Math.max(maxLen, est.getSequence().length());
		}
		// Save max len as it it will never change.
		this.maxESTLen = maxLen;
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
	 * are always a minimum of 3 columns.
	 */
	public int getColumnCount() {
		int colCount = 2;
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
		if (column == 0) {
			return "Index";
		} else if (column == 1) {
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
	 * <li>For column 0, this method returns the logical index of the 
	 * entry in the FASTA file.</li>
	 * 
	 * <li>For column 1, this method always returns the FASTA identifier
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
		if ((row < 0) || (row > estList.getESTs().size()) ||
			(col < 0) || (col > this.getColumnCount())) {
			return "";
		}
		// Obtain the entry for the given row while heeding sorting
		if (sortedIndexs != null) {
			row = sortedIndexs[row];
		}
		EST est = estList.getESTs().get(row);
		if (col == 0) {
			return est.getID();
		} else if (col == 1) {
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
		return estList.getESTs().get(row);
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
	 * @return This method always returns String.class as the class
	 * of the data as this model exposes all the information as
	 * a String.
	 */
	@Override
	public Class<?> getColumnClass(int column) { 
		return String.class;
	}
	
	/**
	 * Interface method to determine if an entry in the JTable is 
	 * editable.
	 * 
	 * @return This method always returns false to indicate that
	 * the data in the FASTA file cannot be edited.
	 */
	@Override
	public boolean isCellEditable(int row, int column) {
		return false;
	}
	
	/**
	 * Obtain the list of ESTs associated with this model.
	 * 
	 * @return The EST list set for use by this model.
	 */
	public ESTList getESTList() { return estList; }
	
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
	private final int maxESTLen;
	
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
	 * A generated serialization UID included just to keep the compiler happy.
	 */
	private static final long serialVersionUID = -3412623208152260577L;
}
