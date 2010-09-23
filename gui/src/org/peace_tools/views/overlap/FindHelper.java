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

import java.util.ArrayList;
import java.util.regex.Pattern;

import org.peace_tools.data.EST;
import org.peace_tools.data.ESTEntry;
import org.peace_tools.data.OverlapModel;
import org.peace_tools.generic.FindEvent;
import org.peace_tools.generic.FindListener;

/**
 * A relatively straightforward helper class to assist with various types
 * of find/search operations.
 *
 * <p>The {@link OverlapView} and its primary sub-components (such as: 
 * {@link HighlightPanel}, {@link OverlapPanel}) perform different forms of
 * find operations to aid the user quickly locate fragments in various related
 * presentations of the data. This class serves to house the different flavors
 * of methods used to achieve the various forms of find operations required by
 * the user. An instance of this class is created by the {@link OverlapView} 
 * and is shared by the various sub-components. This class does not have
 * any GUI components directly associated with it. However, the various GUI
 * components utilize this class for some of their operations.</p>
 * 
 * <p><b>Note</b>This class is meant to be instantiated and used
 * only by the OverlapView class. Consequently it has been made 
 * package private.</p>
 */
class FindHelper implements FindListener {
	/**
	 * The only constructor for this class.
	 * 
	 * The constructor is typically invoked by the {@link OverlapView} class.
	 * The constructor merely initializes the various instance variables to
	 * their default initial value.
	 * 
	 * @param model
	 */
	FindHelper(OverlapModel model) {
		overlapModel = model;
		// Set search parameters to beginning of overlap model.
		prevFindRow  = 0;
		prevFindCol  = 0;
		prevESTIndex = 0;
	}
	
	/**
	 * Interface method to locate entry associated with a given EST index.
	 * 
	 * This method is a convenience method that is primarily used by the
	 * {@link OverlapPanel} to highlight the entry associated with the 
	 * given estIndex.  This method searches from the beginning all the way down
	 * to the end.
	 * 
	 * @param estIndex The index of the EST whose entry is to be located.
	 * 
	 * @return This method returns a non-null entry if it was found.
	 * Otherwise it returns null.
	 */
	public ESTEntry find(int estIndex) {
		final String searchStr = "" + estIndex;
		ESTEntry entry = null;
		// Repeatedly search until we find the correct entry. We need this
		// loop because searching for est index as a string can yield a hit
		// with another EST that has the estIndex string as a substring.
		prevFindRow = 0;
		prevFindCol = -1;
		do {
			prevFindCol++;
			entry =  find(prevFindRow, prevFindCol, true, false, searchStr, false, null);
		} while ((entry != null) && (entry.getEST().getID() != estIndex));
		return entry;
	}
	
	/**
	 * A general purpose search method that is typically used by the 
	 * {@link OverlapView} to handle requests from the generic FindDialog. This
	 * method commences the search from the location where the previous search
	 * terminated.
	 * 
	 * @param searchStr The search string to search. This entry cannot be null or an
	 * empty string.
	 * @param caseSensitive Flag to indicate if the search must be case sensitive. If
	 * this flag is true, then the search is performed in a case sensitive manner.
	 * @param direction The direction in which the search must proceed. If this flag is
	 * true, then the search proceeds towards the bottom of the overlap model. 
	 * On the other hand, if this flag is false, then the search proceeds towards
	 * the top of the model.
	 * @param wrapAround If this flag is true, then the search wraps around when it hits
	 * the logical end of the model in either direction.
	 * 
	 * @return If a matching entry is found then this method returns the entry after
	 * updating {@link #prevFindCol} and {@link #prevFindRow} (which can be accessed via the
	 * {@link #getLastFindCol()} and {@link #getLastFindRow()} API methods). Otherwise it 
	 * returns null without affecting the aforementioned instance variables.
	 */
	public ESTEntry find(String searchStr, boolean caseSensitive, 
			boolean direction, boolean wrapAround) {
		return find(prevFindRow, prevFindCol, direction, wrapAround, 
					searchStr, caseSensitive, null);
	}
	
	/**
	 * A general purpose search method that is typically used by the 
	 * {@link OverlapView} to handle requests from the generic FindDialog. This
	 * method commences the search from the location where the previous search
	 * terminated.
	 * 
	 * @param regExp The regular expression to be used for finding a matching entry.
	 * This entry cannot be null.
	 * @param direction The direction in which the search must proceed. If this flag is
	 * true, then the search proceeds towards the bottom of the overlap model. 
	 * On the other hand, if this flag is false, then the search proceeds towards
	 * the top of the model.
	 * @param wrapAround If this flag is true, then the search wraps around when it hits
	 * the logical end of the model in either direction.
	 * 
	 * @return If a matching entry is found then this method returns the entry after
	 * updating {@link #prevFindCol} and {@link #prevFindRow} (which can be accessed via the
	 * {@link #getLastFindCol()} and {@link #getLastFindRow()} API methods). Otherwise it 
	 * returns null without affecting the aforementioned instance variables.
	 */
	public ESTEntry find(Pattern regExp, boolean direction, boolean wrapAround) {
		return find(prevFindRow, prevFindCol, direction, wrapAround, 
					null, true, regExp);
	}

	/**
	 * Determine the row where the previous find operation found a successful match.
	 * 
	 * @return This method returns the row where the previous find operation
	 * found a successful match.
	 */
	public int getLastFindRow() {
		return this.prevFindRow;
	}

	/**
	 * Determine the row where the previous find operation found a successful match.
	 * 
	 * @return This method returns the column where the previous find operation
	 * found a successful match.
	 */
	public int getLastFindCol() {
		return this.prevFindCol;
	}

	/**
	 * Determine the index of the currently active EST.
	 * 
	 * @return This method returns the index (or ID) of the currently active EST.
	 * This is also the index of the EST that was located by the last successful
	 * search/find operation performed by this class.
	 */
	public int getCurrESTIndex() {
		return this.prevESTIndex;
	}

	@Override
	public boolean find(FindEvent event) {
		ESTEntry entry = (event.isRegExSearch() ? 
				find(event.getRegExPattern(), event.isForwardSearch(), event.isWrapAround()) :
				find(event.getRegExPattern(), event.isForwardSearch(), event.isWrapAround()));
		return (entry != null);
	}
	
	/**
	 * The core search method that is called by various API methods to search 
	 * multiple rows.
	 * 
	 * @param startRow The starting logical row in the overlap model from where
	 * the search must commence.
	 * @param startCol The column (or index) entry within the row from where the
	 * search is to commence.
	 * @param direction The direction in which the search must proceed. If this flag is
	 * true, then the search proceeds towards the bottom of the overlap model. 
	 * On the other hand, if this flag is false, then the search proceeds towards
	 * the top of the model.
	 * @param wrapAround If this flag is true, then the search wraps around when it hits
	 * the logical end of the model in either direction.
	 * @param searchStr The search string to search. If this entry is null then the
	 * regular expression parameter must not be null.
	 * @param caseSensitive Flag to indicate if the search must be case sensitive. If
	 * this flag is true, then the search is performed in a case sensitive manner. This
	 * flag is used only when searching with a simple string (that is searchStr parameter
	 * is not null).
	 * @param regExp The regular expression to be used for search. If this entry is null
	 * then the searchStr parameter must not be null.
	 * @return If a matching entry is found then this method returns the entry after
	 * updating {@link #prevFindCol} and {@link #prevFindRow}. Otherwise it returns null
	 * without affecting the aforementioned instance variables.
	 */
	private ESTEntry find(int startRow, int startCol, boolean direction, boolean wrapAround,
			String searchStr, boolean caseSensitive, Pattern regExp) {
		// Short cut to maximum number of rows for use in this method to streamline code
		final int MaxRow = overlapModel.getMaxRow();
		// Determine the end row based on direction and overlap flags.
		final int endRow = (wrapAround ? startRow : (direction ? MaxRow - 1 : 0));
		// Now determine the ending column on the endRow based on overlap and direction flags
		final int endCol = (wrapAround ? startCol : (direction ? overlapModel.getRowSize(endRow) : 0));
		// Now search through the data until we hit the end point
		int currRow = startRow, currCol = startCol, entryIndex = -1;
		while ((currRow != endRow) || (currCol != endCol)) {
			// Search in the current row in appropriate direction for match
			ArrayList<ESTEntry> rowList = overlapModel.getRow(currRow);
			final int endColInRow = (direction ? rowList.size() : 0);
			entryIndex = find(rowList, Math.min(currCol, endColInRow), 
					Math.max(currCol, endColInRow), searchStr, caseSensitive, regExp);
			if (entryIndex != -1) {
				// Found a match. 
				break;
			}
			// No entries were found on the current row. Onto the next row in appropriate
			// direction.
			currRow += (direction ? 1 : -1);
			if (((currRow < 0) || (currRow >= MaxRow)) && wrapAround) {
				currRow = (currRow < 0) ? MaxRow - 1 : 0;
			}
			// Setup the current column on the currRow depending on the model
			currCol = (direction ? 0 : overlapModel.getRowSize(currRow));
		}
		
		if (entryIndex == -1) {
			// Did not found a valid entry.
			return null;
		}
		// Found a valid entry. Update necessary information and return entry to the
		// caller
		// Return the entry in the current row and column back to the caller.
		prevFindCol  = entryIndex;
		prevFindRow  = currRow;
		final ArrayList<ESTEntry> rowList = overlapModel.getRow(prevFindRow);
		final ESTEntry entry = rowList.get(prevFindCol);
		prevESTIndex = entry.getEST().getID();
		return entry;
	}
	
	/**
	 * The core search method that is called by various methods to search within a 
	 * given row.
	 * 
	 * @param fragmentList The list of entries on a given row to be searched for.
	 * @param startCol The column (or index) entry within the row from where the
	 * search is to commence.
	 * @param endCol The column (or index) entry within the row where the
	 * search is to end.
	 * @param searchStr The search string to search. If this entry is null then the
	 * regular expression parameter must not be null.
	 * @param caseSensitive Flag to indicate if the search must be case sensitive. If
	 * this flag is true, then the search is performed in a case sensitive manner. This
	 * flag is used only when searching with a simple string (that is searchStr parameter
	 * is not null).
	 * @param regExp The regular expression to be used for search. If this entry is null
	 * then the searchStr parameter must not be null.
	 * @return If a matching entry is found then this method returns the column where
	 * the entry was found. Otherwise it returns -1.
	 */
	private int find(ArrayList<ESTEntry> fragmentList,
			int startCol, int endCol,
			String searchStr, boolean caseSensitive, Pattern regExp) {
		for(int col = startCol; (col < endCol); col++) {
			final ESTEntry fragment = fragmentList.get(col);
			final EST est = fragment.getEST();
			boolean found = (searchStr != null) ? est.contains(searchStr, caseSensitive) :
				                                  est.contains(regExp);
			if (found) {
				// This is a matching entry
				return col;
			}
		}
		// No matching entries found.
		return -1;
	}
	
	/**
	 * Instance variable to track the row where the previous find operation 
	 * was completed. This value is updated by various methods to track the
	 * location of the previous operation so that subsequent operations can
	 * commence from that row. The row is the primary entry point into the
	 * {@link #overlapModel} object associated with this class.
	 */
	private int prevFindRow;
	
	/**
	 * Instance variable to track the column where the previous find operation 
	 * was completed. This value is updated by various methods to track the
	 * location of the previous operation so that subsequent operations can
	 * commence from that column. The column is index into a row of the 
	 * {@link #overlapModel} object associated with this class.
	 */
	private int prevFindCol;
	
	/**
	 * The index of the last successfully found EST entry. This value is used
	 * by the {@link OverlapPanel} to suitably highlight the currently active
	 * entry in the overlap view.
	 */
	private int prevESTIndex;
	
	/**
	 * The OverlapModel that provides the necessary information
	 * for searching the core information associated with the various
	 * EST fragments. This instance variable is initialized by the 
	 * constructor and is never changed during the life time of 
	 * this object.
	 */
	private final OverlapModel overlapModel;
}
