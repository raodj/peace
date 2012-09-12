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

package org.peace_tools.decagon.helpers;

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.LayoutManager2;
import java.util.ArrayList;

/**
 * This layout manager organizes components into two columns with each row
 * of varying size.
 * 
 * This layout manager is a custom layout manager that has been developed
 * to avoid unnecessary complications involved in using existing Java
 * layout manager to accomplish laying out input parameters in a two
 * column format. The first column typically contains labels of varying
 * widths, while the second column contains input elements (like spinners,
 * text fields etc.). However, the second column is optional. This
 * layout manager operates in the following manner:
 * 
 * <ol>
 * 
 * <li>Initially the layout does not have any rows. New rows are
 * automatically created as components are added to the layout. The 
 * last row in the layout is always the <i>active</i> row.</li>
 * 
 * <li>When a component is added to a given column and the <i>active</i>
 * row already has a component in that column then a new row is added
 * and the <i>active</i> row is changed to the new row.</li>
 * 
 * <li>Components can be added to the first or the second column.</li>
 * 
 * <li>The component with the widest preferred size in each column 
 * determines the width of the column. If a row has only first column
 * component, then this component occupies the full row.</li>
 * 
 * <li>The overall size of the layout remains the same and the layout
 * does not grow (or shrink) to occupy extra space. In other words,
 * {@link #preferredLayoutSize(Container)}, {@link #maximumLayoutSize(Container)}
 * and {@link #minimumLayoutSize(Container)} all return the same
 * dimensions for this layout.</li> 
 * 
 * </ol>
 *
 * <b>NOTE:</b> Refer to the documentation for LayoutManager2 for 
 * API contract for the various methods implemented by this class.
 */
public class TwoColumnLayout implements LayoutManager2 {

	/**
	 * Constant that is used as a constraint to add an component to the
	 * first column in the <i>active</i> row in this layout manager.
	 */
	public static final String FIRST_COLUMN = "column1";
	
	/**
	 * Constant that is used as a constraint to add an component to the
	 * second column in the <i>active</i> row in this layout manager.
	 */
	public static final String SECOND_COLUMN = "column2";

	/**
	 * A simple zero dimension object that is used to streamline the
	 * math in {@link Row#getDimensions(Dimension)} method below.
	 */
	private static final Dimension ZERO = new Dimension(0, 0);
	
	/**
	 * A internal class to encapsulate the two components in a given row.
	 * 
	 * This class is used as a convenience class to merely encapsulate
	 * the two components in each row of the TwoColumnLayout. Each entry
	 * has two optional components. These components are directly 
	 * accessed and used by the TwoColumnLayout. However, few additional 
	 * methods are included here to streamline the code.
	 */
	private class Row {
		/**
		 * A convenience constructor.
		 * 
		 * @param c1 The component to be set for the first column in this row.
		 * @param c2 The component to be set for the second column in this row.
		 */
		public Row(Component c1, Component c2) {
			column1 = c1;
			column2 = c2;
		}
		
		/**
		 * Helper method to determine dimensions of this row.
		 * 
		 * This is a convenience method that is used to determine the
		 * height for this row as well as track the widths for the
		 * two columns. Note that this method does not account for
		 * {@link TwoColumnLayout#rowSpacing} or {@link TwoColumnLayout#colSpacing}
		 * values.
		 * 
		 * @param colSizes An array object that is used to track the 
		 * maximum width of column1, maximum width of column2, and the
		 * minimum column width (for rows with just components in column1) in
		 * index positions 0, 1, and 2 respectively.
		 *  
		 * @return This method returns the row height. The row height is
		 * the maximum of the preferred heights of column1 and column2 
		 * objects.
		 */
		public int getDimensions(int[] colSizes) {
			final Dimension col1Dim = ((column1 != null && column1.isVisible()) ? 
					column1.getPreferredSize() : ZERO);
			final Dimension col2Dim = ((column2 != null && column2.isVisible()) ? 
					column2.getPreferredSize() : ZERO);
			// Compute the width of the first column based on column sizes.
			if (column2 == null) {
				// Scenario in which one object occupies the whole row.
				colSizes[2] = Math.max(colSizes[2], col1Dim.width);
			} else {
				// Track the widest component in first column
				colSizes[0] = Math.max(colSizes[0], col1Dim.width);
			}
			// Compute the width of the second column.
			// Track the widest component in the second column
			colSizes[1] = Math.max(colSizes[1], col2Dim.width);
			
			// Return the largest of the two rows.
			return Math.max(col1Dim.height, col2Dim.height);
		}
		
		/**
		 * The component in the first column of this row.
		 */
		public Component column1 = null;
		/**
		 * The component in the second column of this row.
		 */
		public Component column2 = null;
	}

	/**
	 * The list of rows of components maintained by this layout manager.
	 * The last row in this array list is always the <i>active</i>
	 * row. New rows are appended whenever a column to which a component
	 * is to be added already has a component in it.
	 */
	private final ArrayList<Row> componentRows = new ArrayList<Row>();
	
	/**
	 * Space to be included between adjacent rows of the layout.
	 * This value can be set via the {@link #setRowSpacing(int)}
	 * method. The current value can be obtained via the 
	 * {@link #getRowSpacing()} method.
	 */
	private int rowSpacing = 0;

	/**
	 * Space to be included between adjacent columns of the layout.
	 * This value can be set via the {@link #setColSpacing(int)}
	 * method. The current value can be obtained via the 
	 * {@link #getColSpacing() method.
	 */
	private int colSpacing = 0;
	
	/**
	 * The default constructor.
	 * 
	 * The constructor sets {@link #rowSpacing} and {@link #colSpacing}
	 * to zero.
	 */
	public TwoColumnLayout() {
		// nothing to be done for now.
	}
	
	/**
	 * Constructor to set row and column spacing to given values.
	 * 
	 * @param rowSpacing The row spacing ({@link #rowSpacing}) to be set.
	 * @param colSpacing The column spacing ({@link #colSpacing}) to be set.
	 */
	public TwoColumnLayout(final int rowSpacing, final int colSpacing) {
		this.rowSpacing = rowSpacing;
		this.colSpacing = colSpacing;
	}
	
	@Override
	public void addLayoutComponent(String column, Component comp) {
		if (column == null) {
			// Default assume first column.
			column = FIRST_COLUMN;
		}
		// Get the active row. If one does not exist, then create one.
		if (componentRows.isEmpty()) {
			componentRows.add(new Row(null, null));
		}
		Row activeRow  = componentRows.get(componentRows.size() - 1);
		if (FIRST_COLUMN.equals(column)) {
			if (activeRow.column1 != null) {
				// Create a new row because we already have a component
				// in first column in the active row.
				componentRows.add(new Row(comp, null));
			} else {
				activeRow.column1 = comp;
			}
		} else {
			// Second column.
			if (activeRow.column2 != null) {
				// Create a new row because we already have a component
				// in second column in the active row.
				componentRows.add(new Row(null, comp));
			} else {
				activeRow.column2 = comp;
			}
		}
	}

	@Override
	public void removeLayoutComponent(Component comp) {
		for(int rowNum = 0; (rowNum < componentRows.size()); rowNum++) {
			Row r = componentRows.get(rowNum);
			if (r.column1 == comp) {
				r.column1 = null;
			} else if (r.column2 == comp) {
				r.column2 = null;
			}
			if ((r.column1 == null) && (r.column2 == null)) {
				// Empty row. Remove the entry.
				componentRows.remove(rowNum);
				rowNum--; // Account for removal during next iteration.
			}
		}
	}

	/**
	 * Helper method to compute width for first column, second column, and
	 * return total height.
	 * 
	 * This is a refactored utility method that is used to compute
	 * the widths for the first column and second column. It also
	 * handles the situation when a row has only a single component in
	 * the first column in it. The dimensions for each row are
	 * determined via call to {@link Row#getDimensions(int[])} method.
	 *  
	 * @param colWidths An array object that is used to track the 
	 * maximum width of column1, maximum width of column2, and the
	 * minimum column width (for rows with just components in column1) in
	 * index positions 0, 1, and 2 respectively. Note that the 
	 * {@link #colSpacing} is not included in the widths.
	 * 
	 * @return The total height of the layout. The height includes
	 * {@link #rowSpacing} values.
	 */
	private int computeLayoutSize(final int colWidths[]) {
		// Track the required widths for first and second columns and
		// the widest column (for cases where a row has only one component) 
		// in the index positions 0, 1, and 2 respectively.
		int totalHeight = 0;
		// The widest component in each column determines each columns width.
		// The tallest component in each row determines row height.
		for(Row row : componentRows) {
			// Get height of row. If row has all invisible components then
			// row height will be zero and we don't need to account for that row
			final int rowHeight = row.getDimensions(colWidths);
			if (rowHeight > 0) {
				// Track total height and width of each column
				totalHeight += row.getDimensions(colWidths) + rowSpacing;
			}
		}
		// Now account for the widest column for rows that have only 
		// a single component in it.
		if (colWidths[0] + colWidths[1] < colWidths[2]) {
			// The second column is grown to accommodate extra space
			colWidths[1] = colWidths[2] - colWidths[0];
		}
		// Return height of layout.
		return totalHeight - rowSpacing;
	}
	
	@Override
	public Dimension preferredLayoutSize(Container parent) {
		// Track the required widths for first and second columns and
		// the widest column (for cases where a row has only one component) 
		// in the index positions 0, 1, and 2 respectively.
		int colWidths[] = new int[3];
		int totalHeight = computeLayoutSize(colWidths);
		// Account for insets of the parent container
		Insets insets    = parent.getInsets();
		// Finally compute overall dimensions
		final Dimension prefSize = 
			new Dimension(colWidths[0] + colWidths[1] + 
					insets.left + insets.right + colSpacing, 
					totalHeight + insets.top + insets.bottom);
		// Return the overall preferred size.
		return prefSize;
	}

	@Override
	public Dimension minimumLayoutSize(Container parent) {
		return preferredLayoutSize(parent);
	}

	@Override
	public void layoutContainer(Container parent) {
		// Track the required widths for first and second columns and
		// the widest column (for cases where a row has only one component) 
		// in the index positions 0, 1, and 2 respectively.
		int colWidths[] = new int[3];
		computeLayoutSize(colWidths);
		// Next set the bounds of each component in each row.
		final Insets insets = parent.getInsets();
		int top = insets.top, left = insets.left;
		for(Row row: componentRows) {
			// Layout component in first column.
			if (row.column1 != null) {
				final Dimension col1Dim = row.column1.getPreferredSize();
				if (row.column2 == null) {
					// First column can occupy the full row in this case.
					col1Dim.width = colWidths[0] + colWidths[1];
				}
				row.column1.setSize(col1Dim.width, col1Dim.height);
				row.column1.setBounds(left, top, col1Dim.width, col1Dim.height);
			}
			// Layout component in second column.
			if (row.column2 != null) {
				final Dimension col2Dim = row.column2.getPreferredSize();
				row.column2.setSize(col2Dim.width, col2Dim.height);
				row.column2.setBounds(left + colWidths[0] + colSpacing, top, 
						col2Dim.width, col2Dim.height);
			}
			// The tallest of the two columns determines height
			top += row.getDimensions(colWidths) + rowSpacing;
		}
	}

	@Override
	public void addLayoutComponent(Component comp, Object column) {
		if ((column == null) || (column instanceof String)) {
			addLayoutComponent((String) column, comp);
		} else {
			throw new IllegalArgumentException("Cannot add to TwoColumnLayout: " +
					"column must be FIRST_COLUMN or SECOND_COLUMN");
		}
	}

	@Override
	public Dimension maximumLayoutSize(Container parent) {
		return preferredLayoutSize(parent);
	}

	@Override
	public float getLayoutAlignmentX(Container paramContainer) {
		return 0.5f;
	}

	@Override
	public float getLayoutAlignmentY(Container paramContainer) {
		return 0.5f;
	}

	@Override
	public void invalidateLayout(Container paramContainer) {
		// This layout does not cache any information. Consequently,
		// there is nothing to invalidate here.
	}
	
	/**
	 * Get the current row spacing being used by this layout.
	 * 
	 * Row spacing indicates the number of empty pixels between adjacent
	 * rows that are automatically inserted by this layout to improve
	 * presentation.
	 * 
	 * @return The current row spacing used by this layout. The default
	 * value is zero. 
	 */
	public int getRowSpacing() {
		return rowSpacing;
	}


	/**
	 * Set the row spacing to be used by this layout.
	 * 
	 * Row spacing indicates the number of empty pixels between adjacent
	 * rows that are automatically inserted by this layout to improve
	 * presentation.
	 * 
	 * @param rowSpacing the row spacing to be set for this layout.
	 */
	public void setRowSpacing(int rowSpacing) {
		this.rowSpacing = rowSpacing;
	}
	
	/**
	 * Get the current column spacing being used by this layout.
	 * 
	 * Column spacing indicates the number of empty pixels between the
	 * first and second column that are automatically inserted by this 
	 * layout to improve presentation.
	 * 
	 * @return The current column spacing used by this layout. The default
	 * value is zero. 
	 */
	public int getColSpacing() {
		return colSpacing;
	}

	/**
	 * Set the column spacing to be used by this layout.
	 * 
	 * Column spacing indicates the number of empty pixels between
	 * the first and second columns that are automatically inserted 
	 * by this layout to improve presentation.
	 * 
	 * @param colSpacing the column spacing to be set for this layout.
	 */
	public void setColSpacing(int colSpacing) {
		this.colSpacing = colSpacing;
	}
}
