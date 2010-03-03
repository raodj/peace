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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Rectangle;

import javax.swing.JComponent;

import org.peace_tools.data.OverlapModel;

/**
 * The overlap panel that displays the actual overlap image.
 *
 * <p>This panel performs the core task of drawing the various fragments
 * illustrating the overlap. The fragments are drawn row-by-row using
 * the data from a given (specified in the constructor) OverlapModel
 * object. The size and location of the rectangles for each fragment
 * depend on the current scaling factors set for this component (via 
 * call to setScale() method). The scaling factors also control the 
 * preferred size reported by this panel. The color map maintained by
 * this class is used to color the fragments based on their cluster ID.</p>   
 *
 * <p><b>Note</b>This class is meant to be instantiated and used
 * only by the OverlapView class. Consequently it has been made 
 * package private.</p>
 */
class Ruler extends JComponent {
	/**
	 * Enumerations to define the two types of rulers that this class can
	 * be used for. Vertical rulers track the rows in the model while 
	 * horizontal rulers track the columns in the model. 
	 */
	public enum Orientation { HORIZONTAL, VERTICAL };
	
	/**
	 * The only constructor. 
	 * 
	 * @param modelWidth The total width of the model in nucleotides. This
	 * value is typically the {@link OverlapModel#getMaxCol()}.
	 * 
	 * @param modelHeight The total height of the model in number of rows
	 * of nucleotides. This value is typically {@link OverlapModel#getMaxRow()}.
	 * 
	 * @param orientation The orientation that defines the type of ruler that
	 * this instance is to represent. The orientation determines the overall 
	 * behavior of this class.
	 */
	public Ruler(int modelWidth, int modelHeight, Orientation orientation) {
		this.modelSize   = new Dimension(modelWidth, modelHeight);
		this.orientation = orientation;
	}
	
	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		// Paint the component.
		if (this.orientation == Orientation.HORIZONTAL) {
			paintHorizRule(g);
		} else {
			paintVerticalRule(g);
		}			
	}
	
	/**
	 * Helper method to paint the horizontal ruler. 
	 * 
	 * This method is invoked from the main paint component method in order
	 * to draw a horizontal rule with the necessary ticks and labels. This
	 * method assumes that this ruler will always be added to the top
	 * of the overlap view panel.
	 * 
	 * @param g The graphics object to be used for drawing the ticks.
	 */
	private void paintHorizRule(Graphics g) {
		final int MIN_TICK_SPACING = 50;
		// Compute pixels to skip along x-axis. Minimum pixel skip is 50.
		// However, if columns occupy more than 1 pixel try to center ticks
		// on those columns as much as possible.
		final int tickSpacing = (colScale <= 1) ? MIN_TICK_SPACING 
				: (int) (((int) (MIN_TICK_SPACING / colScale) * colScale));
		// Total graphical size of this ruler (in pixels)
		final Dimension mySize = getPreferredSize();
		// Draw ticks in the given clip bounds
		final Rectangle bounds = g.getClipBounds();
		// Compute the bounds of iteration that we need to do.
		final int startTickPos = (bounds.x - (bounds.x % tickSpacing));
		final int endTickPos   = startTickPos + bounds.width + tickSpacing;
		// Setup small font for the tick labels.
		g.setFont(g.getFont().deriveFont(8.0f));
		final FontMetrics fm = g.getFontMetrics();
		// Iterate over bounds and draw ticks and labels
		g.setColor(Color.black);
		// Draw the ticks. Here we use long because when drawing wide rules the 
		// tickValue comes out negative due to int overflow.
		for(int tickPos = startTickPos; (tickPos < endTickPos); tickPos += tickSpacing) {
			// Compute the model value for the tick, given its pixel tick position
			final int tickValue = (int) ((long) tickPos * modelSize.width / mySize.width);
			// Draw a simple black line for the tick
			g.drawLine(tickPos, mySize.height - 5, tickPos, mySize.height);
			// Draw the value centered around the given tickPosition.
			String label   = Integer.toString(tickValue);
			Rectangle labelBounds = fm.getStringBounds(label, g).getBounds();
			// Compute coordinates to center the string.
			int labelX = tickPos - labelBounds.width / 2;
			int labelY = mySize.height - 6 - fm.getDescent();
			// Draw the label.
			g.drawString(label, labelX, labelY);
		}
		// Finally draw a darker line at the bottom to make the ruler
		// look pretty.
		g.setColor(getBackground().darker());
		g.drawLine(bounds.x, RULE_HEIGHT - 1, bounds.x + bounds.width, RULE_HEIGHT - 1);
	}
	
	/**
	 * Helper method to paint the horizontal ruler. 
	 * 
	 * This method is invoked from the main paint component method in order
	 * to draw a horizontal rule with the necessary ticks and labels. This
	 * method assumes that the vertical ruler will always be displayed to
	 * the left of the overlap panel.
	 * 
	 * @param g The graphics object to be used for drawing the ticks.
	 */
	private void paintVerticalRule(Graphics g) {
		final int MIN_TICK_SPACING = 50;
		// Compute pixels to skip along y-axis. Minimum pixel skip is 50.
		// However, if rows occupy more than 1 pixel try to center ticks
		// on those rows as much as possible.
		final int tickSpacing = (rowScale <= 1) ? MIN_TICK_SPACING 
				: (int) (((int) (MIN_TICK_SPACING / rowScale) * rowScale));
		// Total graphical size of this ruler (in pixels)
		final Dimension mySize = getPreferredSize();
		// Draw ticks in the given clip bounds
		final Rectangle bounds = g.getClipBounds();
		// Compute the bounds of iteration that we need to do.
		final int startTickPos = (bounds.y - (bounds.y % tickSpacing));
		final int endTickPos   = startTickPos + bounds.height + tickSpacing;
		// Setup small font for the tick labels.
		g.setFont(g.getFont().deriveFont(8.0f));
		final FontMetrics fm = g.getFontMetrics();
		// Iterate over bounds and draw ticks and labels
		g.setColor(Color.black);
		// Draw the ticks. 
		for(int tickPos = startTickPos; (tickPos < endTickPos); tickPos += tickSpacing) {
			// Compute the model value for the tick, given its pixel tick position
			// Here we use long because when drawing tall rules the tickValue 
			// comes out negative due to integer overflow.
			final int tickValue = (int) ((long) tickPos * modelSize.height / mySize.height);
			// Draw a simple black line for the tick
			g.drawLine(mySize.width - 5, tickPos, mySize.width, tickPos);
			// Draw the value right-justified to the left of the tick.
			String label   = Integer.toString(tickValue);
			Rectangle labelBounds = fm.getStringBounds(label, g).getBounds();
			// Compute coordinates to center the string around tick
			int labelX = mySize.width - 6 - labelBounds.width;
			int labelY = tickPos + labelBounds.height / 2;
			// Draw the label.
			g.drawString(label, labelX, labelY);
		}
		// Finally draw a darker line to the right to make the ruler
		// look pretty.
		g.setColor(getBackground().darker());
		g.drawLine(RULE_WIDTH - 1, bounds.y, RULE_WIDTH - 1, bounds.y + bounds.height);
	}

	/**
	 * Method to set the current scales used by this component.
	 * 
	 * <p>This method must be used to set the scales to be used for painting
	 * rows and columns in this panel. The scales also determine the 
	 * overall size of this component.</p> 
	 * 
	 * <b>Note</b>: The containing component my repaint itself to ensure that
	 * changes in bounds are correctly refreshed on screen.
	 * 
	 * @param rowScale The scale that determines the number of pixels occupied
	 * by each row. This value cannot be zero. However, no checks are performed
	 * on the scale. Consequently, the behavior of this component with a zero
	 * scale is undefined.
	 * 
	 * @param colScale The scale that determines the number of pixels occupied
	 * by each row. This value cannot be zero. However, no checks are performed
	 * on the scale. Consequently, the behavior of this component with a zero
	 * scale is undefined.
	 */
	public void setScale(final double rowScale, final double colScale) {
		this.rowScale = rowScale;
		this.colScale = colScale;
		// Ensure overall size is set to reflect the scale change.
		Dimension prefSize = 
			new Dimension((int) (colScale * modelSize.width),
					(int) (rowScale * modelSize.height));
		// Update the dimension depending on the orientation
		prefSize.width  = (orientation == Orientation.HORIZONTAL) ? prefSize.width  : RULE_WIDTH;
		prefSize.height = (orientation == Orientation.VERTICAL)   ? prefSize.height : RULE_HEIGHT;
		// Set the preferred size for this component.
		this.setPreferredSize(prefSize);
		// Finally re-layout the panel (and its parent).
		invalidate();
	}
	
	/**
	 * The scale to be applies when rendering a single row. This value
	 * represents the maximum height of each row.  Once the row height is
	 * above 3 at least one blank line is allocated for each row to clearly
	 * delineate one row from the other. This value is initialized in the 
	 * constructor and can be changed via a call to the setScale() method 
	 * in this class. Typically this value is in integral units without any
	 * fractions. However, it is maintained as a double to accommodate 
	 * non-integral scales in the future. However, this value should never
	 * be set to zero.  
	 */
	private double rowScale;
	
	/**
	 * The scale to be applies when rendering a single column. This value
	 * represents the maximum width of each column.  This value is 
	 * initialized in the constructor and can be changed via a call to
	 * the setScale() method in this class. Typically this value is in
	 * integral units when columns occupy more than 1 pixel. If a column
	 * occupies fewer than one pixel then this value is a fractional value.
	 * This value should never be set to zero.
	 */
	private double colScale;
	
	/**
	 * The Dimension that provides the necessary information
	 * for rendering the ruler displayed in this class. Horizontal
	 * ruler uses the modelWidth while the vertical ruler uses 
	 * modelHeight value from this instance variable.
	 * This instance variable is initialized by the constructor and
	 * is never changed during the life time of this object.
	 */
	private final Dimension modelSize;

	/**
	 * The orientation being used for this ruler. This value is set
	 * in the constructor and is never changed for this model.
	 */
	private final Orientation orientation;
	
	/**
	 * The preferred width of the vertical ruler. This width is defined
	 * as a constant so that it is easy to set the size of the ruler
	 * appropriately.
	 */
	private static int RULE_WIDTH = 50;
	
	/**
	 * The preferred height of the horizontal ruler. This height is defined
	 * as a constant so that it is easy to set the size of the ruler
	 * appropriately.
	 */
	private int RULE_HEIGHT = 23;
	
	/**
	 * A generated serialization UID just to keep the compiler happy.
	 */
	private static final long serialVersionUID = -292125515405864885L;
}
