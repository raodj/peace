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

package org.peace_tools.generic;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Insets;
import javax.swing.border.AbstractBorder;

/**
 * A custom border that draws 4 lines along the edges of a given
 * component based on either color codes or four custom colors.
 * Each side of the border can be independently controlled.
 */
public class CustomBorder extends AbstractBorder {
	/**
	 * Create a custom border using 4 character codes for
	 * top, left, bottom, right. The character codes are
	 * "d" for darker, 'D' for much-darker, "l" for lighter, 
	 * 'L' for much-lighter, and "s" for same. Any other 
	 * characters result in a null color and a border is not
	 * drawn. 
	 *  
	 * @param colorCodes A string with 4 characters (example:
	 * "ddll" indicating the border colors in the order top, 
	 * left, bottom, right.
	 */
	public CustomBorder(String colorCodes) {
		this.colorCodes = colorCodes;
		this.colors     = null;
	}
	
	/**
	 * Create a custom border with a given set of colors. If one of the
	 * colors is null, then a border is not drawn on that side and the
	 * border does not take up space on that side.
	 * 
	 * @param top Color to be used for the top line.
	 * @param left Color to be used for the left line.
	 * @param bottom Color to be used for the bottom line.
	 * @param right Color to be used for the right line.
	 */
	public CustomBorder(Color top, Color left, Color bottom, Color right) {
		colorCodes = null;
		colors     = new Color[4];
		colors[0]  = top;
		colors[1]  = left;
		colors[2]  = bottom;
		colors[3]  = right;
	}
	
	/**
	 * Paints the border for the specified component with the
	 * specified position and size.
	 * 
	 * @param comp the component for which this border is being painted
	 * @param g the paint graphics
	 * @param x the x position of the painted border
	 * @param y the y position of the painted border
	 * @param width the width of the painted border
	 * @param height the height of the painted border
	 */
	public void paintBorder(Component comp, Graphics g, int x, int y, 
			int width, int height) {
		int w   = width;
		int h   = height;
		Color c = null;
		g.translate(x, y);

		// draw top line.
		if ((c = getColor(0, comp)) != null) {
			g.setColor(c);
			g.drawLine(0, 0, w-1, 0);
		}
		// draw left line.
		if ((c = getColor(1, comp)) != null) {
			g.setColor(c);
			g.drawLine(0, 0, 0, h-1);
		}
		// draw bottom line.
		if ((c = getColor(2, comp)) != null) {
			g.setColor(c);
			g.drawLine(0, h-1, w-1, h-1);
		}
		// draw right line.
		if ((c = getColor(3, comp)) != null) {
			g.setColor(c);
			g.drawLine(w-1, 0, w-1, h-1);
		}
	
		g.translate(-x, -y);
	}

	/**
	 * Returns the insets of the border.
	 * @param c the component for which this border insets value applies
	 */
	public Insets getBorderInsets(Component c) {
		return new Insets(getColor(0, c) != null ? 1 : 0, 
				getColor(1, c) != null ? 1 : 0,
				getColor(2, c) != null ? 1 : 0,
				getColor(3, c) != null ? 1 : 0);
	}

	/**
	 * Reinitialize the insets parameter with this Border's current Insets.
	 * @param c the component for which this border insets value applies
	 * @param insets the object to be reinitialized
	 */
	public Insets getBorderInsets(Component c, Insets insets) {
		insets.top    = getColor(0, c) != null ? 1 : 0;
		insets.left   = getColor(1, c) != null ? 1 : 0;
		insets.bottom = getColor(2, c) != null ? 1 : 0;
		insets.right  = getColor(3, c) != null ? 1 : 0;
		return insets;
	}

	/**
	 * Returns whether or not the border is opaque.
	 * 
	 * @return This method always returns true to indicate that this
	 * border is opaque.
	 */
	public boolean isBorderOpaque() { return true; }

	/**
	 * Helper method to get the color associated with a given side
	 * of the border. The numbers 0, 1, 2, and 3 are top, left, bottom,
	 * and right.
	 * 
	 * @param i The index/direction of border.
	 * @param c The component from which colors are derived.
	 * @return A color code if a border is to be drawn. null otherwise.
	 */
	protected Color getColor(int i, Component c) {
		if (colors != null) {
			return colors[i];
		}
		Color base = c.getBackground();
		switch(colorCodes.charAt(i)) {
		case 'd': return base.darker();
		case 'l': return base.brighter();
		case 'D': return base.darker().darker().darker();
		case 'L': return base.brighter().brighter().brighter();
		case 's': return base;
		}
		// No color for this edge.
		return null;
	}
	

	/**
	 * A generated serialization values. 
	 */
	private static final long serialVersionUID = 3100922135406430414L;

	/**
	 * Create a custom border using 4 character codes for
	 * top, left, bottom, right. The character codes are
	 * "d" for darker and "l" for lighter. Any other 
	 * characters result in a null color and a border is not
	 * drawn. 
	 */
	String colorCodes;
	
	/**
	 * An array of four colors corresponding to top, left, bottom,
	 * right borders
	 */
	Color colors[];
}
