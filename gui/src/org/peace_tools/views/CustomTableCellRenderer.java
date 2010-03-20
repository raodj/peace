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

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.Graphics;

import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

/**
 * A simple renderer that overrides font in default renderer.
 * 
 * This is an custom class that is used by the tree or table views
 * to suitably update the font size and font type for the default
 * cell renderer. The same column renderer can be used for all objects
 * and the font is also reused.
 * 
 * @see org.peace_tools.data.ESTTableModel
 * @see ClusterTreeTableView
 */
public class CustomTableCellRenderer extends DefaultTableCellRenderer {
	/**
	 * The mono-spaced font used to render the EST sequences in
	 * the tree view.
	 */
	private Font monoFont = null;

	/**
	 * <p>Flag to indicate if base pairs in sequences must be highlighted.
	 * This feature increases time for rendering data but can be useful
	 * to luckily identify patterns.</p>
	 * 
	 * <p><b>Note</b>: Currently this feature does not operate correctly.</p> 
	 */
	private boolean highlightForeground = false;
	
	/**
	 * Flag to indicate if background for the  base pairs in sequences
	 * must be colored differently for different base pairs.
	 * This feature increases time for rendering data but can be useful
	 * to quickly identify patterns.
	 * 
	 * <p><b>Note</b>: Currently this feature does not operate correctly.</p>
	 */
	private boolean highlightBackground = false;
	
	/**
	 * The constructor merely initializes the mono space font used
	 * to render EST sequences.
	 */
	public CustomTableCellRenderer() {
		Font defaultFont = super.getFont();
		monoFont = new Font(Font.MONOSPACED, defaultFont.getStyle(),
				defaultFont.getSize() - 2);
	}

	@Override
	public Component getTableCellRendererComponent(JTable table,
			Object value, boolean isSelected, boolean hasFocus, int row,
			int column) {
		// Let the user class set up the renderer with defaults
		super.getTableCellRendererComponent(table, value, isSelected, hasFocus, 
				row, column);
		// Set the font to be used to render
		setFont(monoFont);
		return this;
	}
	
	@Override
	public void paint(Graphics g) {
		if (!highlightForeground) {
			// No highlighting to be done.
			super.paint(g);
			return;
		}
		// Draw the string with highlighting enabled.
		// First figure out the width of one character. Recollect we are dealing
		// with mono spaced font here.
		final int charWidth = g.getFontMetrics().charWidth('0');
		// Paint the text while highlighting the base pair strings encountered.
		final char[] seq = getText().toLowerCase().toCharArray();
		// Draw the characters one at a time.
		int currX          = 0;
		final int minX     = g.getClipBounds().x;
		final int maxX     = g.getClipBounds().x + g.getClipBounds().width;
		final int height   = g.getClipBounds().height;
		final int startY   = g.getClipBounds().y;
		final String Bases = "atcg";
		// Draw background for each character if background highlighting is requested.
		if (this.highlightBackground) {
			for(int i = 0; (i < seq.length); i++) {
				if (currX > maxX) {
					break;
				}
				if (currX + charWidth > minX) {
					int colorCode = Bases.indexOf(seq[i]);
					if (colorCode != -1) {
						g.setColor(BASE_FORE_COLORS[colorCode]);
					} else {
						g.setColor(getBackground());
					}
					g.fillRect(currX, startY, charWidth, height);
				}
				// Onto the next character space
				currX += charWidth;
			}
		}
		// Draw the base pairs as requested.
		final int textY = g.getClipBounds().y + g.getFontMetrics().getHeight();
		g.setColor(getForeground());
		g.drawString(getText(), minX, textY);
	}

	/**
	 * The list of colors to be used to highlight the base pairs. These colors
	 * are primarily to be used as the foreground colors. 
	 */
	private static final Color BASE_FORE_COLORS[] = {
			Color.red, Color.blue, Color.green, Color.MAGENTA
	};
	
	/**
	 * The general serial version UID to enable serialization of this class as
	 * per Java requirements.
	 */
	private static final long serialVersionUID = -9175031053323866144L;
}
