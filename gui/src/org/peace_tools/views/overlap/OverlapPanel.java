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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.util.ArrayList;

import javax.swing.JPanel;

import org.peace_tools.data.EST;
import org.peace_tools.data.ESTEntry;
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
class OverlapPanel extends JPanel {
	/**
	 * The only constructor. 
	 * 
	 * @param overlapModel The overlap model from where the data for 
	 * rendering the fragments is obtained.
	 * 
	 * @param colorMapper The helper object to be used to determine the
	 * color for rendering fragments belonging to various clusters.
	 * 
	 *  @param rowScale The initial scale for the rows that determines the
	 *  pixel size of each row.
	 *  
	 *  @param colScale The initial scale for columns that determines the
	 *  pixel size of each column. 
	 */
	public OverlapPanel(OverlapModel overlapModel, ClusterColorMapper colorMapper, 
			double rowScale, double colScale) {
		super(new BorderLayout(0, 0));
		// Save references to clusterFile for future use
		this.overlapModel = overlapModel;
		this.colorMapper  = colorMapper;
		this.renderInfo   = 1;
		// Save scaling information.
		setScale(rowScale, colScale);
		// Set our default background color to be always white
		setBackground(Color.white);
	}
	
	/**
	 * Set the type of information to be rendered for each fragment.
	 * The valid values that can passed to this method are:
	 * 0, 1, 2, or 3. These values are interpreted as follows:
	 * 
	 * <ol start=0>
	 * <li>Do not render any details (just draws a rectangle).</li>
	 * <li>Show nucleotide sequence information.</li>
	 * <li>Show FASTA identifier/header information.</li>
	 * <li>Show debugging information.</li>
	 * </ol>
	 * 
	 * @param renderInfo A simple integer code indicating the type of 
	 * information to be rendered by this panel. See list above for 
	 * details.
	 */
	public void setRenderInfo(int renderInfo) {
		this.renderInfo = renderInfo;
	}

	/**
	 * Obtain the type of information being rendered for each fragment.
	 * 
	 * @return Returns a simple integer code indicating the type of
	 * information being rendered by this panel. The value values returned
	 * by this method are 0, 1, 2, or 3 and they are to be interpreted
	 * as follows:
	 * 
	 * <ol start=0>
	 * <li>Not rendering any details (just draws a rectangle).</li>
	 * <li>Showing nucleotide sequence information.</li>
	 * <li>Showing FASTA identifier/header information.</li>
	 * <li>Showing debugging information.</li>
	 * </ol>
	 */
	public int getRenderInfo() {
		return renderInfo;
	}
	
	/**
	 * Obtain the current row scaling factor being used.
	 * 
	 * @return The current row scaling factor being used. This value is
	 * integral values.
	 */
	public double getRowScale() {
		return rowScale;
	}

	/**
	 * Obtain the current column scaling factor being used.
	 * 
	 * @return The current column scaling factor being used. This value can
	 * be fractional values that may have to be normalized to display or 
	 * work within the GUI for sliders.
	 */
	public double getColScale() {
		return colScale;
	}
	
	/**
	 * Helper method to translate graphics coordinates to model coordinates.
	 * 
	 * This method is a helper method that is used when painting
	 * fragments. This method is used to translate the graphics coordinates
	 * (in pixel values) to model coordinates (rows and columns in the 
	 * OverlapModel). This method uses the rowScale and colScale values
	 * (set via call to {@link OverlapPanel#setScale(double, double)} method)
	 * to perform the translation. 
	 * 
	 * @param pixelBounds The graphical coordinates (in pixel values) to be 
	 * translated.
	 * 
	 * @return A Rectangle object containing the graphical coordinates translated
	 * to model coordinates.
	 */
	private Rectangle toModelCoordinates(Rectangle pixelBounds) {
		// Compute the area to be painted in row, column
		// values instead of the raw pixel values based on 
		// rowScale and colScale values.
		Rectangle normBounds  = new Rectangle((int)(pixelBounds.x / colScale), 
				(int) (pixelBounds.y / rowScale), 
				(int) (pixelBounds.width / colScale)  + 1,
				(int) (pixelBounds.height / rowScale) + 1);
		// Return the translated value back for further use.
		return normBounds;
	}
	
	/**
	 * Override the default paint method for this panel to render fragments.
	 * 
	 * This method overrides the default paint method in the parent class
	 * to render fragments using the following approach:
	 * 
	 * <ol>
	 * 
	 * <li>The clip bounds of the graphics object is translated to model
	 * coordinates using the toModelCoordinates method.</li>
	 * 
	 * <li>For each row to be rendered, this method obtains the list of fragments
	 * in each row and checks to see if the fragment is to be rendered.</li>
	 *  
	 * </ol>
	 * 
	 * @param g The graphics component to be used for rendering the 
	 * fragments.
	 */
	public void paintComponent(Graphics g) {
		// Let the parent class clear out background
		super.paintComponent(g);
		// Compute model bounds to render fragments.
		final Rectangle modelBounds = toModelCoordinates(g.getClipBounds());
		// System.out.println("Painting " + modelBounds + " within " + g.getClipBounds());
		// Render fragments in a row-by-row fashion.
		for(int rowOffset = 0; (rowOffset < modelBounds.height); rowOffset++) {
			// Track the current row we are working on.
			final int row = modelBounds.y + rowOffset;
			// Get the fragments in the row using row offset
			final ArrayList<ESTEntry> fragmentList = 
				overlapModel.getRow(row);
			if (fragmentList == null) {
				// No fragments to render on this row.
				continue;
			}
			// Check and render fragments that intersect with the model bounds.
			for(ESTEntry fragment: fragmentList) {
				// check if fragment bounds intersects with our model bounds
				if (modelBounds.intersects(fragment.getColumn(), row, 
						fragment.getEST().getSequence().length(), 1)) {
					// OK, render this fragment using helper method.
					paint(g, row, fragment);
				}
			}
		}
	}
	
	/**
	 * Helper method to paint a given fragment using a given graphics handle.
	 * 
	 * This method operates as follows:
	 * 
	 * <ol>
	 * 
	 * <li>It computes the bounds of the fragment in terms of pixel coordinates
	 * given the fragments row, starting column, fragment sequence length, and
	 * the rowScale.</li>
	 * 
	 * <li>Using the cluster color set for the fragment's cluster, it draws a 
	 * rectangle of the appropriate size.</li>
	 * 
	 * <li>It then renders the nucleotide sequence or FASTA identifier information
	 * depending on the configuration.</li>
	 * 
	 * </ol>
	 * 
	 * @param main The graphics handle to be used for rendering the given 
	 * fragment.
	 * 
	 * @param row The logical model row where the fragment is to be rendered.
	 * 
	 * @param fragment The fragment to be rendered. 
	 */
	private void paint(Graphics main, final int row, final ESTEntry fragment) {
		// Translate the logical model coordinates to actual graphical
		// pixel coordinates and crate a new graphics object to make rest of
		// paint logic easier.
		final EST est    = fragment.getEST();
		final int width  = (int) (est.getSequence().length() * colScale);
		final int height = (int) ((rowScale < 3) ? rowScale : (rowScale - 1)); 
		Graphics g = main.create((int) (fragment.getColumn() * colScale),
				(int) (row * rowScale),    // y-coordinate
				width, height);
		// Now use the sub-graphics object to do rest of painting. First render
		// a rectangle to illustrate the area of the fragment. Here we do a bit
		// of fancy graphics if the rectangle is of sufficient size.
		Color bgColor = colorMapper.getColor(fragment.getClusterID());
		if (bgColor == null) {
			bgColor = Color.lightGray;
		}
		g.setColor(bgColor);
		if ((width > 4) && (height > 4)) {
			// Draw a bit of a fancy rectangle.
			g.fill3DRect(0, 0, width, height, true);
		} else {
			// Insufficient size. Do a simple rectangle.
			g.fillRect(0, 0, width, height);
			// No need to render additional information as the space is
			// too small to be meaningful.
			return;
		}
		// When control drops here we have sufficient space to render.
		// Next compute the information that has been requested to be rendered
		// in the rectangles.
		String info = "";
		switch (this.renderInfo) {
		case 1: info = est.getSequence(); break;
		case 2: info = est.getInfo(); break;
		case 3: info = String.format("Row: %d, Col: %d, len: %d, cluster: %d", 
				row, fragment.getColumn(), est.getSequence().length(), 
				fragment.getClusterID());
		default:
		}
		// Ensure that the information that we render never exceeds the
		// space allocated for this fragment.
		if (info.length() > est.getSequence().length()) {
			info = info.substring(0, est.getSequence().length());
		}
		// Set the color and render characters.
		// Determine font size if needed.
		if (font == null) {
			font = this.computeFont(g);
		}
		g.setFont(font);
		// Use black for brighter colors. Otherwise use white.
		if (Math.max(Math.max(bgColor.getRed(), bgColor.getGreen()), bgColor.getBlue()) > 100) {
			g.setColor(Color.black);
		} else {
			g.setColor(Color.white);
		}
		// Get generic font height and ascent to center characters vertically
		final int baseLine = g.getFontMetrics().getAscent() + 
			((int) (rowScale - g.getFontMetrics().getHeight()) / 2) - 1;
		for(int i = 0, xPos = 0; (i < info.length()); i++, xPos += colScale) {
			// Draw a character at the appropriate column and offset
			g.drawString(info.substring(i, i + 1), xPos, baseLine);
		}
	}
	
	/**
	 * Compute the font that will fit letters within each fragment rectangle.
	 * 
	 * This method is invoked from the {@link #paintComponent(Graphics)} method
	 * whenever the {@link #font} instance variable is null. This method uses
	 * the supplied Graphics object to compute the largest mono-spaced font that
	 * will fit within the allocated space for a fragment.
	 * 
	 * @param g The graphics handle to be used for computing the font sizes.
	 * 
	 * @return The font object to be used for rendering characters within the
	 * boundaries of each fragment.
	 */
	private Font computeFont(Graphics g) {
		// Start from a  reasonably small size.
		int size = (int) Math.min(rowScale, colScale) / 2;
		// Grow size while ensuring the font fits within the allocated space.
		// But don't let size get too big.
		for(; (size < 14); size++) {
			Font font      = new Font(Font.MONOSPACED, Font.PLAIN, size);
			FontMetrics fm = g.getFontMetrics(font);
			if ((fm.getHeight() >= rowScale - 1) || (fm.charWidth('*') >= colScale)) {
				// This size is already too big. Use previous size.
				break;
			}
		}
		// When control drops here, the size variable refers to the font that
		// is 1 point bigger than what we care for...
		return new Font(Font.MONOSPACED, Font.PLAIN, size - 1);
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
			new Dimension((int) (colScale * overlapModel.getMaxCol()),
					(int) (rowScale * overlapModel.getMaxRow()));
		this.setPreferredSize(prefSize);
		// Clear out the font to force recompute during next painting
		font = null;
		// Finally re-layout the panel (and its parent).
		invalidate();
	}
	
	/**
	 * The scale to be applies when rendering a single row. This value
	 * represents the maximum height of each row.  Once the row height is
	 * above 3 at least one blank line is allocated for each row to clearly
	 * delineate one row from the other. This value is 
	 * initialized in the constructor and can be changed via a call to
	 * the setScale() method in this class. Typically this value is in
	 * integral units without any fractions. However, it is maintained as
	 * a double to accommodate non-integral scales in the future. However,
	 * this value should never be set to zero.  
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
	 * The OverlapModel that provides the necessary information
	 * for rendering the core information required by this panel.
	 * This instance variable is initialized by the constructor and
	 * is never changed during the life time of this object.
	 */
	private final OverlapModel overlapModel;
	
	/**
	 * Helper to look-up color for a given cluster ID.
	 * 
	 * This interface is used when rendering fragments. The color mapper
	 * is used to look up the colors set for a given cluster, given its
	 * cluster ID. The cluster colors are set and updated indirectly by
	 * the OverlapView class.
	 */
	private final ClusterColorMapper colorMapper;
	
	/**
	 * Variable to indicate the information about a fragment that should
	 * be rendered. The value of this instance variable is either 0, 1, 2, or 3.
	 * These values are interpreted as follows:
	 * 
	 * <ol start=0>
	 * <li>Do not render any details (just draws a rectangle).</li>
	 * <li>Show nucleotide sequence information.</li>
	 * <li>Show FASTA identifier/header information.</li>
	 * <li>Show debugging information.</li>
	 * </ol>
	 */
	private int renderInfo;
	
	/**
	 * The font being used by this panel to render information for 
	 * each fragment. This font is computed based on the scale set for
	 * this panel. The font is computed (in the {@link #computeFont(Graphics)}
	 * method) such that the characters fit within the scales specified 
	 * for this object.  Each time the scale is changed the font is reset to
	 * null (in the {@link #setScale(double, double)} method) to fore 
	 * computing a new font.
	 */
	private Font font;
	
	/**
	 * Generated serialization UID (merely to keep the compiler happy)
	 */
	private static final long serialVersionUID = 5321838842105746304L;
}
