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
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.TreeMap;

import javax.swing.JPanel;
import javax.swing.Timer;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import org.peace_tools.data.EST;
import org.peace_tools.data.ESTEntry;
import org.peace_tools.data.ESTTableModel;
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
 * <p>This class also implements the TableModelListener interface to
 * intercept and process events dispatched by the table in the 
 * {@link HighlightPanel} whenever the user checks/clears an EST
 * entry. This class appropriately updates the list of ESTs to 
 * be highlighted. If there is at least one EST to be highlighted then
 * this class runs a timer to periodically change the background color
 * of highlighted ESTs making them appear to flash.</p> 
 * 
 * <p><b>Note</b>This class is meant to be instantiated and used
 * only by the OverlapView class. Consequently it has been made 
 * package private.</p>
 */
class OverlapPanel extends JPanel implements TableModelListener, ActionListener {
	/**
	 * The only constructor. 
	 * 
	 * @param overlapModel The overlap model from where the data for 
	 * rendering the fragments is obtained.
	 * 
	 * @param colorMapper The helper object to be used to determine the
	 * color for rendering fragments belonging to various clusters.
	 * 
	 * @param findHelper The helper class that tracks the currently active
	 * EST fragment. This information is used to suitably highlight the
	 * currently active fragment. 
	 * 
	 * @param rowScale The initial scale for the rows that determines the
	 * pixel size of each row.
	 *  
	 * @param colScale The initial scale for columns that determines the
	 * pixel size of each column. 
	 */
	public OverlapPanel(OverlapModel overlapModel, ClusterColorMapper colorMapper, 
			FindHelper findHelper, double rowScale, double colScale) {
		super(new BorderLayout(0, 0));
		// Save references to clusterFile for future use
		this.overlapModel = overlapModel;
		this.findHelper   = findHelper;
		this.colorMapper  = colorMapper;
		this.renderInfo   = 1;
		this.setToolTipText("");
		// Save scaling information.
		setScale(rowScale, colScale);
		// Set our default background color to be always white
		setBackground(Color.white);
		// Enable dragging this panel
		setAutoscrolls(true);
		// Create the timer for future reference
		refreshTimer = new Timer(500, this);
		refreshTimer.setActionCommand("refresh");
		// Create the map to maintain list of ESTs to be highlighted
		hilitList = new TreeMap<Integer, Integer>();
		// Finally create the block arrows for highlighting if needed
		if (upDownArrow[0] == null) {
			initializeArrows();
		}
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
		// Render fragments in a row-by-row fashion. While rendering collect
		// information about highlighting arrows to be drawn on top of all the
		// fragments.
		ArrayList<Point> hilitArrowsList = new ArrayList<Point>(10);
		// Add initial blank entries for arrows for the current EST fragment. These
		// will get replaced in the addArrows()  method later on.
		hilitArrowsList.add(null); hilitArrowsList.add(null);
		
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
					paint(g, row, fragment, hilitArrowsList);
				}
			}
		}
		// Draw any highlighting arrows
		paintArrows(g, hilitArrowsList);
	}
	
	/**
	 * This is a helper method that is invoked from the constructor to create two
	 * small arrows. The first arrow points downwards while the second arrow
	 * points upwards. The pair of arrows are used to clearly highlight ESTs
	 * particularly in a large overlap view. Furthermore, the arrows come in
	 * handy when the overlap view is saved as an image.
	 */
	private void initializeArrows() {
		// The following two arrays contain the x and y coordinates of a downward pointing
		// block arrow with its arrow tip at logical (0, 0).
		final int[] xPoints = {0, -4, -2, -2,  2,  2,  4};
		final int[] yPoints = {0, -3, -3, -6, -6, -3, -3};
		// Create the polygon for the downward pointing arrow.
		upDownArrow[0] = new Polygon(xPoints, yPoints, xPoints.length);
		// To create the upward pointing arrow we simply reflect the y-points
		// around the origin.
		for(int i = 0; (i < yPoints.length); i++) {
			yPoints[i] *= -1;
		}
		// Now create the upward pointing arrow
		upDownArrow[1] = new Polygon(xPoints, yPoints, xPoints.length);
	}
	
	/**
	 * Helper method to draw highlighting arrows. 
	 * 
	 * This helper method is called from the {@link #paintComponent(Graphics)} method
	 * to draw any highlighting arrows. This method iterates over the list of entries
	 * in the given arrowList and draws pairs of arrows if the entries are not null.
	 */
	private void paintArrows(Graphics g, ArrayList<Point> arrowsList) {
		for(int i = 0; (i < arrowsList.size()); i += 2) {
			// The first pair of arrows correspond to highlighting the currently selected
			// entry. This arrow is colored in red instead of the usual black.
			g.setColor((i > 0) ? Color.black : Color.red);
			Point topArrow = arrowsList.get(i);
			Point botArrow = arrowsList.get(i+1);
			if (topArrow != null) {
				// If sufficient space is not available for top-arrow then
				// draw a bottom arrow instead.
				if (topArrow.y < ARROW_DIMENSION) {
					paintArrow(g, 1, topArrow.x, botArrow.y);
				} else {
					paintArrow(g, 0, topArrow.x, topArrow.y);
				}
			}
			// Paint bottom arrow
			if (botArrow != null) {
				paintArrow(g, 1, botArrow.x, botArrow.y);
			}
		}
	}
	
	/**
	 * Helper method to actually draw an up or down block arrow at a given
	 * coordinate.
	 * 
	 * @param g The graphics handle to be used to draw the arrow.
	 * @param arrowID The index of the arrow to be drawn. This value is 1 for
	 * up block arrow and 0 for a downwards block arrow.
	 * @param x The x-coordinate where the tip of the arrow needs to be.
	 * @param y The y-coordinate where the tip of the arrow needs to be.
	 */
	private void paintArrow(Graphics g, int arrowID, int x, int y) {
		upDownArrow[arrowID].translate(x, y);
		g.drawPolygon(upDownArrow[arrowID]);
		g.fillPolygon(upDownArrow[arrowID]);
		upDownArrow[arrowID].translate(-x, -y);
	}
	
	/**
	 * Helper method to decide and add appropriate highlighting arrows.
	 * 
	 * This method is invoked from the {@link #paint(Graphics, int, ESTEntry, ArrayList)
	 * method to add an arrow to suitably highlight a given fragment.
	 * 
	 * @param fragX The starting pixel x-coordinate of the fragment.
	 * @param fragY The starting pixel y-coordinate of the fragment.
	 * @param fragWidth The pixel width of the fragment.
	 * @param fragHeight The pixel height to the fragment.
	 * @param hilitArrowsList The list to which the point where the up and down arrows
	 * are to be rendered are to be added.
	 */
	private void addArrows(int fragX, int fragY, int fragWidth, 
			int fragHeight, ArrayList<Point> hilitArrowsList, boolean isCurrEST) {
		Point topArrow = null; // Arrow to draw at the top-left corner
		Point botArrow = null; // Arrow to draw at the bottom-right corner.
		
		// Depending on actual pixel size of the EST we either add just one or two arrows.
		if (fragWidth < 3 * ARROW_DIMENSION) {
			// Space is too small. Draw just a single arrow at the center of the fragment.
			// However, the arrow can be at the top or the bottom.
			final int arrowX = fragX + fragWidth / 2;
			if (fragY < ARROW_DIMENSION) {
				botArrow = new Point(arrowX, fragY + fragHeight + 1);
			} else {
				topArrow = new Point(arrowX, fragY);
			}
		} else {
			// There is sufficient space for two arrows. So add them.
			final int inset = ARROW_DIMENSION / 2;
			topArrow = new Point(fragX + inset, fragY);
			botArrow = new Point(fragX + fragWidth - inset, fragY + fragHeight);
		}
		// Add the pair of arrows to the list at the appropriate location
		if (!isCurrEST) {
			hilitArrowsList.add(topArrow);
			hilitArrowsList.add(botArrow);
		} else {
			hilitArrowsList.set(0, topArrow);
			hilitArrowsList.set(1, botArrow);
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
	 * 
	 *  @param hilitArrowsList This vector is populated with a pair of coordinates
	 *  for each EST to be highlighted with arrows.
	 */
	private void paint(Graphics main, final int row, final ESTEntry fragment, 
					   ArrayList<Point> hilitArrowsList) {
		// Translate the logical model coordinates to actual graphical
		// pixel coordinates and crate a new graphics object to make rest of
		// paint logic easier.
		final EST est    = fragment.getEST();
		final int width  = (int) (est.getSequence().length() * colScale);
		final int height = (int) ((rowScale < 3) ? rowScale : (rowScale - 1)); 
		final int startX = (int) (fragment.getColumn() * colScale);
		final int startY = (int) (row * rowScale);
		// Create a clipped graphics object to make rendering easier.
		Graphics g = main.create(startX, startY, width, height);
		// Now use the sub-graphics object to do rest of painting. First render
		// a rectangle to illustrate the area of the fragment. Here we do a bit
		// of fancy graphics if the rectangle is of sufficient size.
		Color bgColor = colorMapper.getColor(fragment.getClusterID());
		if (bgColor == null) {
			bgColor = Color.lightGray;
		}
		// Adapt the color for flashing highlighted ESTs.
		if (hilitList.containsKey(est.getID())) {
			bgColor = lightOrDark ? bgColor.brighter().brighter() : bgColor.darker().darker();
			// Add arrows to be drawn using helper method. 
			addArrows(startX, startY, width, height, hilitArrowsList, false);
		} else if (findHelper.getCurrESTIndex() == est.getID()) {
			// This is the currently active EST. Add arrow to highlight this entry
			addArrows(startX, startY, width, height, hilitArrowsList, true);
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
	 * This method is invoked whenever the user clicks on a check box
	 * associated with an EST entry in the table housed by the
	 * {@link HighlightPanel}. This method is called because the
	 * {@link OverlapView#OverlapView(org.peace_tools.core.MainFrame, 
	 * OverlapModel, org.peace_tools.data.ClusterFile)} constructor
	 * adds this as the listener to the table in the HighlightPanel.
	 *  
	 * @param e The event containing the information about the row 
	 * in the table that changed.
	 */
	@Override
	public void tableChanged(TableModelEvent e) {
		if (e.getFirstRow() != e.getLastRow()) {
			// The two rows must be same.
			return;
		}
		// Get the current selection status of the EST changed.
		final int estIndex  = e.getFirstRow();
		boolean estSelected = (Boolean) estTableModel.getValueAt(estIndex, 0);
		Integer refCount    = hilitList.get(estIndex);
		if (refCount == null) {
			refCount = 0;
		}
		// Track number of references.
		refCount = refCount + (estSelected ? 1 : -1);
		// Update the list of ESTs to be highlighted.
		if (refCount <= 0) {
			hilitList.remove(estIndex);
		} else {
			hilitList.put(estIndex, refCount);
		}
		// Finally start or stop the timer as needed.
		if (hilitList.size() == 0) {
			// No ESTs to highlight
			refreshTimer.stop();
			// Repaint once more to ensure that highlighted entries
			// are properly repainted
			repaint();
		} else if (!refreshTimer.isRunning()) {
			refreshTimer.start();
		}
		// Select the corresponding entry in the find helper as well.
		ESTEntry entry = findHelper.find(estIndex);
		if ((entry != null) && (estSelected)) {
			// Scroll to make current selection fully visible.
			scrollToCurrent();
		}
	}
	
	/**
	 * Method to scroll the display to the currently selected entry in the
	 * {@link FindHelper} panel.
	 * 
	 * This is a convenience method that can be used to scroll this panel
	 * to display the currently selected entry. This method uses the
	 * information in the {@link #findHelper} object and ensures that the
	 * corresponding fragment is currently visible.
	 * 
	 */
	public void scrollToCurrent() {
		// Get the currently selected fragment.
		final int row    = findHelper.getLastFindRow();
		final int col    = findHelper.getLastFindCol();
		// Get the ESTEntry at the specified row and column from the overlap model.
		final ESTEntry entry = overlapModel.getEntry(row, col);
		final int width  = (int) (entry.getEST().getSequence().length() * colScale);
		final int height = (int) ((rowScale < 3) ? rowScale : (rowScale - 1)); 
		final int startX = (int) (entry.getColumn() * colScale);
		final int startY = (int) (row * rowScale);
		// Dispatch request to scroll to ensure the current selection is visible.
		this.scrollRectToVisible(new Rectangle(startX, startY, width, height));
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if ("refresh".equals(e.getActionCommand())) {
			// Keep toggling the value that is used in the paint() method
			// that renders highlighted ESTs in different colors based on
			// the value of this flag.
			lightOrDark = !lightOrDark;
			repaint();
		}
	}
		
	@Override
	public String getToolTipText(MouseEvent me) {
		// Compute model bounds to render fragments.
		final Rectangle mouseBounds = new Rectangle(me.getX(), me.getY(), 1, 1);
		final Rectangle modelBounds = toModelCoordinates(mouseBounds);
		// Get the fragments in the row using row offset
		final ArrayList<ESTEntry> fragmentList = overlapModel.getRow(modelBounds.y);
		if (fragmentList == null) {
			// No fragments to render on this row.
			return null;
		}
		// Check if some fragment that intersect with the model bounds.
		for(ESTEntry fragment: fragmentList) {
			// check if fragment bounds intersects with our model bounds
			if (modelBounds.intersects(fragment.getColumn(), modelBounds.y, 
						fragment.getEST().getSequence().length(), 1)) {
				// OK, format and return information about this fragment.
				final EST est = fragment.getEST();
				final String info = wrapLine(est.getInfo(), 50);
				final String seq  = wrapLine(est.getSequence(), 60);
				return String.format(TOOL_TIP_FORMAT, est.getID(), info, seq);
			}
		}
		// No fragments under mouse pointer
		return null;
	}
	
	/**
	 * Helper method to wrap long lines to a given size. This method is
	 * used by the {@link #getToolTipText(MouseEvent)} method to format
	 * long headers and nucleotide sequences down to size to make them fit
	 * nicely into a suitable popup box.
	 *  
	 * @param line The line to be wrapped using HTML breaks.
	 * @param width The length of the lines after which the lines are to be wrapped.
	 * @return The line with suitable HTML breaks placed within it.
	 */
	private String wrapLine(final String line, final int width) {
		StringBuilder sb = new StringBuilder(line.length() + 10);
		for(int currPos = 0; (currPos < line.length()); currPos += width) {
			if (currPos > 0) {
				sb.append("<br>");
			}
			sb.append(line.substring(currPos, Math.min(currPos + width, line.length())));
		}
		return sb.toString();
	}
	/**
	 * Method used by {@link OverlapView} to setup the table model to be used
	 * by this class for obtaining EST selection information.
	 * 
	 * @param tableModel The table model to be used for obtaining selection
	 * information.
	 */
	public void setESTTableModel(ESTTableModel tableModel) {
		estTableModel = tableModel;
		estTableModel.addTableModelListener(this);
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
	 * The helper class that tracks the currently active EST fragment
	 * and assists in search operations. This class is used to determine
	 * the currently active entry so that it can be suitably highlighted
	 * in the overlap view.
	 */
	private final FindHelper findHelper;

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
	 * This map is used to maintain the list of ESTs that must be highlighted
	 * by this panel. Entries in this list are managed by the 
	 * {@link #tableChanged(TableModelEvent)} method whenever the user
	 * checks/clears a given EST entry in the {@link HighlightPanel}.
	 */
	private TreeMap<Integer, Integer> hilitList;
	
	/**
	 * This timer is used to periodically schedule a refresh of the current
	 * contents displayed in this overlap panel. The timer is started
	 * whenever at least one EST needs to be highlighted and it stopped
	 * when none of the ESTs need to be highlighted. The timer is started
	 * and stopped in the {@link #tableChanged(TableModelEvent)} method.
	 * The timer triggers calls to the {@link #actionPerformed(ActionEvent)}
	 * method.
	 */
	private Timer refreshTimer;
	
	/**
	 * This flag is used to increase or decrease the background color of an
	 * highlighted EST. This value is toggled each time the refreshTimer
	 * requests a repaint ensuring that the ESTs are prominiently 
	 * highlighted.
	 */
	private boolean lightOrDark;
	
	/**
	 * The table model that actually contains information about the status of
	 * the ESTs being selected (or not selected).
	 */
	private ESTTableModel estTableModel;
	
	/**
	 * The string that is formatted to obtain a HTML style tool tip for display
	 * when the user hovers over an EST entry in the overlap panel.
	 */
	private static final String TOOL_TIP_FORMAT = "<html>" +	
		"<table>" +
			"<tr><td><b>ID/Index:</b></td><td>%d</td></tr>" + 
			"<tr><td valign=\"top\"><b>Information:</b></td><td>%s</td></tr>" +
			"<tr><td valign=\"top\"><b>Sequence:</b></td><td><font size=\"-3\" face=\"monospace\">%s</font></td></tr>" +
		"</table>" + 
		"</html>";
	
	/**
	 * Polygons that are used to draw small arrows around highlighted ESTs to make
	 * the highlighting show up a bit better when the overlap view is saved as 
	 * an image. The arrows are created by the {@link #initializeArrows()} method
	 * (that is called from the constructor) and used by the
	 * {@link #paintAll(Graphics)} method (called from the {@link #paintComponent(Graphics)}
	 * method). The first entry in the array contains a downward pointing block arrow and
	 * the second entry contains an upward pointing block arrow. 
	 */
	private static final Polygon upDownArrow[] = new Polygon[2];
	
	/**
	 * The width or height of the up/down arrows used for highlighting. Ideally
	 * this value would be calculated based on the size of the arrows. However,
	 * for convenience, this constant is being defined.
	 */
	private static final int ARROW_DIMENSION = 6;
	
	/**
	 * Generated serialization UID (merely to keep the compiler happy)
	 */
	private static final long serialVersionUID = 5321838842105746304L;
}
