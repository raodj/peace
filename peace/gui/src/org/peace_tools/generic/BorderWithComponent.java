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

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JComponent;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;

/**
 * A border that displays a component in the border at the top-left corner.
 * 
 *  This is a custom border component that can be used to display any
 *  swing component as part of the border. This class is based on
 *  the example shown at http://www.javalobby.org/java/forums/t33048.html.
 * 
 */
public class BorderWithComponent implements Border, MouseListener {
	/**
	 * The component to be displayed in the border. This value cannot be
	 * null.
	 */
	private final Component comp;

	/**
	 * The rectangle that identifies the bounds of the {@link #comp} component
	 * present in the border. This value is setup in the {@link #paintBorder(Component, Graphics, int, int, int, int)}
	 * method and is used to dispatch mouse events in {@link #dispatchEvent(MouseEvent)}
	 * method in this class.
	 */
	private Rectangle compBounds;
	
	/**
	 * The parent component that actually owns this border. This value
	 * is set in the constructor and is never changed during the life-time
	 * of this component. This value cannot be null.
	 */
	private final JComponent container;

	/**
	 * The border to be actually drawn by this component. This value
	 * cannot be null.
	 */
	private final Border border;

	/**
	 * The pixel offset to the left to indent the component a bit on
	 * the border.
	 */
	private static final int LEFT_OFFSET = 5;
	
	/**
	 * The constructor for creating a border with a component.
	 * 
	 * @param comp The component to be drawn in the border. This component
	 * cannot be null.
	 * 
	 * @param container The container that logically owns this border component.
	 * 
	 * @param border The actual border to be drawn around this component.
	 */
	public BorderWithComponent(Component comp, JComponent container, Border border) { 
		this.comp      = comp; 
		this.container = container; 
		this.border    = border; 
		container.addMouseListener(this); 
	} 

	/**
	 * Method to draw this border component.
	 * 
	 * @param c The component around which the border is to be drawn.
	 * 
	 * @param g The graphics object to be used to draw the actual border.
	 * 
	 *  @param x The x-coordinate (relative to component c) where the
	 *  border is to be drawn.
	 *  
	 *  @param y The y-coordinate (relative to component c) where the 
	 *  border is to be drawn.
	 *  
	 *  @param width The width of the component for which the border is to
	 *  be drawn.
	 *  
	 *  @param height The height of the component for which the border is
	 *  to be drawn.
	 */
	@Override
	public void paintBorder(Component c, Graphics g, int x, int y, int width,
			int height) {
		Insets borderInsets = border.getBorderInsets(c); 
		Insets insets = getBorderInsets(c); 
		int temp = (insets.top-borderInsets.top)/2; 
		border.paintBorder(c, g, x, y+temp, width, height-temp); 
		Dimension size = comp.getPreferredSize(); 
		compBounds = new Rectangle(LEFT_OFFSET, 0, size.width, size.height); 
		SwingUtilities.paintComponent(g, comp, (Container)c, compBounds); 
	}

	/**
	 * This is a helper method that is invoked from the various methods
	 * associated with MouseListener to provide mouse events to {@link #comp}.
	 * 
	 * This method is a helper method that is essentially used to dispatch
	 * a mouse event from the border to the underlying {@link #comp}.
	 * 
	 * @param me The mouse event that is to be dispatched to the underlying
	 * component.
	 */
	private void dispatchEvent(MouseEvent me){ 
		if((compBounds !=null) && (compBounds.contains(me.getX(), me.getY()))) { 
			Point pt = me.getPoint(); 
			pt.translate(-LEFT_OFFSET, 0); 
			comp.setBounds(compBounds); 
			comp.dispatchEvent(new MouseEvent(comp, me.getID(),
					me.getWhen(), me.getModifiers(), pt.x, pt.y, 
					me.getClickCount(), me.isPopupTrigger(), me.getButton()));
			
			if(!comp.isValid()) {
				container.repaint();
			}
		} 
	} 

	@Override
	public void mouseClicked(MouseEvent e) {
		dispatchEvent(e);
	}

	 @Override
	 public void mousePressed(MouseEvent e) {
		 dispatchEvent(e);
	 }

	 @Override
	 public void mouseReleased(MouseEvent e) {
		 dispatchEvent(e);
	 }

	 @Override
	 public void mouseEntered(MouseEvent e) {
		 dispatchEvent(e);
	 }

	 @Override
	 public void mouseExited(MouseEvent e) {
		 dispatchEvent(e);
	 }

	 @Override
	 public Insets getBorderInsets(Component c) {
		 Dimension size = comp.getPreferredSize(); 
		 Insets insets = border.getBorderInsets(c); 
		 insets.top = Math.max(insets.top, size.height); 
		 return insets; 
	 }

	 @Override
	 public boolean isBorderOpaque() {
		 return true;
	 }
}
