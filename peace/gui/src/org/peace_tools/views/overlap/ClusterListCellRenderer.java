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
import java.awt.Component;
import java.awt.Graphics;

import javax.swing.BorderFactory;
import javax.swing.DefaultListCellRenderer;
import javax.swing.Icon;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.border.Border;

import org.peace_tools.data.ClusterNode;

/**
 * A custom cell renderer to efficiently display large cluster list.
 *
 * <p>This cell renderer is meant to be used to render entries in a
 * JList. Specifically it is used to render cluster entries in the
 * cluster list owned by the OverlapView class. This cell renderer
 * utilizes the base class for most of the tasks and then adds a 
 * icon to display the color currently set for this cluster. The
 * color used for a cluster is obtained via the ClusterColorMapper
 * interface method(s).  Therefore the constructor requires an 
 * ClusterColorMapper object to be passed in as a parameter.</p>   
 *
 * <p><b>Note</b>This class is meant to be instantiated and used
 * only by the OverlapView class. Consequently it has been made 
 * package private.</p>
 */
class ClusterListCellRenderer extends DefaultListCellRenderer {
	/**
	 * The default constructor. 
	 * 
	 * @param overlapPanel The overlap panel from where the color information
	 * is to be obtained for a given cluster entry. 
	 */
	public ClusterListCellRenderer(ClusterColorMapper overlapPanel) {
		this.colorMapper = overlapPanel;
		this.colorIcon   = new CustomIcon();
		this.border      = BorderFactory.createEmptyBorder(1, 5, 1, 5);
	}
	
	/**
	 * Overrides base class method to provide customized color icon.
	 * 
	 * This method overrides the default implementation in the parent class
	 * to add a customized icon to the label returned by this method.
	 * 
	 * @param list The list in which this renderer is currently being used.
	 * @param value The value for which a cell renderer is required. This value
	 * must be a ClusterNode in order for a custom icon to be generated.
	 * @param index The logical index of the value in the list.
	 * @param isSelected Flag to indicate if this value is currently selected
	 * in the list.
	 * @param cellHasFocus Indicates if this cell currently has focus.
	 * @return A JLabel with the necessary information and a icon to indicate
	 * color.
	 */
	@Override
	public Component getListCellRendererComponent(JList list, Object value,
			int index, boolean isSelected, boolean cellHasFocus) {
		Component renderer = super.getListCellRendererComponent(list, value, index, isSelected,
				cellHasFocus);
		if ((value instanceof ClusterNode) && (renderer instanceof JLabel)) {
			// Now we set up the icon to be used for this label to represent the
			// color being used to color the cluster.
			ClusterNode node = (ClusterNode) value;
			// Get the color code for the cluster node.
			Color color = colorMapper.getColor(node.getClusterId());
			if (color == null) {
				color = Color.lightGray;
			}
			colorIcon.setColor(color);
			// Set the icon in the label.
			JLabel label = (JLabel) renderer;
			label.setIcon(colorIcon);
			// Setup some empty space to make things look nicer.
			label.setBorder(border);
		}
		// Return the adapted or default renderer back to the caller.
		return renderer;
	}

	/**
	 * A simple, custom icon class.
	 * 
	 * This icon class provides a custom Icon implementation to render
	 * a simple square box, filled with a specific color. The component
	 * is sufficiently simple and has a very small memory foot print.
	 * It paints the icon each time the 
	 * {@link CustomIcon#paintIcon(Component, Graphics, int, int)}
	 * method is invoked.
	 */
	private class CustomIcon implements Icon {
		/**
		 * Default constructor.
		 * 
		 * The constructor sets the default color to light gray.
		 */
		public CustomIcon() {
			color = Color.lightGray;
		}
		
		/**
		 * Set the color to be used to render the icon.
		 * 
		 * @param color The color to be used to fill the square area
		 * of the icon.
		 */
		public void setColor(Color color) {
			this.color = color;
		}
		
		@Override
		public int getIconHeight() {
			return ICON_SIZE;
		}

		@Override
		public int getIconWidth() {
			return ICON_SIZE;
		}

		@Override
		public void paintIcon(Component c, Graphics g, int x, int y) {
			g.setColor(color);
			g.fill3DRect(x + 1, y + 1, ICON_SIZE - 2, ICON_SIZE - 2, true);
			// Draw a black border to make the icon look nice.
			g.setColor(Color.black);
			g.drawRect(x, y, ICON_SIZE - 1, ICON_SIZE - 1);
		}
		
		/**
		 * The color with which the square area of the icon is to be
		 * painted.
		 */
		private Color color;
	}
	
	/**
	 * The icon that is used in the label created by the 
	 * getListCellRendererComponent method. The same icon is reused 
	 * over-and-over to try and reduce the  memory footprint.
	 */
	private final CustomIcon colorIcon;

	/**
	 * The ClusterColorMapper that provides the coloring information
	 * for rendering the colored icons with each entry.
	 * This instance variable is initialized by the constructor and
	 * is never changed during the life time of this object.
	 */
	private final ClusterColorMapper colorMapper;
	
	/**
	 * A simple (typically empty) border around each entry to make
	 * the list appear more pleasing to the eye.
	 */
	private final Border border;
	
	/**
	 * The size of the square to be used to display the color being
	 * used for rendering the cluster.
	 */
	private static final int ICON_SIZE = 12;
	
	/**
	 * Generated serialization UID (merely to keep the compiler happy)
	 */
	private static final long serialVersionUID = -1180299364821332033L;
}
