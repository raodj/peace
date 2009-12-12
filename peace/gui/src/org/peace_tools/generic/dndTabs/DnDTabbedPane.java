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

package org.peace_tools.generic.dndTabs;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.util.ArrayList;

import javax.swing.Icon;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.plaf.SplitPaneUI;
import javax.swing.plaf.basic.BasicSplitPaneUI;

public class DnDTabbedPane extends JTabbedPane {
	/**
	 * A generated serial version ID.
	 */
	private static final long serialVersionUID = 1414403952944170345L;

	/**
	 * This constant determines the default size that is used for panels
	 * that are docked to the edges (top, left, bottom, right) of a given
	 * DnDTabPanel. Setting this value to 0.5 will cause the panels to
	 * be split 50/50.
	 */
	private static final double EDGE_PANEL_SIZE = 0.25;
	
	/**
	 * A boolean variable that indicates if this tab is permanent and should not
	 * be deleted even if it is empty. This flag is set via call to
	 * setPermanent() method. This value is used by the TabDnDHanlder to decide
	 * if this tab should be removed when empty or not. It access's this property
	 * via the isPermanent() method.
	 */
	private boolean permanent;

	/**
	 * This instance variable is used to track if this DndTabbedPane should
	 * paint a outline on itself to indicate where a tab being dragged
	 * will be docked with respect to this window. This outline provides a
	 * visual cue to the user as to where a tab being dragged will end up.
	 * The visual cue is setup by the TabDnDHandler's methods.
	 */
	private Location dockCue = null;
	
	/**
	 * The list of listeners that have been added for this DnDTabbedPane.
	 * The listeners are added via the addListener() method or removed
	 * via the removeListener() method.
	 * 
	 * @note The list of listeners is static so that all DnDTabbedPane
	 * objects, including those created dynamically when the user 
	 * drags and drops tabs, will use this list of listeners. This is
	 * a bit non-standard implementation but appears to be the best
	 * given the dynamic nature of creation of DnDTabbedPanes
	 */
	private static final ArrayList<DnDTabListener> listeners
		= new ArrayList<DnDTabListener>();
	
	/**
	 * Preferred docking locations for tabs.
	 * 
	 * This enumeration value specifies valid docking locations for a given tab.
	 * These values are interpreted with respect to the tab (\c this) on which
	 * various methods are called. The \c CENTER enumeration specifies that a
	 * tab must be added to the DndTabbedPane itself rather than adjacent to it.
	 */
	public enum Location {
		TOP, BOTTOM, LEFT, RIGHT, CENTER
	};

	/**
	 * The default constructor. This constructor creates an empty TabbedPane
	 * with a default tab placement of JTabbedPane.TOP. It also creates a new
	 * DnD handler for this tabbed pane.
	 * 
	 * @param parent The parent tab from where some of the properties
	 * (such as listeners) are to be inherited by this new tab.
	 */
	public DnDTabbedPane() {
		super(TOP, SCROLL_TAB_LAYOUT);
		new DnDTabHandler(this);
		permanent = false;
		setBorder(null);
	}

	/**
	 * This method may be used to determine if this tabbed pane should be
	 * permanently present or not. Permanent tabbed panes are not removed by the
	 * DnDHandler even if it is empty.
	 * 
	 * @return Returns true if this tabbed pane is permanent.
	 */
	public boolean isPermanent() {
		return permanent;
	}

	/**
	 * This method should be used to set if a tabbed pane should be permanent or
	 * not. If the parameter is true then the tabbed pane is set to be permanent
	 * for all subsequent DnD operations.
	 * 
	 * @param perm
	 *            Flag to indicate if this tabbed pane is permanent or not. True
	 *            value indicates the tab is permanent and should not be removed
	 *            even if it is empty.
	 */
	public void setPermanent(boolean perm) {
		permanent = perm;
	}

	/**
	 * This method may be used to manually add components either directly to
	 * this tabbed pane or around this tabbed pane.
	 * 
	 * @param title
	 *            The title to be used for new tab when created.
	 * @param icon
	 *            The icon to be used if any.
	 * @param c
	 *            The actual component to be added.
	 * @param location
	 *            The preferred location around this pane.
	 * 
	 * @return The newly created DndTabbedPane or this.
	 */
	public DnDTabbedPane createSplitPane(String title, Icon icon, Component c,
			Location location) {
		double gravity = 1 - EDGE_PANEL_SIZE;
		if ((location == Location.TOP) || (location == Location.LEFT)) {
			gravity = EDGE_PANEL_SIZE;
		}
		return createSplitPane(title, icon, c, location, gravity);
	}

	/**
	 * This method may be used to manually add components either directly to
	 * this tabbed pane or around this tabbed pane.
	 * 
	 * @param title
	 *            The title to be used for new tab when created.
	 * @param icon
	 *            The icon to be used if any.
	 * @param c
	 *            The actual component to be added.
	 * @param location
	 *            The preferred location around this pane.
	 *            
	 * @param gravity The size into which the two panes must be split.
	 * 
	 * @return The newly created DndTabbedPane or this.
	 */
	public DnDTabbedPane createSplitPane(String title, Icon icon, Component c,
			Location location, double gravity) {
		return createSplitPane(title, icon, c, location, gravity, 0);
	}

	/**
	 * This method may be used to manually add components either directly to
	 * this tabbed pane or around this tabbed pane.
	 * 
	 * @param title
	 *            The title to be used for new tab when created.
	 * @param icon
	 *            The icon to be used if any.
	 * @param c
	 *            The actual component to be added.
	 * @param location
	 *            The preferred location around this pane.
	 *            
	 * @param gravity The ratio of space to be allocated to each window.
	 * This value is set for JSplitPane.setResizeWeight 
	 * 
	 * @param size The size to be set for the component to be added. If
	 * this size is less than the minimum size of the component c, then
	 * the minimum size of component c will be used. If the size is 
	 * negative then a size is not set and the defaults are used.
	 * 
	 * @return The newly created DndTabbedPane or this.
	 */
	public DnDTabbedPane createSplitPane(String title, Icon icon, Component c,
			Location location, double gravity, int size) {
		// If the location is center there is not a whole lot
		// to be done. We can proceed merrily.
		if (location == Location.CENTER) {
			addTab(title, icon, c);
			// Setup the title rendering component
			final int index = indexOfComponent(c);
			setTabComponentAt(index, new DnDTabButton(this, index));
			return this;
		}
		// Obtain a reference to the parent of this current container
		// for future use a few lines below to ensure we don't loose
		// the reference
		Container parent = getParent();

		// Create a DnDpane to hold the new component.
		DnDTabbedPane newPane = new DnDTabbedPane();
		newPane.addTab(title, icon, c);
		final int index = newPane.indexOfComponent(c);
		newPane.setTabComponentAt(index, new DnDTabButton(newPane, index));

		// Determine the width/height that the new split pane should 
		// allocate to the new component based on our current area.
		Dimension currArea = parent.getSize();
		if ((currArea.width == 0) || (currArea.height == 0)) {
			// If current size is zero or not set then use the preferred size
			// currArea = parent.getPreferredSize();
		}
		// Compute preferred width and height while giving preference to the
		// specified size if the preferred dimensions turn out to be zero.
		int prefWidth  = (int) (gravity * currArea.width);
		prefWidth      = (prefWidth == 0) ? size : prefWidth;
		int prefHeight = (int) (gravity * currArea.height);
		prefHeight     = (prefHeight == 0) ? size : prefHeight;
		
		JSplitPane jsp = setSplitWindow(parent, location, this, newPane, 
				prefWidth, prefHeight); 
 		// Setup initial dimension and gravity
		jsp.setResizeWeight(gravity);

		// Return the newly created Dnd pane.
		return newPane;
	}
	
	/**
	 * This is a utility method that is used by both DnDTabbedPane and DnDTabPos
	 * to create a JSplitPane containing two components organized in a specific
	 * manner to occupy a given preferred width/height depending on the type of
	 * orientation. This method specifically replaces window1 in the parent with
	 * a new split pane containing window1 and window2. This method handles both
	 * cases where the parent can be JSplitPane or just a regular container. If
	 * the parent is JSplitPane this method suitably replaces window1 with the
	 * newly created split pane. In the second case, the newly created split
	 * pane is just directly added to the parent.
	 * 
	 * @param parent
	 *            The parent component that contains window1
	 * @param location The location to decide where window2 should be placed. This
	 * location cannot be CENTER.
	 * @param window1 The component to be replaced with a new split window containing
	 * window1 and window2.
	 * @param window2 The second component to be placed in the newly created split
	 * pane.
	 * @param prefWidth The preferred width of window1 if the location requires a 
	 * vertical split.
	 * @param prefHeight The preferred height of window1 if the location requires a 
	 * horizontal split.
	 * @return The newly created split window for further customization.
	 */
	public static JSplitPane setSplitWindow(Container parent, Location location,
			Component window1, Component window2, int prefWidth, int prefHeight) {
		// First create a JSplitPane with suitable orientation.
		int orientation = JSplitPane.HORIZONTAL_SPLIT;
		if (location.equals(Location.TOP) || location.equals(Location.BOTTOM)) {
			orientation = JSplitPane.VERTICAL_SPLIT;
		}
		JSplitPane jsp = new JSplitPane(orientation);
		jsp.setMinimumSize(new Dimension(60, 60));
		jsp.setDividerSize(2);
		// The next 5 lines of code attempt to remove an annoying
		// border around the divider in the split pane. Currently,
		// the code removes only one of the 2 border lines. Unsure
		// how to remove the other line so that the divider looks
		// nice and clean.
		SplitPaneUI ui = jsp.getUI();
		if (ui instanceof BasicSplitPaneUI) {
			BasicSplitPaneUI basicUI = (BasicSplitPaneUI) ui;
			basicUI.getDivider().setBorder(null);
		}
		// Add the two windows to the split pane.
		if (orientation == JSplitPane.HORIZONTAL_SPLIT) {
			jsp.setLeftComponent(location.equals(Location.LEFT) ? 
					window2	: window1);
			jsp.setRightComponent(location.equals(Location.LEFT) ? 
					window1	: window2);
			jsp.setDividerLocation(prefWidth);
		} else {
			jsp.setTopComponent(location.equals(Location.TOP) ? 
					window2	: window1);
			jsp.setBottomComponent(location.equals(Location.TOP) ? 
					window1	: window2);
			jsp.setDividerLocation(prefHeight);
		}
		jsp.setBorder(null);
		
	    // Here the logic for adding/moving windows is different depending
	    // on whether the parent is a split window or not.
	    if (parent instanceof JSplitPane) {
	    	JSplitPane splitParent = (JSplitPane) parent;
	    	// Replace window1 in the parent with new split pane.
	    	if (splitParent.getTopComponent() == null) {
	    		splitParent.setTopComponent(jsp);
	    	} else {
	    		assert ( splitParent.getBottomComponent() == window1 );
	    		splitParent.setBottomComponent(jsp);
	    	}
	    } else {
	        // The parent is not a splitter window. So the only
	        // thing left to do is to add the split window to the parent.
	    	// But first remove window1 from the parent.
	    	parent.remove(window1);
	        parent.add(jsp);
	    }
	    // Ensure that the parent relay's out components
	    parent.validate();
		// Finally return the newly created split pane.
		return jsp;
	}

	/**
	 * Delete (and remove) a given tab component.
	 * @param componentRemoved
	 */
	public void deleteTab(Component componentRemoved) {
		removeTab(componentRemoved);
		// Remote to listeners that the component has been removed
		notifyListeners(componentRemoved);
	}
	
	/**
	 * Remove a tab (or component) from this DnDTabbedPane. When the component
	 * is removed, and the pane is empty, then the pane automatically removes
	 * itself from its parent.
	 * 
	 * @param componentRemoved
	 *            The component to be removed from this pane.
	 */
	public void removeTab(Component componentRemoved) {
		// Remove the tab from ourselves.
		remove(componentRemoved);
		// Do we have any tabs left or is it permanent?
		if ((getTabCount() > 0) || (isPermanent())) {
			// Yes. we do. Good. nothing further to do.
			return;
		}
		// When control drops here, that means we have zero tabs
		// left and it is time to remove ourselves from our parent.
		Container parent = getParent();
		parent.remove(this);
		// If our parent is a JSplit pane we need to remove it too
		// because it has only 1 component and should no longer exist.
		if (parent instanceof JSplitPane) {
			// Our parent was a JSplitPane which needs to be
			// nuked now because it has only 1 component.
			JSplitPane jsp = (JSplitPane) parent;
			Component other = (jsp.getLeftComponent() != null) ? jsp
					.getLeftComponent() : jsp.getRightComponent();
			parent = jsp.getParent();
			parent.remove(jsp);
			parent.add(other);
		}
		parent.validate();
	}

	/**
	 * Helper method to report listeners that a component has been removed.
	 * 
	 * This method is invoked from the removeTab() method after a 
	 * component has been removed from the tab. This method notifies
	 * all registered listeners that the component has been removed.
	 * 
	 * @param The component that has been removed.
	 */
	protected static void notifyListeners(Component comp) {
		for(DnDTabListener listener: listeners) {
			listener.tabDeleted(comp);
		}
	}

	/**
	 * Add a listener to be notified when a tab in this panel is deleted.
	 * 
	 * @note The listener is not recursively propagated to all associated 
	 * tabs. However, the listener is inherited by any new DnDTabbedPane
	 * that may be created from this pane due to user dragging tabs.
	 * 
	 * @param listener The listener to be added.
	 */
	public static void addListener(DnDTabListener listener) {
		listeners.add(listener);
	}

	/**
	 * Remove a listener to be notified from this tab.
	 * 
	 * @note The listener is not recursively removed from all associated
	 * tabs. 
	 * 
	 * @param listener The listener to be removed from the set of 
	 * listeners to be notified.
	 */
	public static void removeListener(DnDTabListener listener) {
		listeners.remove(listener);
	}

	/**
	 * This method can be used to determine the current docking cue
	 * associated with this DndTabPane. The dock cue is used to draw
	 * a rectangle on this tab indicating where a tab will be docked
	 * to provide visual cues to the user. This method is currently 
	 * used only by the TabDnDHandler.
	 * 
	 * @return The current dock cue set for this tab. If a dock cue is
	 * not set, then this method returns null. 
	 */
	public Location getDockCue() { 
		return dockCue;
	}
	
	/**
	 * This method must be used to set/clear the current docking cue
	 * associated with this DndTabPane. The dock cue is used to draw
	 * a rectangle on this tab indicating where a tab will be docked
	 * to provide visual cues to the user. This method is currently 
	 * used only by the TabDnDHandler.
	 * 
	 * @return The dock cue to be set for this tab. If a dock cue is
	 * set to null then a cue is not drawn by this tab. 
	 */
	public void setDockCue(Location cue) {
		dockCue = cue;
	}
	
	/**
	 * This method overrides the default method in the parent class.
	 * This method first calls the parent's paintChildren() method.
	 * After that it checks and draws the cue box if the docCue
	 * instance variable is not null.
	 * 
	 * @param g The graphics object to be used for drawing the cue
	 * box.
	 */
	@Override
	public void paintChildren(Graphics g) {
		super.paintChildren(g.create());
		if (dockCue == null) {
			// No cues have been set. So nothing further to do.
			return;
		}
		// Draw a cue at the appropriate location.
		Rectangle cueBox = new Rectangle(0, 0, getWidth(), getHeight());
		cueBox.x = cueBox.y = 0;
		switch (dockCue) {
		case BOTTOM: cueBox.y += (int) (cueBox.height * (1- EDGE_PANEL_SIZE));
		case TOP: cueBox.height *= EDGE_PANEL_SIZE;
		break;
		case RIGHT: cueBox.x += (cueBox.width * (1 - EDGE_PANEL_SIZE));
		case LEFT: cueBox.width *= EDGE_PANEL_SIZE;
		}
		// If the graphics object is Graphics2D then we setup a minor
		// effect of drawing the cue rectangle using dotted lines
		// instead of solid lines.
		if (g instanceof Graphics2D) {
			float Dashes[] = {1.0f, 1.0f, 1.0f, 1.0f};
			BasicStroke lineStyle = 
				new BasicStroke(1.0f, BasicStroke.CAP_ROUND, 
						BasicStroke.JOIN_ROUND, 0.0F, Dashes, 1.0f);
			Graphics2D g2d = (Graphics2D) g;
			g2d.setStroke(lineStyle);
		}
		// Draw four concentric rectangles as the cue box.
		g.setColor(new Color(0, 0, 0, 128));
		for(int rep = 0; (rep < 4); rep++) {
			g.drawRect(cueBox.x, cueBox.y,cueBox.width, cueBox.height);
			cueBox.x++;
			cueBox.y++;
			cueBox.width -= 2;
			cueBox.height -= 2;
		}
	}
}
