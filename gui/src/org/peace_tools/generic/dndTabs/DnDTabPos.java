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

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.util.ArrayList;

import javax.swing.Icon;
import javax.swing.JOptionPane;
import javax.swing.JSplitPane;

/**
 * Class to remember position of a DnDTabPane.
 * 
 * This class provides a reusable mechanism to remember the relative position of
 * a DnDTabPane. This information is handy in the following cases:
 * 
 * <ul>
 * 
 * <li>When hiding and showing windows, it is handy to restore tabs in the same
 * location (as much as possible) where the user last placed it.</li>
 * 
 * <li>The same concept can be applied to saving and resotring positions of
 * windows in a persistent session file. Note that persisting to a session file
 * is delegated to the actual application classes.</li>
 * 
 * </ul>
 */
public class DnDTabPos {
	/**
	 * The default constructor.
	 * 
	 * The constructor has no special tasks to perform other than to initialize
	 * instance variables to their default initial value.
	 */
	public DnDTabPos() {
		prevName = "";
		prevIcon = null;
		visible  = true;
		path     = new ArrayList<DnDTabbedPane.Location>();
	}

	/**
	 * Determine the visibility information in this tab.
	 * 
	 * This method can be used to determine if the window (if any) logically
	 * associated with this object is visible or not.
	 * 
	 * @return This method returns true if the window is visible. Otherwise it
	 *         returns false.
	 */
	public boolean isVisible() {
		return visible;
	}

	/**
	 * Set the visibility information in this tab.
	 * 
	 * This method can be used to set if the window (if any) logically
	 * associated with this object is visible or not. Note that the visiblity
	 * flag is really independent of the window's visibility. It is used by the
	 * application to operate on its own semantics which may be independent of
	 * wxWiget's semantics.
	 * 
	 * @param visible
	 *            If this parameter is true then the window is assumed to be
	 *            visible. Otherwise the window is assumed to be hidden.
	 */
	public void setVisible(boolean visible) {
		this.visible = visible;
	}

	/**
	 * Set initial information for restoring a window.
	 * 
	 * This method can be used to initialize (or seed) this tab position helper
	 * with some information to restore a window.
	 * 
	 * @param icon
	 *            The icon to be used when restoring a window using this object.
	 * 
	 * @param name
	 *            The name of the tab to be used when restoring a window using
	 *            this object.
	 * 
	 * @param size
	 *            The preferred size of the window when restoring.
	 */
	public void setInfo(Icon icon, String name, Dimension size) {
		this.prevIcon = icon;
		this.prevName = name;
		this.prevSize = size;
	}

	/**
	 * Method to save the relative position of a given child window.
	 * 
	 * This method must be used to have this DndTabPos object to remember the
	 * relative position of a given child window. This method assumes that the
	 * child window is currently added to some DndTabPanel and is not orphaned
	 * (that is its parent pointer is still valid). Any previous position
	 * information stored in this class is lost.
	 * 
	 * @param desktop
	 * 			  The central desktop that contains all the panels including
	 *            the permanent panel. The same panel should be used when
	 *            restoring the window's position via the restorePosition() method. 
	 * @param tab
	 *            The tab (or child window) whose relative position must
	 *            determined and remembered by this class.
	 * 
	 * @return This method returns true if the relative position of this tab was
	 *         successfully saved. If the relative position could not be saved,
	 *         then this method returns false.
	 */
	public boolean savePosition(Container desktop, Component tab) {
	    // Clear the path from the permanent tab to the given tab.
	    path.clear();
	    // Ensure that the top-level desktop is not null.
	    if (desktop == null) {
	    	// Cannot save position without a permanent desktop
	    	return false;
	    }
	    // Get the first child in the desktop. Either it is a
	    // DnDTabbedPane (if windows are not split) or is a JSplitPane
	    // if windows are split.
	    assert(desktop.getComponentCount() == 1);
	    Component topChild = desktop.getComponent(0);
	    // Now have the recursive helper method do rest of the work.
	    boolean retVal = savePosition(topChild, tab);
	    if (!retVal) {
	    	// Huh! the window was not at all found. That is bad.
	    	JOptionPane.showMessageDialog(desktop, 
	    			"Panel whose position is to saved was not " +
	    			"found in the hierarchy!", "Position save error",
	    			JOptionPane.WARNING_MESSAGE);
	    }
	    return retVal;
	}

	/**
	 * Recursive helper method to save relative position.
	 * 
	 * This method is invoked from the overloaded savePosition() method with a
	 * suitable starting point. This method recursively calls itself until the
	 * DndTabPanel that contains the window whose location is to be remembered
	 * is reached. Once it finds the window, it saves the path traced as the
	 * recursion unwinds. Consequently, the path is saved in the reverse order
	 * of how the recursion proceeded.
	 * 
	 * @param currWin
	 *            The current window (may it be a DndTabPanel or a
	 *            wxSplitWindow) from where the search for tab must proceed.
	 * 
	 * @param tab
	 *            The actual window that we are searching for.
	 * 
	 * @return This method returns true if the position of the tab was found and
	 *         the path was saved. Otherwise this method returns false.
	 */
	boolean savePosition(Component currWin, Component tab) {
		// Check if currWindow is a DndTabPanel or if it is a
		// wxSplitterPane. Based on these two types the actions are a bit
		// different.
		assert ((currWin != null) && (tab != null));

		if (currWin instanceof DnDTabbedPane) {
			DnDTabbedPane dndPanel = (DnDTabbedPane) currWin;
			// Check if this dndPanel has the window we are looking
			// for. If so, it is a success.
			if (dndPanel.indexOfComponent(tab) != -1) {
				// Found it! First save information about it.
				saveInfo(dndPanel, tab);
				// Return true to indicate things went fine.
				return true;
			}
			// Not found and no further recursions.
			return false;
		}
	    // This must be a splitter pane. Search both halves..
	    JSplitPane splitWin = (JSplitPane)currWin;
	    // First search the top/left of the split window for the tab.
	    DnDTabbedPane.Location dir = DnDTabbedPane.Location.CENTER;
	    if (savePosition(splitWin.getLeftComponent(), tab)) {
	        // Found it in the top/left. Save the information.
	        dir = (splitWin.getOrientation() == JSplitPane.HORIZONTAL_SPLIT) ?
	             DnDTabbedPane.Location.LEFT : DnDTabbedPane.Location.TOP;
	    } else if (savePosition(splitWin.getRightComponent(), tab)) {
	        // Found it in the bottom/right. Save the information.
	        dir = (splitWin.getOrientation() == JSplitPane.HORIZONTAL_SPLIT) ?
	            DnDTabbedPane.Location.RIGHT : DnDTabbedPane.Location.BOTTOM;
	    } else {
	    	// Window was not found in this hierarchy.
	        return false;
	    }
	    
	    // Save the direction traversed
	    path.add(dir);
	    // Return true to indicate we saved the path
	    return true;
	}

	
	/**
	 * Helper method to save information about the given tab.
	 * 
	 * This is a helper method that is used to save information about a given
	 * child tab. The information saved includes the tab title text and the icon
	 * (if any) associated with the tab. This method updates the prevName and
	 * prevIcon instance variables.
	 * 
	 * @param owner
	 *            The owner tab panel from where the information about the child
	 *            tab is to be saved.
	 * @param tab
	 *            The child tab whose information is to be saved.
	 */
	private void saveInfo(DnDTabbedPane owner, Component tab) {
		// Find the index of the tab in the dnd panel in order to obtain
		// the necessary information.
		int idx = owner.indexOfComponent(tab);
		if (idx != -1) {
			// Found the index. Save the info.
			prevName = owner.getTitleAt(idx);
			// Save the icon for the tab (if any)
			prevIcon = owner.getIconAt(idx);
			// Save the size information.
			prevSize = tab.getSize();
		}
	}

	/**
	 * Primary interface to restore position of a tab.
	 * 
	 * This method is the complementary method to savePosition()  method
	 * in this class. This method performs the task of repositioning a 
	 * given tab within a DnDTabbedPane at its former position. This method
	 * does best effort to restore position of the tab. If the window
	 * layout has significantly changed, then this method places the window
	 * at at position at this closest to its pervious position without
	 * significantly impacting the current window layout.
	 * 
	 * @param desktop
	 * 			  The central desktop that contains all the panels including
	 *            the permanent panel. The same panel should be used when
	 *            saving the window's position via the savePosition() method.
	 *            
	 * @param tab The component whose position is to be restored.
	 */
	public void restorePosition(Container desktop, Component tab) {
	    DnDTabbedPane.Location dir[] = new DnDTabbedPane.Location[1];
	    dir[0] = DnDTabbedPane.Location.CENTER;
	    
	    // Get the first child in the desktop. Either it is a
	    // DnDTabbedPane (if windows are not split) or is a JSplitPane
	    // if windows are split.
	    assert(desktop.getComponentCount() == 1);	
	    Component currWin = getPosition(desktop.getComponent(0), dir);
	    // When control drops here, dir has the direction in which tab
		// must be placed with respect to currWindow. There are two
		// cases: (i) currWindow is a DndTabPanel and everything is good;
		// (ii) currWindow is a wxSplitterWindow in which case we need to
		// suitably re-split it to include the child at the correct
		// location.
	    if (currWin == null) {
	    	// This is no good. Log an error
	    	return;
	    }
	    // This if handles the case where currWin is a DnDTabbedPane
	    if (currWin instanceof DnDTabbedPane)  {
	    	DnDTabbedPane owner = (DnDTabbedPane) currWin;
	        // OK! the original location was found. Simple case.
	    	boolean verticalSplit = dir[0].equals(DnDTabbedPane.Location.BOTTOM) ||
	    		dir[0].equals(DnDTabbedPane.Location.TOP);
	        owner = owner.createSplitPane(prevName, prevIcon, tab, dir[0], 0.0, 
	        		verticalSplit ? prevSize.height : prevSize.width);
	        // Ensure the tab is visible and selected.
	        tab.setVisible(true);
	        owner.setSelectedComponent(tab);
	        return;
	    }
	    // The original location was not found. So have to place the
	    // window in this logical spot by suitably inserting a splitter
	    // window instead of this window.
	    JSplitPane splitWin = (JSplitPane) currWin;
	    // Create a splitWindow to hold a new DndTabPanel and currWin.
	    Container parent = splitWin.getParent();
	    // First create a DndTabPanel to hold this tab.
	    DnDTabbedPane sibling = new DnDTabbedPane();
	    sibling.addTab(prevName, prevIcon, tab);
	    final int index = sibling.indexOfComponent(tab);
	    sibling.setTabComponentAt(index, new DnDTabButton(sibling, index));

	    // Ensure direction is not CENTER in this case.
	    dir[0] = (dir[0].equals(DnDTabbedPane.Location.CENTER)) ? 
	    		DnDTabbedPane.Location.RIGHT : dir[0];
	    // Create a new split pan and configure the new split window
	    // depending on direction and add the new split window suitably
	    // using the helper method defined in the tabbed pane class.
	    DnDTabbedPane.setSplitWindow(parent, dir[0], splitWin, sibling,
	                               prevSize.width, prevSize.height, 
	                               splitWin.getSize());
	}
	
    /**
	 * Traverse saved position information as much as possible.
	 * 
	 * This is a helper method that was introduced to keep the code clutter in
	 * the restorePosition() to a minimum. This method does best effort to
	 * traverse the saved position information, while tracking the sequence of
	 * windows encountered.
	 * 
	 * @param permTab The top-level permanent tab within which the position
	 * is to be determined.
	 * 
	 * @param direction The last direction that was effectively traversed
	 * using the saved position information is set as the first element in
	 * this array. The caller must allocate this array to contain at least 1 
	 * element.
	 * 
	 * @return The last window (split window or DndTabPanel as the case maybe)
	 * that was encountered during the traversal.
	 */
	protected Component getPosition(Component permTab, 
			DnDTabbedPane.Location direction[]) {
	    // If the path is empty, then simply add the window to the
	    // center panel, which is the default location.
	    if (path.size() == 0) {
	        // Just return the center panel after setting direction
	        direction[0] = DnDTabbedPane.Location.CENTER;
	        return permTab;
	    }
	    
	    // Now keep traversing the saved path in reverse until we cannot
	    // proceed any further.
	    Component currWin = permTab;
	    int idx;
	    DnDTabbedPane.Location curr = DnDTabbedPane.Location.CENTER;
	    for(idx = path.size() - 1; (idx >= 0); idx--) {
	    	// Update current movement being considered.
	    	curr = path.get(idx);
	    	if (!(currWin instanceof JSplitPane)) {
	    		// A non-splitter window encountered. Cannot proceed
	            // further.
	            break;
	    	}
	    	JSplitPane splitWin = (JSplitPane) currWin;
	        // Update the currWindow depending on the direction of
	        // movement and move only if possible.
	        if (((curr.equals(DnDTabbedPane.Location.TOP)) || 
	        	 (curr.equals(DnDTabbedPane.Location.BOTTOM))) &&
	            (splitWin.getOrientation() != JSplitPane.VERTICAL_SPLIT)) {
	            // We wan't to traverse up/down but can't.  So stop here.
	            break;
	        }
	        if (((curr.equals(DnDTabbedPane.Location.LEFT)) || 
		        	 (curr.equals(DnDTabbedPane.Location.RIGHT))) &&
		            (splitWin.getOrientation() != JSplitPane.HORIZONTAL_SPLIT)) {
		            // We wan't to traverse left/right but can't.  So stop here.
		            break;
		    }
	        // OK, keep traversing in the correct direction.
	        if ((curr.equals(DnDTabbedPane.Location.TOP)) || 
	        	(curr.equals(DnDTabbedPane.Location.LEFT))) {
	            currWin = splitWin.getTopComponent();
	        } else {
	            currWin = splitWin.getBottomComponent();
	        }
	    }
	    // Update direction if we have reached our final destination
	    // successfully.
	    if ((currWin instanceof DnDTabbedPane) && (idx < 0)) {
	        // We have successfully reached our final destination
	    	curr = DnDTabbedPane.Location.CENTER;
	    }
	    // Update last direction that we tried to traverse
	    direction[0] = curr;
	    // Return the current/final panel
	    return currWin;
	}

	/**
	 * The name for the tab when its position is restored.
	 * 
	 * This instance variable contains the name to be used for the tab when its
	 * position is restored. This information is updated in teh savePosition()
	 * method.
	 */
	String prevName;

	/**
	 * The icon (if any) to be used for the tab.
	 * 
	 * This instance variable contains a pointer to the icon to be used when
	 * restoring the tab. This pointer is updated in the savePosition() method.
	 */
	Icon prevIcon;

	/**
	 * The preferred width / height of the tab.
	 * 
	 * This instance variable is used to maintain the previous preferred width
	 * or height of the window. This information is updated in the
	 * savePosition() method. It is used only in the case that a new splitter
	 * window needs to be created when the tab is restored.
	 */
	Dimension prevSize;

	/**
	 * An additional boolean to track if window is visible.
	 * 
	 * This flag is used to track if the window associated with this DndTabPos
	 * object is currently visible or not. Note that this flag is really
	 * independent of the window's visibility. It is used by the application to
	 * operate on its own semantics which may be independent of wxWiget's
	 * semantics.
	 */
	boolean visible;

	/**
	 * Vector to maintain path from center pane to the desired position of a
	 * tab.
	 * 
	 * This vector is used to hold the necessary direction information that is
	 * used to remember the relative position of the given window with respect
	 * to the center pane. Entries in this vector are added by the
	 * savePosition() method and used by the restorePosition() methods.
	 */
	private ArrayList<DnDTabbedPane.Location> path;
}
