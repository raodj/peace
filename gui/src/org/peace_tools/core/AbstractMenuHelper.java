

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

package org.peace_tools.core;

import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.AbstractButton;
import javax.swing.JMenuItem;
import javax.swing.JTable;
import javax.swing.JTree;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.event.TreeSelectionListener;

import org.peace_tools.generic.Utilities;


/**
 * The base class for various menu helpers.
 * 
 * The base class provides the common infrastructure for all menu
 * helpers in the GUI. In addition, it also serves as a common 
 * interface to obtain and use various menu helpers in the core
 * package.
 * 
 * @see MainFrame#getMenuHelper(org.peace_tools.core.AbstractMenuHelper.HelperType)
 */
public abstract class AbstractMenuHelper implements TableModelListener {
	/**
	 * Enumeration to enable convenient identification/classification
	 * of menu helper objects.
	 * 
	 * The following enumerations provide a convenient mechanism to
	 * identify a menu helper directly from the base class. Each 
	 * menu helper in the system has a unique enumerated identifier
	 * associated with it. Typically, each menu in the main menu 
	 * has its own menu helper. However, the system may have additional
	 * menu helpers that may be associated with a sub-menu or possibly
	 * not associated with a menu at all (as they are used in other
	 * spots in the system).
	 * 
	 * @see AbstractMenuHelper#getHelperType()
	 */
	public enum HelperType {
		FILE_MENU, WORKSPACE_MENU, SERVER_MENU, JOB_MENU, VIEW_MENU, HELP_MENU
	}
	
	/**
	 * Enumeration to enable referring to a specific type of action.
	 * 
	 * This enumeration provides an exhaustive list of various actions
	 * that are supported by various helpers in the system. Note that
	 * not all helpers will support all the actions. For example, the
	 * FILE_MENU helper will only support the options associated with
	 * performing tasks classified as file menu. Refer to the
	 * documentation on each derived helper class for details on the
	 * actions supported by each helper.
	 * 
	 * <p><b>Note:</b>  Each menu helper only supports a subset of these actions.
	 * The objective is to achieve a reasonable degree of "separation
	 * of concerns" so that the code is streamlined and GUI performance
	 * is not negatively impacted.</p>
	 */
	public enum ActionType {
		SAVE_WORKSPACE, SWITCH_WORKSPACE, NEW_DATASET, 
		IMPORT_DATASET, EXIT_PEACE, OPEN_DEFAULT_VIEW, 
		OPEN_CLUSTER_SUMMARY, OPEN_AS_TEXT, // FILE_MENU
		
		ADD_SERVER, REMOVE_SERVER, SHOW_MY_JOBS,
		SHOW_ALL_JOBS, CONNECTION_TEST, // SERVER_MENU
		
		COMPUTE_MST, COMPUTE_CLUSTERS, START_JOB_MONITOR, STOP_JOB_MONITOR, 
		SHOW_JOB_DETAILS, ABORT_JOB, REMOVE_JOB, 
		SHOW_JOBS_ON_SERVER, SHOW_MY_JOBS_ON_SERVER, // JOB_MENU

		NEW_MAIN_FRAME, DUPLICATE_EDITOR, 
		SHOW_WORKSPACE_HIERARCHY_BROWSER, SHOW_WORKSPACE_FILE_BROWSER,
		VIEW_JOBS, VIEW_SERVERS, SHOW_USER_LOGS, SHOW_PROG_LOGS, // VIEW_MENU
		
		SHOW_HELP_CONTENTS, SHOW_USER_TUTORIALS, SHOW_PAPERS,
		SHOW_EMAIL, SHOW_UPDATES, SHOW_PROGRAMMER_DOCS, SHOW_NOTES,
		SHOW_BUGS, SHOW_QUICK_START, SHOW_ABOUT_DIALOG // HELP_MENU
	}
	
	public AbstractMenuHelper(HelperType helperType, MainFrame mainFrame) {
		this.helperType = helperType;
		this.mainFrame  = mainFrame;
	}
	
	
	/**
	 * Determine the type of helper that this object represents.
	 * 
	 * This method can be used to determine the type of helper that
	 * this menu helper represents. This information is useful when
	 * dealing with menu helpers and creating tools because certain
	 * menu helpers only implement certain functionality.
	 * 
	 * @return This method returns the type of helper that was 
	 * specified when this menu helper was instantiated.
	 */
	public HelperType getHelperType() {
		return helperType;
	}
	
	/**
	 * Obtain a main menu item for the given action type.
	 * 
	 * This method must be used to create a suitable main menu entry
	 * for the given type of action. This action type must be supported
	 * by the menu helper in order to obtain a valid menu item. If
	 * the action is not supported, then this method returns null.
	 * 
	 * <p><b>Note:</b>  The returned menu item has all the necessary information
	 * already filled-in. Do not modify the action command and other
	 * properties as it will interfere with correct operation of the
	 * menu item.</p>
	 * 
	 * @param actionType The action type for which a main menu item
	 * is to be created.
	 * 
	 * @param mainMenu If this flag is true, then it indicates that
	 * the menu should be created to be used in the main menu. Otherwise
	 * it is assumed that the menu item will be used in a context 
	 * sensitive popup menu.
	 * 
	 * @return A JMenuItem corresponding to the given type of action.
	 * If the specified actionType is not supported by this menu 
	 * helper, then this method returns null.
	 */
	public abstract JMenuItem getMenuItem(ActionType actionType, boolean mainMenu);
	
	/**
	 * The action listener to handle actions created by this menu helper.
	 * 
	 * This method is a convenience method that can be used to obtain
	 * the actual action listener associated with this menu helper. The
	 * action listener is responsible for handling action call backs
	 * associated with the tools created by this menu helper.
	 * 
	 * @return The action listener associated with this menu helper.
	 * This action listener must be used only with the tools and
	 * menu items created by this menu helper.
	 */
	public abstract ActionListener getActionListener();
	
	/**
	 * The listener to be used in JTree to update tools & menu items.
	 * 
	 * This method is a convenience method that can be used to obtain
	 * a tree selection listener associated with this menu helper. The
	 * tree listener must be added to suitable JTree that display
	 * appropriate Workspace entries. The tree selection listener
	 * intercepts events generated by Java when the user clicks on
	 * specific entries in a JTree. This information is used to 
	 * suitable enable/disable the tools generated by this listener.
	 * 
	 * @param tree The JTree to which the tree selection listener is
	 * going to be added. This reference is maintained by the menu
	 * helper to handle tree selection events.
	 * 
	 * @return The tree listener associated with this menu helper.
	 * If this listener does not provide a tree selection listener
	 * then this method returns null.
	 */
	public abstract TreeSelectionListener getTreeSelectionListener(JTree tree);
	
	/**
	 * The listener to be used in JTable to update tools & menu items.
	 * 
	 * This method is a convenience method that can be used to obtain
	 * a list selection listener associated with this menu helper. The
	 * list listener must be added to suitable JTable that display
	 * appropriate Workspace entries. The table selection listener
	 * intercepts events generated by Java when the user clicks on
	 * specific entries in a JTable. This information is used to 
	 * suitable enable/disable the tools generated by this listener.
	 * 
	 * <p><b>Note:</b>  The list selection listeners handle only single selection
	 * model for the table.</p>
	 * 
	 * @param table The JTable to which the list selection listener is
	 * going to be added. This reference is maintained by the menu
	 * helper to handle row selection events.
	 * 
	 * @return The table listener associated with this menu helper.
	 * If this listener does not provide a table selection listener
	 * then this method returns null.
	 */
	public abstract ListSelectionListener getListSelectionListener(JTable table);

	/**
	 * Obtain a tool bar button for the given action type.
	 * 
	 * This method must be used to create a suitable tool bar  entry
	 * for the given type of action. This action type must be supported
	 * by the menu helper in order to obtain a valid menu item. If
	 * the action is not supported, then this method returns null.
	 * 
	 * <p><b>Note:</b>  The returned button has all the necessary information
	 * already filled-in. Do not modify the action command and other
	 * properties as it will interfere with correct operation of the
	 * too bar button.</p>
	 * 
	 * @param actionType The action type for which a tool bar button
	 * is to be created.
	 * 
	 * @param mainToolBar If this flag is true, then it indicates that
	 * the tool will be used in the main tool bar in PEACE. Otherwise 
	 * it is is assumed that the tool will be used in a view specific
	 * tool bar. 
	 * 
	 * @return An abstract button corresponding to the given type of 
	 * action. If the specified actionType is not supported by this 
	 * menu helper, then this method returns null.
	 */
	public abstract AbstractButton getTool(ActionType actionType, boolean mainToolBar);
	
	/**
	 * Helper method to update status of all items for a given action command.
	 * 
	 * This method iterates over the list of items that have been added
	 * to the contextItemList and changes status of those items that have
	 * the same action command as the one specified.
	 * 
	 * @param actionCommand The action command associated with a given set
	 * of items.
	 * 
	 * @param enabled If this flag is true, then the item is enabled. Otherwise
	 * the item is disabled.
	 */
	protected void setEnabled(String actionCommand, boolean enabled) {
		for(AbstractButton item: contextItemList) {
			if (actionCommand.equals(item.getActionCommand())) {
				item.setEnabled(enabled);
				// For HTML menu items, modify their status via
				// HTML color codes.
				if (item instanceof JMenuItem) {
					Utilities.setEnabled((JMenuItem) item, enabled);
				}
			}
		}
	}
	
	
	/**
	 * Default implementation for TableModelListener.
	 * 
	 * This method intercepts changes occurring to the table and 
	 * translates the event to a suitable ListSelectionListener 
	 * call if the derived class provides a suitable 
	 * ListSelectionLister.
	 */
	@Override
	public void tableChanged(TableModelEvent tme) {
		if ((table != null) && (this instanceof ListSelectionListener)) {
			ListSelectionListener lsl = (ListSelectionListener) this;
			ListSelectionEvent    lse = new ListSelectionEvent(tme.getSource(), 
					table.getSelectedRow(), table.getSelectedRow(), false);
			lsl.valueChanged(lse);
		}
	}
	
	/**
	 * A convenience reference to the main frame that logically owns 
	 * this menu helper. This instance variable provides a convenient
	 * reference to the main frame that logically owns this menu helper.
	 * Each instance of the main frame has its own set of menu 
	 * helpers.
	 */
	protected final MainFrame mainFrame;
	
	/**
	 * The list of main menu items created by this helper.
	 * 
	 * This list is used to update the status of menu entries and tool
	 * bar buttons based on the current context. Selected entries are 
	 * added to this list by derived classes when the getMenuItem() 
	 * or getTool() method is called.
	 */
	protected ArrayList<AbstractButton> contextItemList =
		new ArrayList<AbstractButton>();
	
	/**
	 * Reference to the JTree that last requested a selection listener.
	 * This reference is updated to refer to the last JTree that 
	 * requested a selection listener from this menu helper. This
	 * information is used to enable operation of the listener because
	 * currently in the event model it is impossible to obtain reference
	 * to the tree from the selection event.
	 */
	protected JTree tree = null;

	/**
	 * Reference to the JTable that last requested a selection listener.
	 * This reference is updated to refer to the last JTable that 
	 * requested a selection listener from this menu helper. This
	 * information is used to enable operation of the listener because
	 * currently in the event model it is impossible to obtain reference
	 * to the table from the selection event.
	 */
	protected JTable table = null;
	
	/**
	 * The type of helper that the derived class represents. This
	 * value is set when a menu helper object is instantiated and
	 * is never changed during the lifetime of the object.
	 */
	private final HelperType helperType;
}
