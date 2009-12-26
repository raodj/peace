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

package org.peace_tools.core;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.AbstractButton;
import javax.swing.Box;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TreeSelectionListener;

import org.peace_tools.generic.Utilities;

/**
 * The help menu helper for the application.
 * 
 * This class encapsulates the code related to the operations performed by
 * various menu items in the "Help" menu. This class is typically created once
 * from the MainFrame.createMenus() method. The primary motivation for
 * introducing sub-menu handler classes is to improve code organization and
 * minimize code clutter. Note that the HelpMenuHelper is essentially an event
 * handler that is set on the various menu items in this class. This helper
 * class also provides a  createHelpMenu() method that actually creates
 * the "Help" menu.
 */
public class HelpMenuHelper extends AbstractMenuHelper 
implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "Help" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public HelpMenuHelper(MainFrame mainFrame) {
		super(HelperType.HELP_MENU, mainFrame);
	}

	/**
	 * Helper method to create the help menu.
	 * 
	 * This method performs the actual task of creating the "Help" menu. This
	 * method has been introduced to organize all the methods related to the
	 * "Help" menu into a single class. This method is invoked from the MainMenu
	 * class to create the "Help" menu.
	 * 
	 * @param toolbar The toolbar to which frequently used shortcuts can
	 * be added typically in the form of icons. If the toolbar is null, then
	 * shortcuts are not added.
	 */
	public JMenu createHelpMenu(JToolBar toolbar) {
		// Create the actual menu.
		JMenu helpMenu = new JMenu("Help");
	    // First create and add the new file creation menu options.
		helpMenu.add(getMenuItem(ActionType.SHOW_HELP_CONTENTS, true));
		helpMenu.add(getMenuItem(ActionType.SHOW_USER_TUTORIALS, true));
		helpMenu.add(getMenuItem(ActionType.SHOW_PAPERS, true));
		helpMenu.addSeparator();
		
		//----------------------------------------------------------
		helpMenu.add(getMenuItem(ActionType.SHOW_EMAIL, true));
		helpMenu.add(getMenuItem(ActionType.SHOW_UPDATES, true));
		helpMenu.addSeparator();
		//----------------------------------------------------------
		helpMenu.add(getMenuItem(ActionType.SHOW_PROGRAMMER_DOCS, true));
		helpMenu.add(getMenuItem(ActionType.SHOW_NOTES, true));
		helpMenu.add(getMenuItem(ActionType.SHOW_BUGS, true));
		helpMenu.addSeparator();
		//----------------------------------------------------------
		helpMenu.add(getMenuItem(ActionType.SHOW_QUICK_START, true));
		helpMenu.add(getMenuItem(ActionType.SHOW_ABOUT_DIALOG, true));
		
		if (toolbar != null) {
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(getTool(ActionType.SHOW_QUICK_START, true));
			toolbar.add(getTool(ActionType.SHOW_HELP_CONTENTS, true));	
		}
		return helpMenu;
	}
	
	@Override
	public void actionPerformed(ActionEvent event) {
		// These URLs are listed in the same order as the action commands
		// to see look up and cross reference.
		String HelpURLs[] = {"contents.html", "tutorials.html",
				"papers.html", "contact.html", "updates.html", "progdocs.html", 
				"release_notes.html", "bugreport.html"};
		// Map action command to index ids
		int index;
		final String cmd = event.getActionCommand();
		for(index = 0; (index < ActionCmds.length); index++) {
			if (ActionCmds[index].equals(cmd)) {
				// Found a matching command
				break;
			}
		}
		// Do different actions based on the action chosen.
		if (index < HelpURLs.length) {
			// This is online HTML page request. Use helper method to
			// launch browser if possible.
			mainFrame.showHelp("http://www.peace-tools.org/" + HelpURLs[index]);
		} else if (index == ActionCmds.length - 2) {
			// Show the welcome screen if it is not already there.
			ViewFactory vf = mainFrame.getViewFactory();
			vf.createView("installFiles/welcome.html", null,
						  ViewFactory.ViewType.HTML_VIEW, false, false);
		} else if (index == ActionCmds.length - 1) {
			// Launch about dialog box.
			try {
				AboutDialog ad = new AboutDialog(mainFrame);
				ad.pack();
				ad.setLocationRelativeTo(mainFrame); 
				ad.setVisible(true);
			} catch (Exception e) {
				JPanel msg = Utilities.collapsedMessage(ABOUT_ERROR,
						Utilities.toString(e));
				JOptionPane.showMessageDialog(mainFrame, msg,
						"Error Displaying About Box", JOptionPane.ERROR_MESSAGE);
			}
		}
	}

	@Override
	public ActionListener getActionListener() {
		return this;
	}

	@Override
	public ListSelectionListener getListSelectionListener(JTable table) {
		return null;
	}

	@Override
	public JMenuItem getMenuItem(ActionType actionType, boolean mainMenu) {
		int index = actionType.ordinal() - ActionType.SHOW_HELP_CONTENTS.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainMenu ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";
		// Create and return the main menu item
		return Utilities.createMenuItem(Utilities.MENU_ITEM, MenuTitles[index],
				(mainMenu ? MenuSubTitles[index] : null),
				ActionCmds[index], this, IconPath, 
				null, true, false);
	}
	
	@Override
	public AbstractButton getTool(ActionType actionType, boolean mainToolBar) {
		int index = actionType.ordinal() - ActionType.SHOW_HELP_CONTENTS.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainToolBar ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";

		return Utilities.createToolButton(IconPath, null, ActionCmds[index], this, 
				MenuSubTitles[index], true);
	}

	@Override
	public TreeSelectionListener getTreeSelectionListener(JTree tree) {
		// TODO Auto-generated method stub
		return null;
	}
	
	/**
	 * The strings for each menu title created by this helper. The list of
	 * values are organized in the same order as the ordinal values of 
	 * the ActionType enumeration. Preserving the order is important.
	 */
	private static final String MenuTitles[] = {
		"Help Contents", 
		"User Tutorials",
		"Related Publications",
		"Contact via Email",
		"Software Updates",
		"Programmer Documentation",
		"Release Notes",
		"Report Bugs & Issues",
		"Welcome & Quick Start Guide",
		"About PEACE"
	};
	
	/**
	 * The icon file names for each menu title created by this helper. 
	 * The list of values are organized in the same order as the ordinal 
	 * values of the ActionType enumeration. Preserving the order is
	 * important. Note that a prefix directory (such as: images/16x16)
	 * and a suffix extension (.png) is added when tools or menu items
	 * are created.
	 */
	private static final String IconNames[] = {
		"HelpContents", 
		"Help",
		"Papers",
		"Email",
		"Updates",
		"ProgDoc",
		"ReleaseNotes",
		"Bug",
		"WelcomeSrc",
		"PEACE"
	};
	
	/**
	 * The strings for the action commands generated by the various menu items
	 * in the help menu. The list of values are organized in the same order as
	 * the ordinal values of the ActionType enumeration. Preserving the order 
	 * is important.
	 */
	private static final String ActionCmds[] = {
		"HelpContents", 
		"HelpTutorials",
		"HelpPubs",
		"HelpEmail",
		"HelpUpdates",
		"HelpProgDoc",
		"HelpNotes",
		"HelpBug",
		"Welcome",
		"About"
	};
	
	/**
	 * The various sub menu titles that are used in the main menu. The
	 * sub menu titles are used to provide the user with a bit more 
	 * verbose description on the action that will be performed by a given
	 * menu item. The list of values are organized in the same order as 
	 *  the ordinal values of the ActionType enumeration. Preserving the 
	 *  order is important.
	 */
	private static final String MenuSubTitles[] = {
		"Open the main help and documentation page on PEACE website",
		"Open PEACE web site with tutorials on using and working with PEACE",
		"Browse the list of publications and presentations on PEACE",
		"Launches your default mail tool to send email to PEACE developers",
		"Navigates directly to the downloads page on PEACE website",
		"Opens PEACE web site with technical documentation about PEACE",
		"Opens PEACE web site with information about current release",
		"Opens PEACE web site with forms to log bugs and issues with PEACE",
		"The initial welcome screen with quick start guide about PEACE",
		"Credits and copyright information about the PEACE software system"
	};
	
	/**
	 * This message is displayed when the default system browser could
	 * not be launched to display help.
	 */
	private static final String ABOUT_ERROR = "<html>" +
		"Unable to create the About dialog box with the necessary information.<br>" +
		"Possibly your install of PEACE is incomplete or corrupted.<br>" +
		"You will need to reinstall a recent copy of PEACE to resolve this issue." +
		"</html>";
}
