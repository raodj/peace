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

import javax.swing.Box;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JToolBar;

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
public class HelpMenuHelper implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "Help" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public HelpMenuHelper(MainFrame mainFrame) {
		this.mainFrame = mainFrame;
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
		JMenuItem item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Help Contents",
				"Open the main help and documentation page on PEACE website",
				"HelpContents", this, "images/24x24/HelpContents.png", null, true, false);
		helpMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, "User Tutorials",
					"Open PEACE web site with tutorials on using and working with PEACE",
					"HelpTutorials", this, "images/24x24/Help.png", 
					null, true, false);
		helpMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Related Publications",
				"Browse the list of publications and presentations on PEACE",
				"HelpPubs", this, "images/24x24/Papers.png", null, true, false);
		helpMenu.add(item);
		helpMenu.addSeparator();
		//----------------------------------------------------------
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Contact via Email",
				"Launches your default mail tool to send email to PEACE developers",
				"HelpEmail", this, "images/24x24/Email.png", null, true, false);
		helpMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Software Updates",
				"Navigates directly to the downloads page on PEACE website",
				"HelpUpdates", this, "images/24x24/Updates.png", null, true, false);
		helpMenu.add(item);
		helpMenu.addSeparator();
		//----------------------------------------------------------
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Programmer Documentation",
				"Opens PEACE web site with technical documentation about PEACE",
				"HelpProgDoc", this, "images/24x24/ProgDoc.png", 
				null, true, false);
		helpMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Release Notes",
				"Opens PEACE web site with information about current release",
				"HelpNotes", this, "images/24x24/ReleaseNotes.png", 
				null, true, false);
		helpMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Report Bugs & Issues",
				"Opens PEACE web site with forms to log bugs and issues with PEACE",
				"HelpBug", this, "images/24x24/Bug.png", null, true, false);
		helpMenu.add(item);
		helpMenu.addSeparator();
		//----------------------------------------------------------
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Welcome & Quick Start",
				"The initial welcome screen with quick start guide about PEACE",
				"Welcome", this, "images/24x24/WelcomeSrc.png", 
				null, true, false);
		helpMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"About PEACE",
				"Credits and copyright information about the PEACE software system",
				"About", this, "images/24x24/PEACE.png", null, true, false);
		helpMenu.add(item);
		
		if (toolbar != null) {
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(Utilities.createToolButton("images/24x24/WelcomeSrc.png", 
					null, "Welcome", this, 
					"Show welcome message & quick start guide", true));
			toolbar.add(Utilities.createToolButton("images/24x24/Help.png", 
					null, "HelpContents", this, 
					"Launch online PEACE via your default browser", true));
		}
		return helpMenu;
	}
	
	@Override
	public void actionPerformed(ActionEvent event) {
		// These URLs are listed in the same order as the action commands
		// to see look up and cross reference.
		String HelpURLs[] = {"contents.html", "tutorials.html",
				"papers.html", "updates.html", "progdocs.html", 
				"release_notes.html", "bugreport.html", "contact.html"};
		// The action commands.
		String ActionCmds[] = {"HelpContents", "HelpTutorials", 
				"HelpPubs", "HelpUpdates", "HelpProgDoc",
				"HelpNotes", "HelpBug", "HelpEmail", "Welcome", "About"};
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
			vf.createView("../../../installFiles/welcome.html", null,
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

	/**
	 * This message is displayed when the default system browser could
	 * not be launched to display help.
	 */
	private static final String ABOUT_ERROR = "<html>" +
		"Unable to create the About dialog box with the necessary information.<br>" +
		"Possibly your install of PEACE is incomplete or corrupted.<br>" +
		"You will need to reinstall a recent copy of PEACE to resolve this issue." +
		"</html>";
	
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
}
