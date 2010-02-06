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
import java.awt.Desktop;
import java.net.URI;

import javax.swing.JOptionPane;
import javax.swing.JPanel;

/** 
 * A simple class to handle display of help URLs. 
 * 
 * This class was introduced to centralize the generic process of 
 * displaying help information to the user. This class contains
 * static methods that can be directly used to display help from
 * peace-tools.org via the default browser.
 */
public class HelpHandler {	
	/** Displays a given URL.
	 * 
	 * This method is used when the user chooses to view one of
	 * the selected help topics.  This method is actually invoked 
	 * from various classes to display help. The help is actually
	 * redirected to the PEACE web site.  Having users visit a web
	 * site permits help content to be developed independently
	 * and updated consistently.
	 * 
	 * @param parent The parent component that must be used to display
	 * any error dialogs.
	 * 
	 * @param url The complete web site URL that must be displayed in 
	 * the default system browser.
	 */
	public static void showHelp(Component parent, String url) {
		if (Desktop.isDesktopSupported() && 
				(Desktop.getDesktop().isSupported(Desktop.Action.BROWSE))) {
			// Desktop is supported with browse to launch browser
			try {
				Desktop.getDesktop().browse(new URI(url));
			} catch (Exception e) {
				// Log the exception
				ProgrammerLog.log(e);
				// Display error to the user.
				UserLog.log(UserLog.LogLevel.ERROR, 
						"PEACE", e.getMessage());
				JPanel msg = Utilities.collapsedMessage(HELP_ERROR,
						Utilities.toString(e));
				JOptionPane.showMessageDialog(parent, msg,
						"Error Displaying Help", JOptionPane.ERROR_MESSAGE);
			}
		} else {
			JOptionPane.showMessageDialog(parent, 
					"Your Java version does not provide the necessary features\n" +
					"to launch the online help directly from PEACE GUI.\n" +
					"You may view the requested help information via the URL:\n" +
					url, "Java Feature Unavailable", JOptionPane.WARNING_MESSAGE);
		}
	}
	
	/**
	 * This message is displayed when the default system browser could
	 * not be launched to display help.
	 */
	private static final String HELP_ERROR = "<html>" +
		"The default system browser could not be launched to display help.<br>" +
		"You may view the help web site directly from your system browser by<br>" +
		"navigating to the appropriate help section from http://www.peace-tools.org/" +
		"</html>";
	
	/**
	 * Private constructor to ensure this class is never instantiated.
	 * 
	 * All the methods in this class are static and are designed to be
	 * used directly. In other words, this class is not meant to be
	 * instantiated. Consequently, the constructor is private.
	 */
	private HelpHandler() {}
}
