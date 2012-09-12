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

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JTextField;

/**
 * A convenience class to prompt the user to enter password.
 * 
 * <p>This is a customized component that is used to prompt the user to
 * enter password for privileged operations. This class provides a 
 * standard dialog box with user ID and password fields. Typically,
 * the user ID field is un-editable. This generic class is used
 * from couple of different spots in PEACE to obtain password 
 * information from the user.</p>
 * 
 * <p>This component can be used in the following manner:
 * 
 * <code>
 * 
 * </code>
 * </p>
 */
public class PasswordDialog {

	/**
	 * The GUI component that is used to display the user ID value.
	 * Typically this is left immutable unless the user desires
	 * to make it mutable.
	 */
	private final JTextField userID;
	
	/**
	 * The GUI component that is used to permit the user to enter
	 * the password.
	 */
	private final JTextField password;
	
	/**
	 * The top-level panel that contains the various GUI components
	 * used by this dialog box.
	 */
	private final JPanel topPanel;
	
	/**
	 * The only constructor for this class.
	 * 
	 * @param userName The default user ID to be set in the userID GUI
	 * component created by the constructor. This value cannot be null.
	 * 
	 * @param isUserNameEditable If this flag is true then the user ID
	 * can be edited by the user. If it is false, then the user cannot
	 * edit the user ID value.
	 * 
	 * @param serverName The name of the server for which the password
	 * is desired. This value cannot be null.
	 * 
	 * @param purpose An optional purpose for establishing the connection.
	 * The purpose provides additional information to the user. This
	 * value can be null.
	 */
	public PasswordDialog(final String userName, 
			final boolean isUserNameEditable, final String serverName, 
			final String purpose) {
		// Create text fields for user name and password.
		userID = new JTextField(10);
		userID.setEditable(isUserNameEditable);
		userID.setText(userName);
		// Create the special password field.
		password = new JPasswordField(10);

		// Create components by laying them out appropriately
		JPanel credPanel = new JPanel(new GridLayout(2, 2, 0, 3));
		credPanel.add(new JLabel("User id:"));
		credPanel.add(userID);
		Utilities.adjustDimension(userID, 0, 6);
		credPanel.add(new JLabel("Password:"));
		Utilities.adjustDimension(password, 0, 6);
		credPanel.add(password);
		// Another panel to control the size of the grid layout 
		// to ensure it looks decent with a message.
		JPanel msgPanel = new JPanel(new BorderLayout(0, 5));
		msgPanel.add(credPanel, BorderLayout.SOUTH);
		// Add a label indicating server information.
		JLabel subInfo = new JLabel("<html>Enter credentials for <i>" + 
				serverName + "</i></html>");
		msgPanel.add(subInfo, BorderLayout.NORTH);
		// If purpose has been given add purpose information into another
		// panel.
		if (purpose != null) {
			JLabel info = new JLabel(purpose, Utilities.getIcon("images/32x32/Information.png"), 
					JLabel.LEFT);
			JPanel outer = new JPanel(new BorderLayout(0, 10));
			outer.add(info, BorderLayout.CENTER);
			outer.add(msgPanel, BorderLayout.SOUTH);
			// Set message panel to be the outer most one now.
			msgPanel = outer;
		}
		// Update the top-panel with the main container.
		topPanel = msgPanel;
	}
	
	/**
	 * Actually display the dialog box to prompt the user for password.
	 * 
	 * This method must be used to actually prompt the user to enter 
	 * the password.
	 * 
	 * @param parent The parent component (if any) to be used to prompt the
	 * user to enter the password. This value can be null.
	 * 
	 * @return This method returns true if the user pressed the "OK" button.
	 * Otherwise this method returns false.
	 */
	public boolean prompt(Component parent) {
		// Pack all the elements into an array
		Object items[] = { topPanel };
		int result = 
			JOptionPane.showConfirmDialog(parent, items, "Enter Password", 
					JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);
		return (result == JOptionPane.YES_OPTION); 
	}
	
	/**
	 * Obtain the user ID value set by the user.
	 * 
	 * @return The user ID value edited by the user. This method never
	 * returns a null string.
	 */
	public String getUserID() {
		return userID.getText();
	}
	
	/**
	 * Obtain the password entered by the user.
	 * 
	 * @return The password entered by the user. This method never returns
	 * a null string.
	 */
	public String getPassword() {
		return password.getText();
	}
}