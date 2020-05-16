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

package org.peace_tools.core.server;

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

/**
 * This class serves as the last page in a ServerWizard. This page gives a
 * summary information to the user and permits the user the last chance to
 * cancel out of the server addition process. If the user proceeds, then this
 * method adds a new server entry to the workspace. Next it creates a new
 * PEACEInstaller tab and adds it to the center page. The PEACE installer
 * performs the actual task of installation using the newly added server entry.
 * An independent PEACE installer is used because the installation process takes
 * a while and the user can proceed with other operations rather than waiting
 * for the installation to complete.
 * 
 * <p>
 * <b>Note:</b> This page does the same operation for both local and remote
 * servers.
 * </p>
 */
public class ServerAddingEntryWizardPage extends GenericWizardPage {
	/**
	 * The constructor. The constructor sets up the various components on this
	 * wizard page. The components include a set of labels that display information
	 * to the user.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * @param server The in-memory server object that contains the core information
	 *               about the server.
	 */
	public ServerAddingEntryWizardPage(ServerWizard wizard) {
		this.wizard = wizard;

		// Setup the title(s) for this page and border
		setTitle("Server Summary", "Adds server entry & starts installation");
		setBorder(new EmptyBorder(10, 10, 5, 10));

		// Create labels with all the necessary information.
		for (int i = 0; (i < infoFields.length); i++) {
			infoFields[i] = new JTextField();
			infoFields[i].setEditable(false);
			infoFields[i].setBackground(Color.lightGray.brighter());
		}

		// Pack all the labels into a panel with vertical box layout.
		JPanel container = Utilities.createLabeledComponents(STANDARD_MSG, null, 4, true, new JLabel(" "),
				new JLabel("Server's host name:"), infoFields[0], new JLabel(" "), // Blank line
				new JLabel("Install directory:"), infoFields[1], new JLabel(" "), // Blank line
				new JLabel(INFO_MSG, Utilities.getIcon("images/32x32/Information.png"), SwingConstants.LEFT));
		add(container, BorderLayout.CENTER);
	}

	/**
	 * This method is called just before this page is to be displayed. This page
	 * essentially updates the data being displayed in the GUI fields from the data
	 * stored in the in-memory Server object.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// First update the necessary information.
		infoFields[0].setText(wizard.getServer().getName());
		infoFields[1].setText(wizard.getServer().getInstallPath());
	}

	/**
	 * A reference to the wizard dialog that logically owns this page. This
	 * reference is used to enable and disable buttons on this wizard appropriately.
	 */
	private final ServerWizard wizard;

	/**
	 * Field to be displayed to the user in this wizard pane.
	 */
	private JTextField[] infoFields = new JTextField[2];

	/**
	 * A simple standard message to be displayed at the top of this page.
	 */
	private static final String STANDARD_MSG = "<html>You have provided all the necessary information to<br>"
			+ "create a server entry. On clicking the Finish button, a new server<br>"
			+ "entry will be added to the workspace.<br><br>"
			+ "Next, the actual PEACE runtime installation will commence.<br>"
			+ "The installation process will take several (about 5 to 10)<br>"
			+ "minutes depending on the server speed and network bandwidth.</html>";

	/**
	 * A simple standard message to be displayed at the bottom of this page.
	 */
	private static final String INFO_MSG = "<html>The install process will be run in the background<br>"
			+ "and the results will be displayed on a separate<br>window in the workpsace.</html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
