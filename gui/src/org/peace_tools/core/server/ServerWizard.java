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
import javax.swing.JOptionPane;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.MainFrame;
import org.peace_tools.core.PEACEInstaller;
import org.peace_tools.core.session.ServerSession;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.generic.dndTabs.DnDTabbedPane;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as the top-level class for adding a new server entry to the
 * work space. This top-level class merely creates the various pages and adds
 * them to the wizard. Each page performs a specific task required to create a
 * complete Server
 */
public class ServerWizard extends WizardDialog {
	/**
	 * The constructor. The constructor sets up the various pages in this wizard.
	 * Most of the pages are independent classes that are developed separately and
	 * tied together in this wizard.
	 * 
	 * @param title      The title to be set for this wizard.
	 * @param parent     The parent frame that logically owns this wizard.
	 * @param tabbedPane The permanent, center tabbed pane to which the actual
	 *                   installer pane must be added if this wizard completes
	 *                   successfully.
	 */
	public ServerWizard(String title, MainFrame parent, DnDTabbedPane tabbedPane) {
		super(parent);

		// Save reference to the tabbed pane
		this.tabbedPane = tabbedPane;

		setTitle(title);
		setResizable(false);
		// Set up the title image we want to use.
		setTitleBackground("images/peace_wizard_header.png", Color.white);
		// Set up the column image we want to use.
		setSequenceBackground("images/peace_wizard_column.png");

		// Create a dummy server entry that is being edited by this wizard
		server = new Server("", "", "", "", "", Server.getDefaultPollDuration(), true, -1);

		// First setup the overview page.
		createOverview();

		// Create the server type selection page.
		ServerTypeWizardPage stwp = new ServerTypeWizardPage(this, false);
		addPage(stwp);

		// Create the server information entry page.
		ServerInfoWizardPage siwz = new ServerInfoWizardPage(this, false);
		addPage(siwz);

		// Create the software component selection page.
		ServerComponentsSelectionWizardPage scswp = new ServerComponentsSelectionWizardPage(this);
		addPage(scswp);

		// Create the final summary page.
		ServerAddingEntryWizardPage saewp = new ServerAddingEntryWizardPage(this);
		addPage(saewp);
	}

	/**
	 * Helper method to create the overview page. This method was introduced to keep
	 * the code clutter in the constructor to a bare minimum.
	 */
	private void createOverview() {
		JLabel message = new JLabel(OVERVIEW_MSG);
		GenericWizardPage overview = new GenericWizardPage();
		overview.add(message, BorderLayout.NORTH);
		overview.setTitle("Overview", "Overview of tasks in this wizard.");
		overview.setBorder(new EmptyBorder(30, 15, 10, 5));
		addPage(overview);
	}

	/**
	 * Helper method invoked when user clicks cancel button.
	 * 
	 * This is a helper method that is overridden in this class. This method is
	 * invoked when the user clicks the cancel button in the wizard. This method is
	 * used to display a confirmation dialog to ensure that the user really wants to
	 * exit out of the wizard.
	 * 
	 * @return This method returns true if the user wants to quit out of the wizard
	 *         dialog.
	 */
	@Override
	protected boolean cancel() {
		// Check if user really want's to quit.
		int result = JOptionPane.showConfirmDialog(this, "Are you sure you want to exit from this wizard?", "Confirm",
				JOptionPane.YES_NO_OPTION);

		// The user does not want to quit.
		if (result == JOptionPane.NO_OPTION) {
			return false;
		}

		// Yes, the user wants to quit. Ensure no threads are left.
		return super.cancel();
	}

	/**
	 * Obtain the current server entry that is being edited via this server wizard.
	 */
	public Server getServer() {
		return server;
	}

	/**
	 * Obtain the server session that is active within this wizard.
	 */
	public ServerSession getServerSession() {
		return serverSession;
	}

	/**
	 * Update the server session of this wizard.
	 */
	public void setServerSession(ServerSession serverSession) {
		this.serverSession = serverSession;
	}

	/**
	 * This method overrides the final notification method in this class to launch
	 * the actual background installer thread if the wizard successfully completed.
	 */
	@Override
	public void done(boolean success) {
		if (!success) {
			if (serverSession != null) {
				serverSession.disconnect();
			}
			return;
		}

		// First add a new server entry to the global workspace.
		Workspace.get().getServerList().add(server);

		// Save the work space.
		MainFrame mf = (MainFrame) getParent();
		mf.saveWorkspace(null);

		// Create installer and add it to the main frame.
		PEACEInstaller installer = new PEACEInstaller(server, serverSession);
		tabbedPane.createSplitPane("PEACE Installer", Utilities.getIcon("images/16x16/ServerInstalling.png"),
				installer.getPanel(), DnDTabbedPane.Location.CENTER);

		// Now start the installer as a background thread.
		installer.execute();
	}

	/**
	 * Determine if g++ has been detected on the target server. This method is used
	 * by the {@link ServerComponentsSelectionWizardPage} to permit the user to
	 * enable/disable EAST installation.
	 * 
	 * @return This method returns true if g++ is available on the target server.
	 *         Otherwise it returns false.
	 */
	public boolean haveGCC() {
		return gccFlag;
	}

	/**
	 * Set flag to indicate if g++ was detected on the target server. The flag is
	 * set by the {@link ServerInfoWizardPage#checkGCC(String[])} method if g++ was
	 * detected on the target machine.
	 * 
	 * @param flag If this flag is true, that indicates g++ was successfully
	 *             detected.
	 */
	public void setHaveGCC(final boolean flag) {
		this.gccFlag = flag;
	}

	/**
	 * Flag that is used to track if the server on which EAST is to be installed has
	 * g++ on it. The availability of g++ is detected by the
	 * {@link ServerInfoWizardPage#checkGCC(String[])} method which then sets this
	 * flag. The flag is used by the {@link ServerComponentsSelectionWizardPage}
	 * class.
	 */
	private boolean gccFlag;

	/**
	 * This is the server entry that is being currently edited via this wizard. This
	 * object is accessed by the various pages in this wizard and updated as the
	 * user navigates through the wizard pages.
	 */
	private final Server server;

	/**
	 * This is the server session that is active within this wizard. This object is
	 * accessed by the various pages in this wizard and updated as the user
	 * navigates through the wizard pages.
	 */
	private ServerSession serverSession;

	/**
	 * The DnDTabbedPane to which the actual PEACE installer tab must be added. This
	 * value is set in the constructor.
	 */
	private final DnDTabbedPane tabbedPane;

	/**
	 * A static overview message that is displayed in the first overview page
	 * displayed by this wizard to the user.
	 */
	private static final String OVERVIEW_MSG = "<html>" + "This wizard guides you through the process of installing<br>"
			+ "the PEACE computational backend software system on the<br>"
			+ "local computer (aka server) or a remote server. Currently,<br>"
			+ "only Linux servers are supported for remote machines.<br><br>"
			+ "For accessing remote computational clusters, this wizard will<br>"
			+ "prompt you to enter your credentials (user name and password).<br>"
			+ "The user name is persisted in the workspace configuration but<br>"
			+ "<i>the password is not persisted for security reasons</i>. Therefore,<br>"
			+ "you will be prompted to enter password(s) once each time PEACE<br>" + "is run.<br><br>"
			+ "<b>Note</b>: The process of installing PEACE backend on a Linux<br>"
			+ "server takes 3-5 minutes depending on your network,<br>"
			+ "speed of remote server, and load on the machine.<br></html>";

	/**
	 * This is a generated serial version ID.
	 */
	private static final long serialVersionUID = -369000313738268902L;
}
