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

package org.peace_tools.decagon.configMgr;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.FileInfo;
import org.peace_tools.core.PEACEProperties;
import org.peace_tools.core.session.RemoteFileSystemView;
import org.peace_tools.core.session.RemoteServerSession;
import org.peace_tools.core.session.RemoteServerSessionMaker;
import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.decagon.PropertyKeys;
import org.peace_tools.generic.BorderWithComponent;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;

/**
 * Class to set information about shared storage space used by DECAGON.
 * 
 * This page permits the user to provide the information about the
 * shared storage space to be used by DECAGON. The shared storage
 * space can be on the local machine or on a remote machine. This
 * page provides the user with three different alternatives to specify
 * the server using a combo box provided by ServerInfoPanel in this
 * package. 
 * 
 * Once the server has been set the user can specify the shared drive
 * or path on the server. Note that this page checks to ensure that the
 * supplied path is valid and exists on the specified server.
 */
public class SharedStorageConfigPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a text box and
	 * button to specify the MST file (this file must not yet exist
	 * as we are going to generate it from the job), a combo box to
	 * select the server to run the job (this server must have peace
	 * already installed on it), JSpinners for number of nodes and
	 * CPUs per node to run the job. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param helpInfo A piece of HTML fragment loaded from the
	 * DecagonConfigInfo.html file. This information is extracted
	 * from the data file by DecagonConfigurationManager, wrapped
	 * in a suitable GUI component, and passed-in to this method.
	 */
	public SharedStorageConfigPage(DecagonConfigurationManager dcm, JComponent helpInfo) {
		assert(dcm != null);
		this.dcm = dcm;
		// Setup the title(s) for this page and border
		setTitle("Shared Storage", "Configure path to shared storage");
		setBorder(new EmptyBorder(5, 5, 5, 5));

		// Create the custom component o choose the type of server on
		// which the shared location is present.
		serverPanel = new ServerInfoPanel(true, true);
		
		// Create path-text-box and "Browse" button
		sharedPath = new JTextField(20);
		browse     = new JButton("Browse...");
		browse.setActionCommand("Browse");
		browse.addActionListener(this);
		// Create a panel to contain path-text box and "Browse" button.
		JPanel entryPane = new JPanel(new BorderLayout(5, 5));
		entryPane.add(sharedPath, BorderLayout.CENTER);
		entryPane.add(browse, BorderLayout.EAST);
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents("Set path to shared storage space", 
				"The server & path will be checked when the next button is clicked", 
				0, true, entryPane, Box.createVerticalStrut(5), helpInfo);

		// Create top-level panel to contain all the various fields.
		JPanel mainPanel = new JPanel(new BorderLayout(0, 5));
		mainPanel.add(serverPanel, BorderLayout.NORTH);
		mainPanel.add(subPanel, BorderLayout.SOUTH);
		
		// Create the top-level enable/disable check box
		enableSharedStorage = new JCheckBox("Enable Shared Storage", true);
		enableSharedStorage.setFont(enableSharedStorage.getFont().deriveFont(Font.BOLD));
		enableSharedStorage.setActionCommand("enableSharedStorage");
		enableSharedStorage.addActionListener(this);
		// Setup a custom border with the check-box part of the border
		BorderWithComponent border = new BorderWithComponent(enableSharedStorage, mainPanel, 
				BorderFactory.createCompoundBorder(BorderFactory.createEtchedBorder(),
						BorderFactory.createEmptyBorder(5, 10, 5, 10)));
		// Set up the border to make things look good.
		mainPanel.setBorder(border);
		// Add the contents to this page
		add(mainPanel, BorderLayout.CENTER);
		// Setup intial values
		setupInitalValues();
	}

	/**
	 * This method is invoked from the constructor to setup the default
	 * values for this page.
	 */
	private void setupInitalValues() {
		// Get the server panel to update its information
		serverPanel.updateServerInfo(PropertyKeys.STORAGE_SERVER_TYPE, 
				PropertyKeys.STORAGE_SERVER_KNOWN_HOST_KEY,
				PropertyKeys.STORAGE_SERVER_NAME, PropertyKeys.STORAGE_SERVER_PORT,
				PropertyKeys.STORAGE_SERVER_USERID);		
 		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		// Setup current shared server location (if set)
		String currSharedLoc = props.getProperty(PropertyKeys.STORAGE_SERVER_PATH);
		if (currSharedLoc != null) {
			sharedPath.setText(currSharedLoc);
		}
		// Obtain the enable/disable flag values
		final boolean enableFlag = props.getBooleanProperty(PropertyKeys.SHARED_SOTRAGE_ENABLED, true);
		enableSharedStorage.setSelected(enableFlag);
		enableInputs(enableFlag);
	}
	
	/**
	 * Helper method to enable/disable controls in this panel. 
	 * 
	 * This is a helper method that is used to enable or disable the
	 * various pertinent inputs in this wizard page.
	 * 
	 * @param enableFlag If this flag is true, then the controls are
	 * enabled. Otherwise they are disabled.
	 */
	private void enableInputs(boolean enableFlag) {
		serverPanel.enableInputs(enableFlag);
		sharedPath.setEnabled(enableFlag);
		browse.setEnabled(enableFlag);
	}
	
	/**
	 * Method to handle clicking of "Browse" button. This method essentially
	 * enables and disables the various inputs depending on the check box.
	 * 
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		final String cmd = event.getActionCommand();
		if (cmd.equals("Browse")) {
			try {
				final Server srvr = serverPanel.getServerEntry();
				JFileChooser jfc  = null;
				// Create a local or remote file system browser depending on 
				// local or remote server option.
				if (!srvr.isRemote()) {
					jfc = new JFileChooser();
				} else {
					// Get helper to establish remote server session in background
					// to ensure GUI remains responsive
					RemoteServerSessionMaker rssm = 
						new RemoteServerSessionMaker(dcm, srvr, 
								"Connecting to browse remote filesystem");
					rssm.start();
					RemoteServerSession rss = rssm.getSession();
					if (rss != null) {
						RemoteFileSystemView rfsv = new RemoteFileSystemView(rss);
						jfc = new JFileChooser(rfsv);
					} else {
						// Could not connect to server. User has already been 
						// notified of the problem.
						return;
					}
				} 
				jfc.setDialogTitle("Choose shared storage directory");
				jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

				if (jfc.showDialog(this, "Use Directory") == JFileChooser.APPROVE_OPTION) {
					// Copy the chosen directory to the mst file entry
					String wsPath = jfc.getSelectedFile().getAbsolutePath();
					// Set the selected item in this path.
					sharedPath.setText(wsPath);
				}
			} catch (Exception e) {
				ProgrammerLog.log(e);
				JPanel msg = Utilities.collapsedMessage("Error attempting to browse for shared storage location", 
						Utilities.toString(e));
				JOptionPane.showMessageDialog(this, msg, 
						"Unable to browse", JOptionPane.ERROR_MESSAGE);
			}
		} else if (cmd.equals("enableSharedStorage")) {
			// Enable/disable inputs.
			final boolean enabled = enableSharedStorage.isSelected();
			enableInputs(enabled);
		}
	}
	
	/**
	 * This method validates the shared directory is valid.
	 * 
	 * This method is invoked just before the user navigates to the
	 * adjacent page. This method checks to ensure that the shared
	 * storage path is valid before user navigates to the next page.
	 * 
	 * @param dialog The wizard dialog that owns this page. This parameter
	 * is currently unused.
	 * 
	 * @param currPage The zero-based logical index position of this page
	 * in the set of pages in the wizard. This parameter is not used 
	 * by this method.
	 * 
	 * @param nextPage The zero-based logical index of the subsequent page
	 * to which the user is attempting to switch. If this page is less than
	 * the currPage (the user clicked on "Prev" button) then this method
	 * does not perform any checks and freely permits the user to go back.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		if (!enableSharedStorage.isSelected()) {
			// The user has disabled the use of shared storage. So nothing
			// further to check below.
			return true;
		}
		try {
			// Ensure that the path/file is valid.
			final Server srvr = serverPanel.getServerEntry();
			// Create a server session validate path.
			ServerSession session = null;
			if (srvr.isRemote()) {
				// Create remote server session via background thread to ensure
				// the GUI remains responsive.
				RemoteServerSessionMaker rssm = 
					new RemoteServerSessionMaker(dcm, srvr, 
							"Verifying & validating shared storage directory");
				rssm.start();
				session = rssm.getSession();
			} else {
				session = SessionFactory.createSession(this, srvr);
				session.connect();
			}
			// Get the information about the shared storage location
			final FileInfo fi = session.fstat(sharedPath.getText());
			if (!fi.isDirectory() || !fi.canRead()) {
				// The shared location is not a directory or cannot be read.
				throw new IOException("The specified shared storage path is not a directory or it is not readable");			
			}
			JOptionPane.showMessageDialog(this, "The specified shared storage path successfully validated", 
					"Shared Storage is Good", JOptionPane.INFORMATION_MESSAGE);
		} catch (IOException e) {
			ProgrammerLog.log(e);
			JPanel msg = Utilities.collapsedMessage(INVALID_SHARED_SPACE_MSG, Utilities.toString(e));
			JOptionPane.showMessageDialog(this, msg, 
					"Invalid shared storage path", JOptionPane.ERROR_MESSAGE);			
			return false;
		}
		// The shared storage space is good.
		return true;
	}
	
	/**
	 * Returns summary information as a HTML fragment. 
	 * 
	 * This method is invoked by the top-level wizard for generating 
	 * summary information about the changes to be committed by this
	 * wizard page. 
	 * 
	 * @return A HTML-sub-fragment providing information about the
	 * configuration to be committed by this class.
	 */
	protected String getSummary() {
		StringBuilder summary = new StringBuilder(512);
		summary.append("<b>Use of Shared Storage : <i>" + 
				(enableSharedStorage.isSelected() ? "Enabled" : "Disabled") + "</i></b><br>");
		if (enableSharedStorage.isSelected()) {
			summary.append(serverPanel.getSummary());
			summary.append(Utilities.HTML_4SPACES);
			summary.append("Shared storage path: ");
			summary.append(sharedPath.getText());
			summary.append("<br>");
		} else {
			summary.append(Utilities.HTML_4SPACES);
			summary.append("Additional details are not applicable");
		}
		return summary.toString();
	}
	
	/**
	 * Convenience method to commit user-entered properties to the global
	 * properties. This method is invoked by the top-level Wizard once
	 * the user has verified the configuration.
	 */
	protected void commitProperties() {
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		// Setup the enable/disable flag values
		final boolean enabled = enableSharedStorage.isSelected();
		props.setProperty(PropertyKeys.SHARED_SOTRAGE_ENABLED, enabled);
		serverPanel.commitProperties(enabled, PropertyKeys.STORAGE_SERVER_TYPE, 
					PropertyKeys.STORAGE_SERVER_KNOWN_HOST_KEY,
					PropertyKeys.STORAGE_SERVER_NAME, PropertyKeys.STORAGE_SERVER_PORT,
					PropertyKeys.STORAGE_SERVER_USERID);
		// Setup path or clear the shared storage location
		if (enabled) {
			props.setProperty(PropertyKeys.STORAGE_SERVER_PATH,	sharedPath.getText());
		} else {
			props.removeProperty(PropertyKeys.STORAGE_SERVER_PATH);
		}
	}
	
	/**
	 * The top-level check-box that enables or disables the items in
	 * this wizard-page. If the check-box is checked, then other GUI
	 * items are enabled and permit the user to configure path to shared
	 * storage location. Otherwise the use of shared location is 
	 * disable.
	 */
	private JCheckBox enableSharedStorage;
	
	/**
	 * A custom component in this package that permits the user to conveniently
	 * select the server where the shared storage is located. The custom
	 * component provides three different mechanisms to specify the location
	 * of the shared storage space.
	 */
	private final ServerInfoPanel serverPanel;
	
	/**
	 * Field to read/display the path to the shared storage space.
	 */
	private JTextField sharedPath;

	/**
	 * The browse button to be enabled the user to choose the
	 * directory where the file is to be stored.
	 */
	private JButton browse;

	/**
	 * Reference to the wizard that actually owns this page. This is used
	 * when creating sub-dialogs etc. in this page.
	 */
	private final DecagonConfigurationManager dcm;

	/**
	 * A simple message that is displayed to the user if the shared
	 * storage space appears to be invalid.
	 */
	private static final String INVALID_SHARED_SPACE_MSG = 
		"<html>" +
		"The specified shared storage path is invalid.<br>" +
		"Please verify that the following requirements are met:<ul>" +
		"<li>Ensure that server information (hostname/IP, port) etc. are valid</li>"+
		"<li>Ensure the path is a directory and you have read privileges</li>" +
		"</ul>Please appropriately edit the necessary information."+
		"</html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1193172876744546755L;
}
