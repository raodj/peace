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
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.PEACEProperties;
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
 * Class to configure a SSH-tunnel to access a MySQL data base.
 * 
 * This page permits the user to configure an optional SSH-tunnel
 * for accessing the MySQL data base used by DECAGON for collating
 * results from various experiments. The SSH-tunnel is typically
 * used when the MySQL data base is behind a firewall on a machine
 * for added security. In such scenarios, first an SSH connection
 * is made to a given server, and upon successful authentication
 * a tunnel from the local machine to the MySQL database-port is
 * setup over the SSH-connection. The secure tunnel is then used
 * to access and interact with MySQL server. Note that the SSH-tunnel
 * is optional.
 */
public class MySQLTunnelConfigPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include a
	 * ServerInfoPanel object that permits the user to select the
	 * server to be used for establishing the SSH-tunnel.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param helpInfo A piece of HTML fragment loaded from the
	 * DecagonConfigInfo.html file. This information is extracted
	 * from the data file by DecagonConfigurationManager, wrapped
	 * in a suitable GUI component, and passed-in to this method.
	 */
	public MySQLTunnelConfigPage(DecagonConfigurationManager dcm, 
			JComponent helpInfo) {
		assert(dcm != null);
		this.dcm = dcm;
		// Setup the title(s) for this page and border
		setTitle("MySQL SSH Tunnel", "Configure SSH-tunnel to access MySQL database");
		setBorder(new EmptyBorder(5, 5, 5, 5));

		// Create the custom component to choose the server to be used
		// for SSH-tunnel.
		serverPanel = new ServerInfoPanel(false, true);
		// Create an informational label to provide summary information
		// to the user.
		JLabel info = new JLabel(SSH_TUNNEL_INFO_MSG, 
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		// Create top-level panel to contain the various fields
		JPanel mainPanel = new JPanel(new BorderLayout(0, 5));
		mainPanel.add(info, BorderLayout.NORTH);
		mainPanel.add(serverPanel, BorderLayout.CENTER);
		mainPanel.add(helpInfo, BorderLayout.SOUTH);
		
		// Create the top-level enable/disable check box
		enableSSHTunnel = new JCheckBox("Use SSH Tunnel", true);
		enableSSHTunnel.setFont(enableSSHTunnel.getFont().deriveFont(Font.BOLD));
		enableSSHTunnel.setActionCommand("enableSSHTunnel");
		enableSSHTunnel.addActionListener(this);
		// Setup a custom border with the check-box part of the border
		BorderWithComponent border = new BorderWithComponent(enableSSHTunnel, mainPanel, 
				BorderFactory.createCompoundBorder(BorderFactory.createEtchedBorder(),
						BorderFactory.createEmptyBorder(5, 10, 5, 10)));
		// Set up the border to make things look good.
		mainPanel.setBorder(border);
		// Add the contents to this page
		add(mainPanel, BorderLayout.CENTER);
		// Setup initial values
		setupInitialValues();
	}

	/**
	 * This method is invoked from the constructor to setup the default
	 * values for this page.
	 */
	private void setupInitialValues() {
		// Get the server panel to update its information
		serverPanel.updateServerInfo(PropertyKeys.MYSQL_TUNNEL_SERVER_TYPE,
				PropertyKeys.MYSQL_TUNNEL_SERVER_KNOWN_HOST_KEY,
				PropertyKeys.MYSQL_TUNNEL_SERVER_NAME, PropertyKeys.MYSQL_TUNNEL_SERVER_PORT,
				PropertyKeys.MYSQL_TUNNEL_SERVER_USERID);
		
 		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		// Obtain the enable/disable flag values
		final boolean enableFlag = props.getBooleanProperty(PropertyKeys.MYSQL_TUNNEL_ENABLED, true);
		serverPanel.enableInputs(enableFlag);
		enableSSHTunnel.setSelected(enableFlag);
	}
		
	/**
	 * Method to handle clicking of "Browse" button. This method essentially
	 * enables and disables the various inputs depending on the check box.
	 * 
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		// Enable/disable inputs.
		final boolean enabled = enableSSHTunnel.isSelected();
		serverPanel.enableInputs(enabled);
	}
	
	/**
	 * Convenience method to obtain the Server entry to be used to create 
	 * a SSH-tunnel.
	 * 
	 * This method returns the Server entry configured by the user to 
	 * establish a SSH-tunnel to access a MySQL data base. 
	 * 
	 * @return The server entry to be used to setup SSH-tunnel. If the
	 * user has disabled tunneling then this method returns null.
	 */
	public Server getTunnelServer() {
		return (!enableSSHTunnel.isSelected() ? null : 
			serverPanel.getServerEntry());
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
		if (!enableSSHTunnel.isSelected()) {
			// The user has disabled SSH tunnel. No need to check connection
			return true;
		}
		try {
			// Ensure that the path/file is valid.
			final Server srvr = serverPanel.getServerEntry();
			// Create a server session validate login to the server.
			ServerSession session = null;
			if (srvr.isRemote()) {
				// Create remote server session via background thread to ensure
				// the GUI remains responsive.
				RemoteServerSessionMaker rssm = 
					new RemoteServerSessionMaker(dcm, srvr, 
							"Verifying & validating SSH-tunnel server");
				rssm.start();
				session = rssm.getSession();
				if (session != null) {
					session.disconnect();
				} else {
					// Error connecting via SSH. Show message.
					JOptionPane.showMessageDialog(this, SSH_TUNNEL_ERROR_MSG, 
							"Invalid SSH-Tunnel Configuration", JOptionPane.ERROR_MESSAGE);			
					return false;
				}
			} else {
				session = SessionFactory.createSession(this, srvr);
				session.connect();
				session.disconnect();
			}
			
			JOptionPane.showMessageDialog(this, SSH_TUNNEL_SUCCESS_MSG, 
					"SSH-Tunneling Successful", JOptionPane.INFORMATION_MESSAGE);
		} catch (IOException e) {
			ProgrammerLog.log(e);
			JPanel msg = Utilities.collapsedMessage(SSH_TUNNEL_ERROR_MSG, Utilities.toString(e));
			JOptionPane.showMessageDialog(this, msg, 
					"Invalid SSH-Tunnel Configuration", JOptionPane.ERROR_MESSAGE);			
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
		summary.append("<b>Use of SSH-tunnel for MySQL : <i>" + 
				(enableSSHTunnel.isSelected() ? "Enabled" : "Disabled") + "</i></b><br>");
		if (enableSSHTunnel.isSelected()) {
			summary.append(serverPanel.getSummary());
		} else {
			summary.append(Utilities.HTML_4SPACES);
			summary.append("Additional details are not applicable<br>");
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
		final boolean enabled = enableSSHTunnel.isSelected();
		props.setProperty(PropertyKeys.MYSQL_TUNNEL_ENABLED, enabled);
		serverPanel.commitProperties(enabled, PropertyKeys.MYSQL_TUNNEL_SERVER_TYPE,
					PropertyKeys.MYSQL_TUNNEL_SERVER_KNOWN_HOST_KEY,
					PropertyKeys.MYSQL_TUNNEL_SERVER_NAME, PropertyKeys.MYSQL_TUNNEL_SERVER_PORT,
					PropertyKeys.MYSQL_TUNNEL_SERVER_USERID);
	}
	
	/**
	 * The top-level check-box that enables or disables the items in
	 * this wizard-page. If the check-box is checked, then other GUI
	 * items are enabled and permit the user to configure path to shared
	 * storage location. Otherwise the use of shared location is 
	 * disable.
	 */
	private JCheckBox enableSSHTunnel;
	
	/**
	 * A custom component in this package that permits the user to conveniently
	 * select the server where the shared storage is located. The custom
	 * component provides three different mechanisms to specify the location
	 * of the shared storage space.
	 */
	private final ServerInfoPanel serverPanel;
	
	/**
	 * Reference to the wizard that actually owns this page. This is used
	 * when creating sub-dialogs etc. in this page.
	 */
	private final DecagonConfigurationManager dcm;

	/**
	 * A simple message that is displayed to the user at the top of this
	 * wizard page to provide abbreviated information.
	 */
	private static final String SSH_TUNNEL_INFO_MSG = 
		"<html><font size=\"-2\">" +
		"Provide information about the host to be used to<br>" +
		"establish SSH-tunnel to MySQL database. The actual<br>" + 
		"database information is setup on the next page.<br>" +
		"<b>Note: SSH-tunnel is optional.</b>" +
		"</font></html>";
	
	/**
	 * A simple message that is displayed to the user if the SSH-Tunnel
	 * configuration failed to establish connection with the specified
	 * server.
	 */
	private static final String SSH_TUNNEL_ERROR_MSG = 
		"<html>" +
		"The specified SSH-tunnel configuration is invalid.<br>" +
		"Please verify that the following requirements are met:<ul>" +
		"<li>Ensure that server information (hostname/IP, port) etc. are valid</li>"+
		"<li>Ensure you have typed in the correct password.</li>" +
		"</ul>Please appropriately edit the necessary information."+
		"</html>";

	/**
	 * A simple message that is displayed to the user if the SSH-Tunnel
	 * configuration worked as expected.
	 */
	private static final String SSH_TUNNEL_SUCCESS_MSG = 
		"<html>" +
		"The specified SSH-tunnel configuration is good.<br>" +
		"In the next page the actual MySQL database connection<br>" + 
		"will be gathered and database connections will be <br>" +
		"established via the SSH-tunnel." +
		"</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8629449833492205394L;
}
