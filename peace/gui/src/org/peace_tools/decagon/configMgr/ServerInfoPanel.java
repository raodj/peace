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
import java.awt.CardLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.xml.datatype.Duration;

import org.peace_tools.core.PEACEProperties;
import org.peace_tools.core.job.ServerComboBox;
import org.peace_tools.decagon.PropertyKeys;
import org.peace_tools.decagon.PropertyKeys.ServerEntryType;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

/**
 * A custom component that provides a convenient GUI to let user specify
 * server information.
 * 
 * This component is a custom component used by the DECAGON configuration
 * manager to obtain server information from the user. The class provides
 * the following three different strategies (that can be appropriately 
 * configured) for specifying remote server information:
 * 
 * <ol>
 * 
 * <li>The local machine or <code>localhost</code>. This setup requires
 * no other additional configuration information.</li>
 * 
 * <li>A pre-existing server configuration in the workspace.</li>
 * 
 * <li>A custom server given the host name (or IP address), the SSH port,
 * and user ID.</li>
 *  
 * </ol>
 */
public class ServerInfoPanel extends JPanel implements ActionListener {
	/**
	 * The one and only constructor for this class.
	 * 
	 * The constructor creates and lays out the various GUI sub-components
	 * for this class.
	 * 
	 * @param includeLocalHost If this flag is true then local host is 
	 * included as list of options. Otherwise local host is not included.
	 * 
	 * @param includeCustomHost If this flag is true then custom host is
	 * included as list of options. Otherwise custom host is not included.
	 */
	public ServerInfoPanel(final boolean includeLocalHost, final boolean includeCustomHost) {
		super(new BorderLayout(0, 5));
		// Create the combo-box that provides the list of known/pre-configured servers
		knownServerList = new ServerComboBox();
		// Populate the list with the list of known servers.
		if (ServerComboBox.getServerList(false, false, this.knownServerList) == 0) {
			// No valid server entries were found. This combo-box should be disabled.
			knownServerList.setEnabled(false);
		}
		// Create the combo box option to choose the type of server.
		serverType = createServerTypeSelectionPanel(includeLocalHost,
				(knownServerList.getItemCount() > 0));

		// Create the panel with custom server information entry fields.
		customInfoPanel = createCustomServerInfoPanel();
		// Create our local-host informational label just as a place-holder.
		localHostInfo = new JLabel(LOCAL_HOST_INFO, 
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		localHostInfo.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEtchedBorder(),
				BorderFactory.createEmptyBorder(2, 5, 2, 5)));
		// Create place holder for server type specific configuration where only
		// one option is displayed at a time (due to card layout) and add the
		// three server-type-specific GUI components.
		serverOptionPanel = new JPanel(new CardLayout());
		add(serverOptionPanel, BorderLayout.SOUTH);
		if (includeLocalHost) {
			serverOptionPanel.add(localHostInfo, SERVER_TYPE_LIST[0]);
		}
		if (knownServerList.getItemCount() > 0) {
			serverOptionPanel.add(knownServerList, SERVER_TYPE_LIST[1]);
		}
		if (includeCustomHost) {
			serverOptionPanel.add(customInfoPanel, SERVER_TYPE_LIST[2]);
		}
		// Initial default setting to the first entry in the combo-box.
		setServerType(0); 
	}

	/**
	 * Helper method to create a Server entry based on current user selections.
	 * 
	 * This is a helper method that is used to create a suitable Server entry
	 * based on the current user's selection.
	 *  
	 * @return A server entry based on the current selections made by the user.
	 */
	public Server getServerEntry() {
		Server srvr       = null;
		Duration pollTime = Server.getDefaultPollDuration();

		if (SERVER_TYPE_LIST[0].equals(serverType.getSelectedItem())) {
			// Local host case.
			srvr = new Server("SharedServer", "localhost",
					"Shared server for DECAGON", "", "", pollTime, false, 22);
		} else if (SERVER_TYPE_LIST[1].equals(serverType.getSelectedItem())) {
			// Use an existing server entry.
			srvr = (Server) knownServerList.getSelectedItem();
		} else if (SERVER_TYPE_LIST[2].equals(serverType.getSelectedItem())) {
			final int portNum = (Integer) port.getValue();
			// Check if existing custom server entry is good.
			if ((customServer == null) || 
					(!customServer.getName().equals(customHost.getText())) || 
					(customServer.getPort() != portNum) || 
					(!customServer.getUserID().equals(userName.getText()))) {
				// The existing custom server entry is not good. 
				// Create custom server entry and use it to create session.
				customServer = new Server("SharedServer", customHost.getText(), 
						"Shared server for DECAGON", userName.getText(),
						"", pollTime, true, portNum);
				srvr = customServer;
			}
			// Now customServer is ready for use.
			srvr = customServer;
		}
		return srvr;
	}

	/**
	 * Helper method to appropriately enable/disable GUI components
	 * managed by this component.
	 * 
	 * This is a convenience method that can be used to appropriately enable/
	 * disable various input GUI components managed by this class.
	 *  
	 * @param enableFlag If this flag is true then the components are
	 * enabled permitting the user to edit values.
	 */
	public void enableInputs(final boolean enableFlag) {
		Utilities.setEnabled(this, enableFlag);
		knownServerList.setEnabled(enableFlag && knownServerList.getItemCount() > 0);
	}

	/**
	 * This method can be used to obtain the JPanel containing custom  
	 * server information.
	 * 
	 * <p>This method returns a JPanel that contains the host name, port, 
	 * and user ID information. The GUI components are labeled and
	 * appropriately organized in the JPanel and can be readily used.</p>
	 * 
	 * <p>NOTE: This method sets the chosen input type to be a custom
	 * server configuration.</p>
	 * 
	 * @return The custom panel with GUI components for obtaining
	 * information about host name, port, and user ID.
	 */
	public JPanel getCustomInfoPanel() {
		knownServerList.setSelectedItem(SERVER_TYPE_LIST[2]);
		setServerType(2);
		return customInfoPanel;
	}

	/**
	 * Helper method to create the server selection combo-box.
	 * 
	 * This is a helper method that is called only once from the constructor.
	 * This method was introduced to streamline the code in the constructor.
	 * This method creates the server-type selection combo-box and adds it
	 * to the top of this component.
	 * 
	 * @param includeLocalHost If this flag is true then local host is 
	 * included as list of options. Otherwise local host is not included.
	 * 
	 * @param includeKnownHost If this flag is true then known host option
	 * is included as list of options.
	 * 
	 * @return This method returns the combo-box created by this method.
	 */
	private JComboBox createServerTypeSelectionPanel(final boolean includeLocalHost,
			final boolean includeKnownHost) {
		// Create the combo box option to choose the type of server.
		// The following method call sets
		JComboBox srvrTypeList = new JComboBox();
		if (includeLocalHost) {
			// We make an assumption that local host is entry 0.
			srvrTypeList.addItem(SERVER_TYPE_LIST[0]);
		}
		if (includeKnownHost) {
			// We make an assumption that local host is entry 1.
			srvrTypeList.addItem(SERVER_TYPE_LIST[1]);
		}
		srvrTypeList.addItem(SERVER_TYPE_LIST[2]);
		
		srvrTypeList.setActionCommand("ServerType");
		srvrTypeList.addActionListener(this);
		// Layout the top-level server-type selection component
		JPanel serverTypePanel  = new JPanel();
		serverTypePanel.setLayout(new BoxLayout(serverTypePanel, BoxLayout.X_AXIS));
		serverTypePanel.add(new JLabel("Select server Type: "));
		serverTypePanel.add(srvrTypeList);
		// Add combo-box with title to the top of this component
		add(serverTypePanel, BorderLayout.NORTH);
		// return the newly created component for further use
		return srvrTypeList;
	}

	/**
	 * Helper method to create the GUI to enter information about a custom server.
	 * 
	 * This method is invoked only once from the constructor to create the
	 * GUI components used to enter information about a custom server. This
	 * method creates the host name, port, and userID entries. It wraps the
	 * GUI elements in a suitable panel and returns it back to the caller.
	 *  
	 * @return A panel containing custom server information entry boxes.
	 */
	private JPanel createCustomServerInfoPanel() {
		// Create the various input fields that we are going to use.
		customHost = new JTextField(15);
		port       = new JSpinner(new SpinnerNumberModel(22, 1, 65535, 1));
		userName   = new JTextField(9);
		// Wrap the above components into a panel.
		JPanel customSrvrInfo = new JPanel();
		customSrvrInfo.setLayout(new BoxLayout(customSrvrInfo, BoxLayout.X_AXIS));
		customSrvrInfo.setAlignmentX(0);
		// Add components to the custom server information panel
		customSrvrInfo.add(Utilities.createLabeledComponents("Host name or IP:", null, 2, false, customHost));
		customSrvrInfo.add(Box.createHorizontalStrut(5));
		customSrvrInfo.add(Utilities.createLabeledComponents("SSH Port:", null, 0, false, port));
		customSrvrInfo.add(Box.createHorizontalStrut(5));
		customSrvrInfo.add(Utilities.createLabeledComponents("User ID:", null, 2, false, userName));
		// Return the customized panel back to the caller.
		return customSrvrInfo;
	}

	/**
	 * Helper method to change GUI element settings based on server type.
	 * 
	 * This a utility method that is used internally to enable/disable
	 * various GUI elements based on the current server type. The valid
	 * strings that can be passed this method and the corresponding
	 * action performed are:
	 * 
	 * <ul>
	 * 
	 * <li> When serverType is 0 (localHost) this method disables
	 * known server combo-box and the input fields for custom server.</li>
	 * <li> When serverType is 1 (<code>knownHost</code>) this method enables
	 * known server combo-box but disables the input fields for custom server.</li>
	 * <li> When serverType is 2 (<code>customHost</code>) this method disables
	 * known server combo-box but enables the input fields for custom server.</li>
	 * </ul>
	 * 
	 * @param serverTypeID The type of server setup currently desired.
	 */
	private void setServerType(final int serverTypeID) {
		// Update the server option appropriately.
		CardLayout cl = (CardLayout)(serverOptionPanel.getLayout());
		String  page  = (String) serverType.getItemAt(serverTypeID);
		cl.show(serverOptionPanel, page);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		this.setServerType(serverType.getSelectedIndex());
	}

	/**
	 * Returns summary information about the server configuration
	 * as a HTML fragment. 
	 * 
	 * This method is a convenience method that can be used obtain server
	 * summary information from this component. The information provided here
	 * is sufficiently generic.
	 * 
	 * @return A HTML-sub-fragment providing information about the
	 * configuration to be committed by this class.
	 */
	protected String getSummary() {
		final Server tmpSrvr  = this.getServerEntry();
		return String.format(SERVER_SUMMARY_INFO, serverType.getSelectedItem(),
				tmpSrvr.getID(), tmpSrvr.getName(), tmpSrvr.getPollTime(), 
				tmpSrvr.getDescription(), tmpSrvr.getUserID());
	}

	/**
	 * Convenience method to commit user-entered properties to the global
	 * properties. This method is invoked by the top-level Wizard once
	 * the user has verified the configuration.
	 * 
	 * @param setProperty If this flag is true then the values are set 
	 * in the global PEACEProperties. If this flag is false then the
	 * values are removed from the global PEACEPropreties.
	 * 
	 *  @param serverTypePropKey The property key to be used to obtain
	 * the server type from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#STORAGE_SERVER_TYPE} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_TYPE} etc.
	 * 
	 * @param knownHostPropKey The property key to be used to obtain
	 * the known host ID value from the global PEACEProperties. This
	 * parameter has values similar to: {@value PropertyKeys#STORAGE_SERVER_KNOWN_HOST_KEY}
	 * or {@link PropertyKeys#MYSQL_SERVER_KNOWN_HOST_KEY}.
	 * 
	 * @param serverNamePropKey The property key to be used to obtain
	 * the server name from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#STORAGE_SERVER_NAME} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_NAME} etc.
	 * 
	 * @param serverPortPropKey The property key to be used to obtain
	 * the SSH-port number from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#STORAGE_SERVER_PORT} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_PORT} etc.
	 * 
	 * @param userIDPropKey The property key to be used to obtain
	 * the user ID (to be used for logging in via SSH) from the global 
	 * PEACEProperties. This value is similar to: 
	 * {@link PropertyKeys#STORAGE_SERVER_USERID} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_USERID} etc.
	 */
	protected void commitProperties(final boolean enabled, final String serverTypePropKey, 
			final String knownHostPropKey,
			final String serverNamePropKey, final String serverPortPropKey,
			final String userIDPropKey) {
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		if (enabled) {
			// Get a temporary server entry to ease creation of server
			final Server tmpSrvr = getServerEntry();
			// Infer the server type.
			final String srvrKind = (String) serverType.getSelectedItem();
			ServerEntryType srvrType = ServerEntryType.LOCAL_HOST;
			if (SERVER_TYPE_LIST[1].equals(srvrKind)) {
				srvrType = ServerEntryType.KNOWN_HOST;
			} else if (SERVER_TYPE_LIST[2].equals(srvrKind)) {
				srvrType = ServerEntryType.CUSTOM_HOST;
			}
			props.setProperty(serverTypePropKey, srvrType.toString());
			// Set other properties about the server.
			props.setProperty(knownHostPropKey, tmpSrvr.getID());
			props.setProperty(serverNamePropKey, tmpSrvr.getName());
			props.setProperty(serverPortPropKey, "" + tmpSrvr.getPort());
			props.setProperty(userIDPropKey, tmpSrvr.getUserID());
		} else {
			// Clear out the properties
			props.removeProperty(serverTypePropKey);
			props.removeProperty(knownHostPropKey);
			props.removeProperty(serverNamePropKey);
			props.removeProperty(serverPortPropKey);
			props.removeProperty(userIDPropKey);
		}
	}

	/**
	 * Convenience method to update information displayed in this 
	 * component. 
	 * 
	 * This method is typically used by owning component to update
	 * the data just before a page is to be displayed. The data is
	 * obtained from the global PEACEProperties configuration management
	 * class. The various keys (string defined in PropertyKeys class)
	 * to access the values in PEACEProperties is passed-in as the 
	 * parameter to this method. Depending on the class that is used
	 * this object, these parameters will vary. However, the set of 
	 * parameters passed to this method must be consistent in order for
	 * the data to be meaningful.
	 *
	 * @param serverTypePropKey The property key to be used to obtain
	 * the server type from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#STORAGE_SERVER_TYPE} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_TYPE} etc.
	 * 
	 * @param knownHostPropKey The property key to be used to obtain
	 * the known host ID value from the global PEACEProperties. This
	 * parameter has values similar to: {@value PropertyKeys#STORAGE_SERVER_KNOWN_HOST_KEY}
	 * or {@link PropertyKeys#MYSQL_SERVER_KNOWN_HOST_KEY}.
	 * 
	 * @param serverNamePropKey The property key to be used to obtain
	 * the server name from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#STORAGE_SERVER_NAME} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_NAME} etc.
	 * 
	 * @param serverPortPropKey The property key to be used to obtain
	 * the SSH-port number from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#STORAGE_SERVER_PORT} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_PORT} etc.
	 * 
	 * @param userIDPropKey The property key to be used to obtain
	 * the user ID (to be used for logging in via SSH) from the global 
	 * PEACEProperties. This value is similar to: 
	 * {@link PropertyKeys#STORAGE_SERVER_USERID} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_USERID} etc.
	 */
	public void updateServerInfo(final String serverTypePropKey, 
			final String knownHostPropKey,
			final String serverNamePropKey, final String serverPortPropKey,
			final String userIDPropKey) {
		final int knownSrvrCount = knownServerList.getItemCount();
		knownServerList.setEnabled(knownSrvrCount > 0);
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();

		// Setup current values if they are known and valid.
		final String srvrTypeStr = props.getProperty(serverTypePropKey); 
		final ServerEntryType srvrType = (srvrTypeStr != null ? 
				ServerEntryType.valueOf(ServerEntryType.class, srvrTypeStr) : 
					ServerEntryType.LOCAL_HOST);
		// Load the server ID key from the list.
		final String serverID = props.getProperty(knownHostPropKey);

		// Flag to determine if a custom-server entry was found in the loop below.
		boolean entryFound = false;
		if (ServerEntryType.KNOWN_HOST.equals(srvrType) && (serverID != null)) {
			// Select server entry from the known server list.
			for(int i = 0; (i < knownServerList.getItemCount()); i++) {
				Server srvr = (Server) knownServerList.getItemAt(i);
				if (srvr.getID().equals(serverID)) {
					knownServerList.setSelectedIndex(i);
					entryFound = true;
					break;
				}
			}
			// Check to see if the entry was found. If not report a warning
			if (!entryFound) {
				JOptionPane.showMessageDialog(this, BAD_SERVER_ENTRY,
						"Server Entry Invalid", JOptionPane.WARNING_MESSAGE);
			} else {
				// Found the existing server entry
				serverType.setSelectedItem(SERVER_TYPE_LIST[1]);
			}
		}

		if (ServerEntryType.LOCAL_HOST.equals(srvrType)) {
			serverType.setSelectedItem(SERVER_TYPE_LIST[0]);
		} else if (ServerEntryType.CUSTOM_HOST.equals(srvrType) || (!entryFound)) { 				
			// Custom server. Setup host name, port, and user ID
			// Load the server name, port, and host etc. from global properties
			final String  portStr = props.getProperty(serverPortPropKey, "22");
			customHost.setText(props.getProperty(serverNamePropKey));
			port.setValue(new Integer(portStr));
			userName.setText(props.getProperty(userIDPropKey));
			serverType.setSelectedItem(SERVER_TYPE_LIST[2]);
		}
		// Setup the selection appropriately.
		setServerType(serverType.getSelectedIndex());
	}

	/**
	 * The list of different types of server entries that are supported by
	 * this class. This list is used to setup the various types of servers
	 * in the constructor.
	 */
	private static final String[] SERVER_TYPE_LIST = {"Local machine", 
		"Pre-configured Server", "Custom/Other Server"};

	/**
	 * The combo box that permits user to select the type of server
	 * entry to be used. This combo box has server entries of the 
	 * following type in the given order:
	 * 
	 * <ol>
	 * <li>Local host (optional, may or maynot be present)</li>
	 * <li>Known/pre-configured Server</li>
	 * <li>Custom/Other server</li>
	 * </ol>
	 */
	private final JComboBox serverType;

	/**
	 * Field to read/display the IP or name of the custom host for the
	 * user to enter.
	 */
	private JTextField customHost;

	/**
	 * Field to read/display the port number set for a custom server. 
	 * This entry is enabled only for remote servers. The default
	 * port number is 22 (for ssh).
	 */
	private JSpinner port;

	/**
	 * Field to read/display the user/login ID. This information is
	 * meaningful for custom host options only.
	 */
	private JTextField userName;

	/**
	 * The list of known servers that are pre-configured on this workspace.
	 * This combo-box is enabled and used when the user indicates that the
	 * shared space is present on a known server. This list is populated
	 * in the constructor.
	 */
	private ServerComboBox knownServerList;

	/**
	 * A panel that wraps the custom server information such as: host name,
	 * port, and user ID. This panel is used to ease showing/hiding the
	 * GUI elements based on the server type selected by the user.
	 */
	private final JPanel customInfoPanel;

	/**
	 * A simple label that provides the user with information indicating
	 * that no additional selection is necessary in this situation.
	 */
	private final JLabel localHostInfo;

	/**
	 * This panel serves as a place holder to display different GUI controls
	 * based on server type. This panel is created in the constructor and
	 * depending on the server type it contains one of the following: 
	 * localHostInfo, knownServerList, or customInfoPanel
	 */
	private final JPanel serverOptionPanel;

	/**
	 * A simple message that is displayed to the user if an known server
	 * in the current configuration could not be found.
	 */
	private static final String LOCAL_HOST_INFO = 
		"<html>" +
		"Local host refers to this machine on which PEACE is running.<br>" +
		"No additional server information is needed." +
		"</html>";

	/**
	 * A custom server entry that is reused as much as possible to 
	 * avoid having to create server entries frequently. This entry
	 * is managed by the {@link #getServerEntry()} method.
	 */
	private Server customServer;

	/**
	 * A simple message that is displayed to the user if an known server
	 * in the current configuration could not be found.
	 */
	private static final String BAD_SERVER_ENTRY = 
		"<html>" +
		"The existing shared storage server entry is no longer valid.<br>" +
		"The existing entry will be converted to a custom server entry.<br>" +
		"It is best if you reconfigure the shared storage to be used." +
		"</html>";

	/**
	 * A HTML fragment that is formatted to fill-in the necessary information.
	 */
	private static final String SERVER_SUMMARY_INFO = 
		Utilities.HTML_4SPACES + "Server Type: %s<br>" +
		Utilities.HTML_4SPACES + "Server ID: %s<br>" +
		Utilities.HTML_4SPACES + "Host name: %s<br>" +
		Utilities.HTML_4SPACES + "Port: %d<br>" +
		Utilities.HTML_4SPACES + "Description: %s<br>" +
		Utilities.HTML_4SPACES + "User ID: %s<br>";

	/**
	 * A generated serialization UID value to keep the compiler happy. 
	 */
	private static final long serialVersionUID = -621472316559913205L;
}
