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


package org.peace_tools.decagon;

import java.awt.Component;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.xml.datatype.Duration;

import org.peace_tools.core.PEACEProperties;
import org.peace_tools.core.session.RemoteServerSessionMaker;
import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.decagon.PropertyKeys.ServerEntryType;
import org.peace_tools.generic.Log.LogLevel;
import org.peace_tools.generic.PasswordDialog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * Convenience class to encapsulate connection information to MySQL database
 * and facilitate connecting to the database.
 * 
 * This class serves as a central hub to connect to a MySQL database used
 * by DECAGON to store results from running a given genome-assembler for
 * further analysis.
 * 
 */
public class MySQLConnector {

	/**
	 * Obtain the process-wide unique singleton instance.
	 * 
	 * <p>This method must be used to obtain the process-wide singleton instance
	 * of this class to perform various operations supported by it.</p>
	 * 
	 * <p><b>NOTE: This instance should be used for performing all
	 * routine operations.</b></p>
	 * 
	 * @return The process-wide unique singleton instance of this class.
	 * This method never returns a null object.
	 */
	public static MySQLConnector get() {
		return uniqueInstance;
	}

	/**
	 * Constructor for creating temporary connector for testing MySQL 
	 * configuration.
	 * 
	 * This constructor is a convenience constructor that is used by the
	 * MySQLConfigPage to validate MySQL configuration. 
	 *  
	 * @param tunnelServer The MySQL tunnel configuration to be used for
	 * testing. This parameter can be null if tunneling is not to be used. 
	 * 
	 * @param dbServer The host/server information on which MySQL database
	 * is running. This value cannot be null.
	 * 
	 * @param dbName The actual MySQL database to which this connector is
	 * associated. This parameter cannot be null but can be an empty string
	 * if just connection testing is to be performed.
	 */
	public MySQLConnector(final Server tunnelServer, final Server dbServer, 
			final String dbName) {
		this.sshTunnelServer = tunnelServer;
		this.mySqlServer     = dbServer;
		this.dbName          = dbName;
		this.sqlConnection   = null;
	}

	/**
	 * Constructor for creating temporary connector to an alternative database.
	 * 
	 * This constructor is a convenience constructor that is used by 
	 * MySQLDBConfigPage to create temporary connections to a newly created
	 * database.
	 * 
	 * @param src The source connector to be re-purposed to connect to a
	 * different database.
	 * 
	 * @param dbName The name of the database to which a new connection
	 * is to be established. This parameter cannot be null but can be an 
	 * empty string if just connection testing is to be performed.
	 */
	public MySQLConnector(MySQLConnector src, final String dbName) {
		this.sshTunnelServer = src.sshTunnelServer;
		this.mySqlServer     = src.mySqlServer;
		this.dbName          = dbName;
		this.sqlConnection   = null;
	}

	/**
	 * Method to establish a connection to the MySQL database (if needed).
	 * 
	 * This method must be used to establish connection to the MySQL 
	 * database. This method creates necessary tunnel session if so
	 * configured.
	 * 
	 * @param caller The GUI component that is actually calling this method.
	 * This component is used to prompt for password if needed. This
	 * parameter can be null.
	 * 
	 * @param purpose A brief message describing the purpose for the
	 * connection. This is optional and can be null.
	 * 
	 * @return This method returns true if the connection to the database
	 * was successfully established. It returns false if the user cancels
	 * the operation.
	 * 
	 * @throws Exception This method throws various exceptions when errors
	 * occur when attempting to connect to the database.
	 */
	public boolean connect(final Component caller, 
			final String purpose) throws Exception {
		if ((sqlConnection != null) && (sqlConnection.isValid(1))) {
			// Already connected and have valid connection.
			return true;
		}
		// Clear out current connection
		sqlConnection = null;
		// The URL is completed in the if-else block below.
		String dbURL = "jdbc:mysql://"; 
		// First establish SSH-tunnel if needed.
		if (sshTunnelServer != null) {
			tunnelSession = SessionFactory.createSession(caller, sshTunnelServer);
			RemoteServerSessionMaker rssm = 
				new RemoteServerSessionMaker(caller, sshTunnelServer, 
				"Creating SSH-tunnel to access MySQL database");
			rssm.start();
			tunnelSession = rssm.getSession();
			if (tunnelSession == null) {
				// The user must have canceled the SSH-tunnel. Can't proceed further
				return false;
			}
			// Enable port forwarding.
			localPort = tunnelSession.forwardPort(-1, mySqlServer.getName(), mySqlServer.getPort());
			// Cut some logs for troubleshooting.
			UserLog.log(LogLevel.INFO, "DECAGON", "Established port forwarding from " + 
					"localhost:" + localPort + " to " + mySqlServer.getName() + 
					":" + mySqlServer.getPort());
			// Add the necessary information.
			dbURL += "localhost:" + localPort; 
		} else {
			// No port forwarding using the mySqlServer information directly to
			// build the JDBC URL
			dbURL += mySqlServer.getName() + ":" + mySqlServer.getPort();
		}
		// Now complete the URL by adding the optional database name.
		dbURL += "/" + dbName;
		// Load the driver. The following call should never fail 
		// if the mysql-connector-java-5.1.18-bin.jar is in class path.
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		// Try to connect to the MySQL database (no more than 3 retries)
		for(int retries = 3; (retries > 0); retries--) {
			try {
				if (mySqlServer.getPassword() == null) {
					// Prompt the user to enter password.
					PasswordDialog pwdDialog = new PasswordDialog(mySqlServer.getUserID(), false, 
							mySqlServer.getName(), purpose);
					if (!pwdDialog.prompt(caller)) {
						// The user canceled the login
						return false;
					}
					// Save the password supplied by the user.
					mySqlServer.setPassword(pwdDialog.getPassword());
				}
				// Now establish the connection.
				final String password = (mySqlServer.getPassword().length() == 0 ? null : 
					mySqlServer.getPassword());
				// The following call can throw a SQL exception.
				sqlConnection = DriverManager.getConnection (dbURL, mySqlServer.getUserID(), password); 
				UserLog.log(LogLevel.INFO, "DECAGON", "Connection to MySQL database successfully established");
				return true;
			} catch (SQLException exp) {
				sqlConnection = null;
				final String expMsg = exp.getMessage();
				final String errMsg = String.format(SQL_CONNECT_ERROR, 
						Utilities.wrapStringToHTML(expMsg, 75)); 
				JPanel msgPanel = Utilities.collapsedMessage(errMsg, Utilities.toString(exp));
				JOptionPane.showMessageDialog(caller, msgPanel, 
						"Error connecting to Database", JOptionPane.ERROR_MESSAGE);
				// Check to see if it is worth while continuing.
				if ((expMsg == null) || (expMsg.indexOf("Access denied") == -1)) {
					// This is an error message that is not due to invalid password.
					// So bail out.
					return false;
				} else {
					// Retry prompting user for password.
					mySqlServer.setPassword(null);
				}
			}
		}
		// When control drops here that means connection to MySQL was 
		// not established even after 3 tries
		return false;
	}

	/**
	 * Loads the current settings and configuration to access
	 * a MySQL database (if necessary).
	 * 
	 * @param caller The GUI component that is triggering this server
	 * entry creation. This object is used to display error/warning
	 * messages if needed. This parameter can be null.
	 */
	public boolean loadCurrentSettings(Component caller) {
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		final boolean tunnelEnabled = props.getBooleanProperty(PropertyKeys.MYSQL_TUNNEL_ENABLED, false);
		final boolean mySqlEnabled  = props.getBooleanProperty(PropertyKeys.MYSQL_ENABLED, false);

		if ((tunnelEnabled) && (sshTunnelServer == null)) {
			// Create server entry from the global properties setting.
			sshTunnelServer = getServer(PropertyKeys.MYSQL_TUNNEL_SERVER_TYPE, 
					PropertyKeys.MYSQL_TUNNEL_SERVER_KNOWN_HOST_KEY, 
					PropertyKeys.MYSQL_TUNNEL_SERVER_NAME, 
					PropertyKeys.MYSQL_TUNNEL_SERVER_PORT,
					PropertyKeys.MYSQL_TUNNEL_SERVER_USERID, caller);
		}

		if ((mySqlEnabled) && (mySqlServer == null)) {
			// Create server entry from the global properties setting.
			sshTunnelServer = getServer(PropertyKeys.MYSQL_TUNNEL_SERVER_TYPE, 
					"", PropertyKeys.MYSQL_SERVER_NAME, PropertyKeys.MYSQL_SERVER_PORT,
					PropertyKeys.MYSQL_SERVER_USERID, caller);
		}

		return true;
	}

	/**
	 * Helper method to ensure a string is not empty.
	 * 
	 * This is a convenience method that checks to ensure a given string
	 * object is not null and is not an empty string.
	 * 
	 * @param str The string to be validated.
	 * 
	 * @return This method returns true if the parameter is a non-null, 
	 * string with at least one valid (not a white space) character.
	 */
	private boolean isNonEmpty(final String str) {
		return ((str != null) && (str.trim().length() > 0));
	}

	/**
	 * Helper method to create a Server entry using values from
	 * global properties.
	 * 
	 * This method is invoked from {@link #loadCurrentSettings()}
	 * method to create server entries for {@link #mySqlServer} or
	 * {@link #sshTunnelServer}. The required data is loaded from 
	 * the global (cross-workspace) properties. The various keys 
	 * (string defined in PropertyKeys class) to access the values
	 * in PEACEProperties is passed-in as the parameter to this 
	 * method. Depending on the class that is used this object, 
	 * these parameters will vary. However, the set of parameters 
	 * passed to this method must be consistent in order for the 
	 * data to be meaningful.
	 *
	 * @param serverTypePropKey The property key to be used to obtain
	 * the server type from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#MYSQL_TUNNEL_SERVER_TYPE} or
	 * <code>null</code> (for MySQL server).
	 * 
	 * @param serverKnownHostPropKey The property key to be used to obtain
	 * the known host ID value from the global PEACEProperties. This
	 * parameter has values similar to: {@value PropertyKeys#STORAGE_SERVER_KNOWN_HOST_KEY}
	 * or {@link PropertyKeys#MYSQL_SERVER_KNOWN_HOST_KEY}.
	 * 
	 * @param serverNamePropKey The property key to be used to obtain
	 * the server name from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#MYSQL_SERVER_NAME} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_NAME} etc.
	 * 
	 * @param serverPortPropKey The property key to be used to obtain
	 * the SSH-port number from the global PEACEProperties. This value
	 * is similar to: {@link PropertyKeys#MYSQL_SERVER_PORT} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_PORT} etc.
	 * 
	 * @param userIDPropKey The property key to be used to obtain
	 * the user ID (to be used for logging in via SSH) from the global 
	 * PEACEProperties. This value is similar to: 
	 * {@link PropertyKeys#MYSQL_SERVER_USERID} or
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_USERID} etc.
	 * 
	 * @param caller The GUI component that is triggering this server
	 * entry creation. This parameter can be null.
	 */
	public Server getServer(final String serverTypePropKey,
			final String serverKnownHostPropKey, 
			final String serverNamePropKey, final String serverPortPropKey,
			final String userIDPropKey, final Component caller) {
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		// Get the type of server we are dealing with and convert it to an
		// enumeration for further operations (default is local_host)
		final String srvrTypeStr       = props.getProperty(serverTypePropKey);
		final ServerEntryType srvrType = (srvrTypeStr != null ? 
				Enum.valueOf(ServerEntryType.class, srvrTypeStr) : ServerEntryType.LOCAL_HOST);
		Server srvr       = null;
		Duration pollTime = Server.getDefaultPollDuration();

		// Extract additional information to be used in the switch below.
		final String  srvrName = props.getProperty(serverNamePropKey);
		final String  portStr  = props.getProperty(serverPortPropKey);
		final String  userID   = props.getProperty(userIDPropKey);
		final String  serverID = props.getProperty(serverKnownHostPropKey, "DECAGONServer");

		// Create appropriate server entry 
		switch (srvrType) {
		case LOCAL_HOST:
			// Local host case.
			srvr = new Server("SharedServer", "localhost",
					"Shared server for DECAGON", "", "", pollTime, false, 22);
			break;
		case KNOWN_HOST:
			if (((srvr = Workspace.get().getServerList().getServer(serverID)) == null) || 
					(!srvr.getName().equals(srvrName) || !srvr.getUserID().equals(userID) || 
							!("" + srvr.getPort()).equals(portStr))) {
				// The supplied server entry is not valid. Report error to user and
				// check to see if the user would like to proceed further.
				final String msg = String.format(INVALID_KNOWN_HOST_ENTRY, serverID,
						srvrName, portStr, userID);
				final int choice = JOptionPane.showConfirmDialog(caller, msg, 
						"Inconsistent Server Entry", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
				if (choice == JOptionPane.NO_OPTION) {
					// The user does not want to proceed with this entry.
					// Skip it.
					break;
				}
			} else {
				// Found a valid server entry. Nothing further to be done.
				break;
			}
		case CUSTOM_HOST:
			// Ensure all the required information for server entry is set.
			if (isNonEmpty(srvrName) && isNonEmpty(portStr) && isNonEmpty(userID)) {
				final Integer portNum  = new Integer(portStr);
				srvr = new Server(serverID, srvrName, "Server entry for DECAGON", 
						"", userID, pollTime, true, portNum);
			} else {
				// Some information for custom host is missing or is not valid.
				// Show error to the user.
				final String msg = String.format(INVALID_CUSTOM_HOST_ENTRY, srvrName, portStr, userID);
				JOptionPane.showMessageDialog(caller, msg, "Invalid custom host", JOptionPane.ERROR_MESSAGE);		
			}
		}
		return srvr;
	}

	/**
	 * Obtain information about the current MySQL server.
	 * 
	 * This method must be used to obtain the The server entry 
	 * that encapsulates information about the machine that houses the
	 * actual MySQL database. This information is used to connect to 
	 * MySQL (without creating a server session).
	 * 
	 * @return A server object that encapsulates the necessary information
	 * to connect the MySQL database.
	 */
	public Server getMySQLServer() {
		return mySqlServer;
	}

	/**
	 * Obtain the information about the server to be used as SSH-tunnel to
	 * access the MySQL database.
	 * 
	 * @return This method returns null if a SSH-tunnel is not to be used.
	 * Otherwise it returns a Server object that contains the necessary
	 * information to establish a SSH-tunnel.
	 */
	public Server getSSHTunnelServer() {
		return sshTunnelServer;
	}

	/**
	 * Obtain the existing connection to MySQL (if any).
	 * 
	 * This method returns a valid connection to MySQL only if the
	 * {@link #connect(Component, String)} method has been called
	 * and returns true. 
	 * 
	 * @return The existing connection (if any) to MySQL. If a connection
	 * has not been established, then this method returns null.
	 */
	public Connection getConnection() {
		return sqlConnection;
	}

	/**
	 * Obtain name of the database set for this connector.
	 * 
	 * @return The database set for this connector. If the connector is
	 * not set to a specific database then this method returns an
	 * empty string.
	 */
	public String getDBName() {
		return dbName;
	}

	/**
	 * The default and only constructor for this class.
	 * 
	 * This class is implemented as a singleton pattern. To ensure that
	 * this class is instantiated only once, the constructor has been
	 * made private. 
	 * 
	 * @see #get()
	 */
	private MySQLConnector() {
		mySqlServer     = null;
		sshTunnelServer = null;
		dbName          = null;
		tunnelSession   = null;
		localPort       = -1;
	}

	/**
	 * The server entry that encapsulates information about the machine
	 * that houses the actual MySQL database. This information is used
	 * to connect to MySQL (without creating a server session). Even 
	 * though this object is not strictly used in the sense of a regular
	 * server, it is re-purposed here to encapsulate the necessary
	 * information.
	 */
	private Server mySqlServer;

	/**
	 * Information about the server to be used as an SSH-tunnel to access
	 * the MySQL database. This server entry is optional and is used
	 * only when a SSH-tunnel has been enabled. If an SSH-tunnel is not
	 * being used then this object is null.
	 */
	private Server sshTunnelServer;

	/**
	 * An error message that is formatted and displayed to the user
	 * whenever a known-host entry is not found in the current
	 * workspace. 
	 */
	private static final String INVALID_KNOWN_HOST_ENTRY = "<html>" +
	"The configured server (with ID: %s) was not<br>" +
	"found in the workspace or is inconsistent<br>" + 
	"This entry has been configured to refer to:<br>" +
	"&nbsp; &nbsp; &nbsp; Host: %s (port: %s)<br>" +
	"&nbsp; &nbsp; &nbsp; User ID: %s<br>" +
	"It would be best to reconfigure this server entry.<br><br>" +
	"<b>Nevertheless, would you like to proceed with<br>" +
	"with the the above information for now?</b>" + 
	"</html>";

	/**
	 * An error message that is formatted and displayed to the user
	 * whenever a custom-host entry is not found to be consistent.
	 */
	private static final String INVALID_CUSTOM_HOST_ENTRY = "<html>" +
	"The custom server entry does not have complete<br>" +
	"information set for it and is not usable.<br>" + 
	"This entry has been configured to refer to:<br>" +
	"&nbsp; &nbsp; &nbsp; Host: %s (port: %s)<br>" +
	"&nbsp; &nbsp; &nbsp; User ID: %s<br>" +
	"You need to reconfigure this server entry." +
	"</html>";

	/**
	 * The SSH-tunnel session that is being used to access the MySQL
	 * database. This entry is initialized to null. A suitable
	 * SSH-tunnel is created by the {@link #connect(Component)} method.
	 */
	private ServerSession tunnelSession;

	/**
	 * This value is set to the local port that was used to forward
	 * database connection to the MySQL server. This value is meaningful
	 * only when tunneling is used and {@link #tunnelSession} is valid.
	 */
	private int localPort;

	/**
	 * The actual name of the database in MySQL to which this connector
	 * is interacting with.
	 */
	private String dbName;

	/**
	 * The actual JDBC SQL connection to the MySQL database. This connection
	 * is initialized to null in the constructors and setup to a valid
	 * connection by the {@link #connect(Component)} method.
	 */
	private Connection sqlConnection;

	/**
	 * Just a simple error message that is displayed when connection
	 * to the MySQL database could not be successfully established.
	 */
	private static final String SQL_CONNECT_ERROR = "<html>" + 
	"Unable to connect to MySQL database.<br>" +
	"Please verify username and password are valid.<br>" +
	"The error reported was:<ul><li><i>%s</i></li></ul></html>";

	/**
	 * The process-wide unique singleton instance of this class.
	 * 
	 * @see #get()
	 */
	private static final MySQLConnector uniqueInstance = new MySQLConnector(); 
}
