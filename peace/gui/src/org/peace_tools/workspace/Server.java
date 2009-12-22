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

package org.peace_tools.workspace;

import java.io.PrintWriter;
import java.util.Calendar;

import javax.xml.datatype.DatatypeConfigurationException;
import javax.xml.datatype.DatatypeFactory;
import javax.xml.datatype.Duration;

import org.w3c.dom.Element;

/**
 * Class to encapsulate information regarding a Server entry in a PEACE
 * work space configuration file. A Server entry represents either a single
 * stand alone machine or the head node of a supercomputing cluster on which Jobs
 * can be run. A server entry encapsulates all the information needed to access
 * the server and run jobs on it. In addition, it also provides the necessary
 * infrastructure for marshaling and un-marshaling data for persisting the
 * information in the work space configuration file.
 */
public class Server {
	/**
	 * Different enumerations defining the last known operational status of a
	 * given Server entry. These enumerations were introduced to reflect those
	 * used in the XML and to ensure that the code is overall more readable.
	 */
	public enum ServerStatusType {
		/**
		 * This status indicates that the GUI is making attempt to install
		 * the runtime subsystem of PEACE onto the server. This is a 
		 * relatively long-running process and the server can be in this
		 * state for 5-10 minutes initially. After that a server never
		 * transitions to this state.
		 */
		INSTALLING,
		
		/**
		 * This status indicates that the installation process of the server
		 * failed and the server is unusable.
		 */
		INSTALL_FAILED,
		/**
		 * The server is good to go and is ready for further use.
		 */
		GOOD,
		/**
		 * This status indciates that the GUI is in the process of 
		 * uninstalling the runtime subsystem of PEACE from the remote
		 * machine. In this state, the server is not usable.
		 */
		UNINSTALLING,
		/**
		 * This state indicates that the uninstall attempt on the server
		 * has failed. The server is now in an undefined state and the
		 * user must try to clean up the entry manually.
		 */
		UNINSTALL_FAILED,
		/**
		 * The previous attempt to talk to the server failed and the user
		 * needs to diagnose this server.
		 */
		CONNECT_FAILED
	};
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * Server entry. This method is typically  used to create a suitable
	 * Server entry (either LocalServer or RemoteServer) using data from
	 * a given DOM element.
	 * 
	 * @param serverNode The DOM element to be used for creating the server
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created server entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static Server create(Element serverNode) throws Exception {
		// First obtain the type of this server entry.
		String type = serverNode.getAttribute("type");
		final boolean remote = "remote".equals(type);
		// First extract the necessary information from the DOM tree.
		String ID     = DOMHelper.getStringValue(serverNode, "ID");
		String name   = DOMHelper.getStringValue(serverNode, "Name");
		String desc   = DOMHelper.getStringValue(serverNode, "Description", false);
		desc          = (desc != null) ? desc : "";
		String status = DOMHelper.getStringValue(serverNode, "Status");
		String userId = null;
		if (remote) {
			// Only remote servers have userID
			userId = DOMHelper.getStringValue(serverNode, "UserID");
			
		}
		// Determine the install path.
		String path   = DOMHelper.getStringValue(serverNode, "InstallPath");
		Duration pollTime = null;
		// Obtain polling time as string first, if it exists, and convert it
		if (DOMHelper.hasElement(serverNode, "PollTime")) {
			String pollStr= DOMHelper.getStringValue(serverNode, "PollTime");
			// Now convert polling time from string to a duration.
			DatatypeFactory codec = DatatypeFactory.newInstance();
			pollTime = codec.newDuration(pollStr);
		}
		// Finally create the server object to encapsulate the data.
		Server srvr = new Server(ID, name, desc, userId, path, 
				pollTime, remote);
		srvr.setStatus(ServerStatusType.valueOf(status.toUpperCase()));
		return srvr;
	}
	
	/** 
	 * The constructor. The constructor merely initializes all the instance
	 * variables using the supplied parameters. The parameter values are
	 * either read from a configuration file (via the createServer method)
	 * or obtain from the user (via a suitable GUI dialog).
	 * 
	 * @param ID A unique identifier for this Server entry. For new Server
	 * entries  this value is obtained via the ServerList.reserveServerID() method.
	 * 
	 * @param name The domain name (or IP address) to be used for accessing 
	 * this server. For local machine, this value is simply set to null.
	 * 
	 * @param description A user-assigned description for this server entry.
	 * The description can be anything the user chooses to assign.
	 * 
	 * @param userID The login user ID to be used for accessing remote 
	 * clusters. For the local machine, this value is set to null.
	 * 
	 * @param installPath The location on the Server where PEACE is installed
	 * and the necessary runtime components of PEACE are located.
	 * 
	 * @param pollTime The delay between successive checks for job status
	 * on the server.
	 * 
	 * @param remote This flag indicates if this server entry represents
	 * a local server or a remote server.
	 */
	public Server(String ID, String name, String description,
			String userID, String installPath, Duration pollTime,
			boolean remote) {
		this.ID          = ID;
		this.name        = name.trim();
		this.description = description.trim();
		if (userID != null) {
			this.userID  = userID.trim();
		}
		this.installPath = installPath.trim();
		this.pollTime    = pollTime;
		this.password    = null;
		this.remote      = remote;
		this.status      = ServerStatusType.CONNECT_FAILED; 
	}
	
	/** 
	 * Returns the workspace-unique ID assigned for this server entry.
	 * 
	 * @return Return the workspace-unique ID assigned for this server entry.
	 * This value is used to make cross references to this server entry in
	 * Job and DataFile entries.
	 */
	public String getID() { return ID; }
	
	/**
	 * Returns the server's domain name (or IP address) set for for
	 * this server entry.
	 *  
	 * @return The server's domain name (or IP address).
	 */
	public String getName() { return name; }

	/**
	 * Change the server's domain name (or IP address). Changing the name
	 * is meaningful only for remote entries. For local server (the 
	 * same machine), the name is always null. 
	 * 
	 * @note Changing the server name does not impact any connections that
	 * may be currently open for this server. Note that this method is 
	 * overriden in the LocalServer class to ignore server name changes.
	 * 
	 * @param name The new domain name (or IP address) to be set for
	 * this server entry.
	 */
	public void setName(String name) { this.name = name; }

	/**
	 * Obtain the user-specified description for this entry.
	 * 
	 * @return This method returns the user-specified description for
	 * this entry. This value is always non-null (but could be an 
	 * empty string).
	 */
	public String getDescription() { return description; }
	
	/**
	 * Change the description set for this server entry.
	 * 
	 * @param description The new description to be set for this entry.
	 * If the new description is null, then the description is reset to
	 * an empty string.
	 */
	public void setDescription(String description) {
		this.description = description;
		if (this.description == null) {
			this.description = "";
		}
	}
	
	/**
	 * Obtain the user ID to be used for logging onto this server.
	 * 
	 * @return The userID to be used for logging onto this server. For
	 * local server, the userID is null.
	 */
	public String getUserID() { return userID; }
	
	/**
	 * Set the user ID to be used for this server entry.  
	 * 
	 * @note Calling this method on a local server entry has no effect.
	 * 
	 * @param userID The new user ID to be used for logging onto the
	 * remote server.
	 */
	public void setUserID(String userID) {
		this.userID = userID;
	}
	
	/**
	 * Obtain the directory where the runtime files associated with
	 * PEACE are installed on the server.
	 * 
	 * @return An absolute path indicating the directory where PEACE
	 * runtime files are installed.
	 */
	public String getInstallPath() {
		return installPath;
	}
	
	/**
	 * Set the directory where the runtime files associated with
	 * PEACE are installed on the server.
	 * 
	 * @param path The new install path to be set for this entry.
	 */
	public void setInstallPath(String path) {
		this.installPath = path;
		if (this.description == null) {
			this.description = "";
		}
	}
	
	/**
	 * Obtain the delay between successive polling efforts for this server.
	 * 
	 * @return The delay in seconds between two successive polling/
	 * status checks on this server.
	 */
	public long getPollTime() {
		return pollTime.getTimeInMillis(Calendar.getInstance()) / 1000L;
	}
	
	/**
	 * Change the delay between successive polling efforts on status checks
	 * on this server.
	 * 
	 * @param seconds The delay in seconds to be set between successive
	 * polling efforts (for job status etc.) on this server.
	 */
	public void setPollTime(long seconds) {
		try {
			DatatypeFactory codec = DatatypeFactory.newInstance();
			pollTime = codec.newDuration(seconds * 1000L);
		} catch (DatatypeConfigurationException e) {
			// Cut log entry in programmer log.
			e.printStackTrace();
		}
	}
	
	/**
	 * Set the password to be used for logging on to a remote server.
	 * This method does not have any effect for a local server.
	 * 
	 * @param password The password to be used for logging onto a 
	 * remote server. This value is never persisted but is retained
	 * in memory until this workspace is active.
	 */
	public void setPassword(String password) {
		this.password = password;
	}
	
	/**
	 * Return the password to be used for logging on to a remote 
	 * server. If a password is not set then this method returns
	 * null. An empty password is indicated by the empty string ("").
	 * The password is meaningful only for a remote server.
	 * 
	 * @return The password to be used for logging onto a 
	 * remote server. This value is never persisted but is retained
	 * in memory until this work space is active.
	 */
	public String getPassword() {
		return password;
	}
	
	/**
	 * Method to determine if this Server entry represents a local or a remote
	 * server. 
	 * 
	 * @return This method returns true if this Server object is a remote server.
	 * For a local server entry, this method returns false.
	 */
	public boolean isRemote() { return remote; }
	
	/**
	 * Method to set if this server entry represents a local or a 
	 * remote server.
	 * 
	 * @param remote The parameter must be true to set it as a remote
	 * server. Otherwise the entry is set as a local server.
	 */
	public void setRemote(boolean remote) {
		this.remote = remote;
	}
	
	/**
	 * Change the status for this server.
	 * 
	 * This method also notifies all work space listeners about the
	 * change in the status of this server entry.
	 * 
	 * @param status The new status value to be set for this server entry.
	 */
	public void setStatus(ServerStatusType status) {
		this.status = status;
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(this, WorkspaceEvent.Operation.UPDATE);
		Workspace.get().fireWorkspaceChanged(we);
	}
	
	/**
	 * Obtain the current status set for this server. The job status value
	 * is the last known status value for this server.
	 * 
	 * @return The current status for this server.
	 */
	public ServerStatusType getStatus() { return status; }
	
	/**
	 * Set the unique ID for this server. 
	 * 
	 * This method is typically called only once for each new server
	 * entry from the ServerList.add() method.
	 * 
	 * @param serverID The unique ID value to be set for the server.
	 */
	protected void setID(String serverID) {
		this.ID = serverID;
	}
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent ServerList node in the DOM tree.
	 * 
	 * @param serverList The DOM element corresponding to the ServerList
	 * node that contains this entry.
	 */
	public final void marshall(Element serverList) {
		// Create a top-level server entry for this server
		Element server = DOMHelper.addElement(serverList, "Server", null);
		// Add the type attributes for this server.
		server.setAttribute("type", isRemote() ? "remote" : "local");
		// Add new sub-elements for each value.
		DOMHelper.addElement(server, "ID", ID);
		DOMHelper.addElement(server, "Name", name);
		DOMHelper.addElement(server, "Description", description);
		if (userID != null) {
			DOMHelper.addElement(server, "UserID", userID);
		}
		DOMHelper.addElement(server, "InstallPath", installPath);
		if (pollTime != null) {
			DOMHelper.addElement(server, "PollTime", pollTime.toString());
		}
		DOMHelper.addElement(server, "Status", status.toString().toLowerCase());
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t";
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		// Create a top-level server entry for this server
		out.printf("%s<Server type=\"%s\">\n", Indent, 
				(isRemote() ? "remote" : "local"));
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "ID", ID);
		out.printf(STR_ELEMENT, "Name", name);
		out.printf(STR_ELEMENT, "Description", DOMHelper.xmlEncode(description));
		// UserID is optional. So check and write only if not null
		if (userID != null) {
			out.printf(STR_ELEMENT, "UserID", userID);
		}
		out.printf(STR_ELEMENT, "InstallPath", installPath);
		if (pollTime != null) {
			out.printf(STR_ELEMENT, "PollTime", pollTime.toString());
		}
		// Include status information.
		out.printf(STR_ELEMENT, "Status", status.toString().toLowerCase());
		// Close the server tag
		out.printf("%s</Server>\n", Indent);
	}
	
	@Override
	public String toString() {
		return name;
	}
	
	/**
	 * A unique identifier for this Server entry. For new Server entries 
	 * this value is obtained via the ServerList.reserveServerID() method.
	 */
	private String ID;
	
	/** The domain name (such as: redhawk.hpc.muohio.edu) or IP address
	 * (such as: 134.53.13.131) to be used for accessing the server. 
	 * For local machine, this value is simply set to null.
	 */
	private String name;
	
	/**
	 * A user-assigned description for this server entry.
	 * The description can be anything the user chooses. This is meant
	 * to be meaningful only to the user.
	 */
	private String description;
	
	/**
	 * The login user ID to be used for accessing remote 
	 * clusters. For the local machine, this value is set to null.
	 */
	private String userID;
	
	/**
	 * The location on the Server where PEACE is installed
	 * and the necessary runtime components of PEACE are located.
	 */
	private String installPath;
	
	/**
	 * The delay between successive checks for job status
	 * on the server.
	 */
	private Duration pollTime;
	
	/**
	 * This instance variable indicates if this server entry 
	 * represents a local or a remote server.
	 */
	private boolean remote;
	
	/**
	 * The current operational status of this server. This value
	 * changes periodically as the server is installed, used, and
	 * uninstalled.
	 */
	private ServerStatusType status;
	
	/**
	 * This is a transient field that is never persisted (for security
	 * purposes). Typically, it is set once (each time the GUI is run)
	 * when an attempt is made to access the server.
	 */
	private transient String password;	
}
