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

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.workspace.JobBase.JobType;
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
		 * This status indicates that the GUI is in the process of 
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
	 * Enumerations to provide a more convenient mechanism for 
	 * referring to the actual type of the server that this
	 * server session has been connected to. 
	 */
	public enum OSType {
		UNIDENTIFIED, LINUX, UNIX, WINDOWS
	}

	/**
	 * Enumerations to provide a more convenient mechanism for 
	 * referring to the type of executable program that this
	 * installed on this server. These programs are installed
	 * as a part of the standard installation-package included
	 * as a part of the PEACE-GUI distribution. 
	 */
	public enum EXEKind {
		PEACE_CLUSTER_MAKER,
		EAST_EXE,
		BATON_EXE,
		WIN_LAUNCHER,
		JOB_RUNNER
	}
	
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
		// Find and save the OS type if available.
		OSType osType = OSType.UNIDENTIFIED;
		if (serverNode.hasAttribute("os")) {
			String osStr = serverNode.getAttribute("os");
			osType = OSType.valueOf(osStr.toUpperCase());
		}
		// First extract the necessary information from the DOM tree.
		String ID     = DOMHelper.getStringValue(serverNode, "ID");
		String name   = DOMHelper.getStringValue(serverNode, "Name");
		int    port   = -1;
		if (remote) {
			// Get the port number for the remote server.
			port = DOMHelper.getIntValue(serverNode, "Port");
		}
		String desc   = DOMHelper.getStringValue(serverNode, "Description", true);
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
				pollTime, remote, port);
		srvr.setStatus(ServerStatusType.valueOf(status.toUpperCase()));
		srvr.setOSType(osType);
		// Setup flag to indicate if server has EAST installed in it.
		if (DOMHelper.hasElement(serverNode, "hasEAST")) {
			String hasEASTStr = DOMHelper.getStringValue(serverNode, "hasEAST");
			boolean hasEAST   = Boolean.parseBoolean(hasEASTStr);
			srvr.setEASTInstalled(hasEAST);
		}
		// Setup flag to indicate if server has DECAF installed in it.
		if (DOMHelper.hasElement(serverNode, "hasDECAF")) {
			String hasDECAFStr = DOMHelper.getStringValue(serverNode, "hasDECAF");
			boolean hasDECAF   = Boolean.parseBoolean(hasDECAFStr);
			srvr.setDECAFInstalled(hasDECAF);
		}
		// Setup flag to indicate if server has AMOS installed in it.
		if (DOMHelper.hasElement(serverNode, "hasAMOSTools")) {
			String hasAMOSFStr = DOMHelper.getStringValue(serverNode, "hasAMOSTools");
			boolean hasAMOS    = Boolean.parseBoolean(hasAMOSFStr);
			srvr.setHasAMOSTools(hasAMOS);
		}
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
	 * 
	 * @param port The port number (meaningful only for remote servers)
	 * over which secure connections are to be established. The default
	 * value is 22.
	 */
	public Server(String ID, String name, String description,
			String userID, String installPath, Duration pollTime,
			boolean remote, int port) {
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
		this.hasEAST     = false;
		this.hasDECAF    = false;
		this.hasAMOStools= false;
		this.osType      = OSType.UNIDENTIFIED;
		this.port        = port;
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
	 * Change the server's domain name (or IP address). Changing the name is
	 * meaningful only for remote entries. For local server (the same machine),
	 * the name is always null.
	 * 
	 * <p>
	 * <b>Note:</b> Changing the server name does not impact any connections
	 * that may be currently open for this server. Note that this method is
	 * overridden in the LocalServer class to ignore server name changes.
	 * </p>
	 * 
	 * @param name The new domain name (or IP address) to be set for 
	 * this server entry.
	 */
	public void setName(String name) { this.name = name; }

	/**
	 * Returns the server's port number set for this server entry.
	 * The value returned by this method is meaningful only for
	 * remote server entries. Local servers have this value set
	 * to -1.
	 * 
	 * @return The server's port number over which connections to this
	 * remote server is to be established. Typically this value is 22
	 * (for ssh)
	 * 
	 * @see #isRemote()
	 */
	public int getPort() { return port; }

	/**
	 * Change a remote server's port number. Changing the port number
	 * meaningful only for remote entries. For local server (the same machine),
	 * the port number is -1 and this value is ignored.
	 * 
	 * <p>
	 * <b>Note:</b> Changing the port number does not impact any connections
	 * that may be currently open for this server. Note that this method is
	 * overridden in the LocalServer class to ignore port number changes.
	 * </p>
	 * 
	 * @param port The new port number to be set for this server entry.
	 */
	public void setPort(int port) { this.port = port; }

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
	 * <p><b>Note:</b>  Calling this method on a local server entry 
	 * has no effect.</p>
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
	 * <p>Method to determine if this server contains EAST (a MST-based
	 * assembly software) installed on this server. This flag is used
	 * to determine the list of servers on which EAST can be run.</p>
	 * 
	 * <b>Note:</b> that this flag is meaningful only if the status 
	 * of the server is GOOD.
	 * 
	 * @return This method returns true if the Server object has
	 * EAST installed on it. Otherwise (and by default) it returns false.
	 */
	public boolean hasEASTInstalled() { return hasEAST; }
	
	/**
	 * Method to set if this server has a valid install of EAST. This
	 * flag is set when a new server entry is added to a Workspace 
	 * to indicate if it has EAST installed. 
	 * 
	 * @param hasEAST If this flag is set to true then it indicates that
	 * this server entry has the EAST assembler installed on it.
	 */
	public void setEASTInstalled(boolean hasEAST) {
		this.hasEAST = hasEAST;
	}

	/**
	 * <p>Method to determine if this server has AMOS tools installed on 
	 * it. This flag is used to determine if the server can convert
	 * ACE file formats generated by PEACE to SAM file format.</p>
	 * 
	 * <b>Note:</b> that this flag is meaningful only if the status 
	 * of the server is GOOD and the server has EAST installed on it.
	 * 
	 * @return This method returns true if the Server object has
	 * AMOS-tools installed on it. Otherwise (and by default) it returns false.
	 */
	public boolean hasAMOSTools() { return hasAMOStools; }
	
	/**
	 * Method to set if this server has a valid install of AMOS-tools. This
	 * flag is set when a new server entry is added to a Workspace 
	 * to indicate if it has AMOS installed. This method is currently
	 * used by the ServerInfoWizardPage 
	 * 
	 * @param hasAMOStools If this flag is set to true then it indicates that
	 * this server entry has AMOS-tools installed on it.
	 */
	public void setHasAMOSTools(boolean hasAMOStools) {
		this.hasAMOStools = hasAMOStools;
	}

	/**
	 * <p>Method to determine if this server contains DECAF (a Distributed
	 * Empirical Comparison and Analysis Framework) installed on this server. 
	 * This flag is used to determine the list of servers on which empirical
	 * comparative analysis can be run.</p>
	 * 
	 * <b>Note:</b> that this flag is meaningful only if the status 
	 * of the server is GOOD.
	 * 
	 * @return This method returns true if the Server object has
	 * DECAF installed on it. Otherwise (and by default) it returns false.
	 */
	public boolean hasDECAFInstalled() { return hasDECAF; }
	
	/**
	 * Method to set if this server has a valid install of DECAF. This
	 * flag is set when a new server entry is added to a Workspace 
	 * to indicate if it has DECAF installed. 
	 * 
	 * @param hasDECAF If this flag is set to true then it indicates that
	 * this server entry has DECAF installed and ready to use.
	 */
	public void setDECAFInstalled(boolean hasDECAF) {
		this.hasDECAF = hasDECAF;
	}

	/** Determine the type of OS on the server.
	 * 
	 * This method can be used to determine the type of operating system
	 * (OS) installed on the server. This information is often used to
	 * appropriately interact with the server to start jobs and perform
	 * other operations. The OS type is typically set once, when a new
	 * server entry is added to the workspace. It is persisted and loaded
	 * from the workspace configuration file.
	 * 
	 * @return The type of the OS installed on the remote server. 
	 */
	public OSType getOSType() {
		return osType;
	}
	
	/** Set the type of OS running on the server.
	 * 
	 * This method is typically used just once to set the type of OS
	 * installed/running on the server. This method is used when a
	 * new Server entry is being added to the workspace.
	 * 
	 * @param osType The type of OS running on the server. This value
	 * is typically determined from {@link ServerSession#getOSType()} 
	 * method.
	 *
	 */
	public void setOSType(OSType osType) {
		this.osType = osType;
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
		server.setAttribute("os", osType.toString().toLowerCase());
		// Add new sub-elements for each value.
		DOMHelper.addElement(server, "ID", ID);
		DOMHelper.addElement(server, "Name", name);
		if (isRemote()) {
			// Add a port entry for remote servers
			DOMHelper.addElement(server, "Port", port);
		}
		DOMHelper.addElement(server, "Description", description);
		if (userID != null) {
			DOMHelper.addElement(server, "UserID", userID);
		}
		DOMHelper.addElement(server, "InstallPath", installPath);
		if (pollTime != null) {
			DOMHelper.addElement(server, "PollTime", pollTime.toString());
		}
		DOMHelper.addElement(server, "Status",  status.toString().toLowerCase());
		DOMHelper.addElement(server, "hasEAST", Boolean.toString(hasEAST));
		DOMHelper.addElement(server, "hasDECAF", Boolean.toString(hasDECAF));
		DOMHelper.addElement(server, "hasAMOSTools", Boolean.toString(hasAMOStools));
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
		final String INT_ELEMENT = Indent + "\t" + "<%1$s>%2$d</%1$s>\n";
		// Create a top-level server entry for this server
		out.printf("%s<Server type=\"%s\" os=\"%s\">\n", Indent, 
				(isRemote() ? "remote" : "local"), 
				osType.toString().toLowerCase());
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "ID", ID);
		out.printf(STR_ELEMENT, "Name", name);
		if (isRemote()) {
			out.printf(INT_ELEMENT, "Port", port);
		}
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
		// Include flag to indicate if EAST assembler is installed on this server
		out.printf(STR_ELEMENT, "hasEAST", Boolean.toString(hasEAST));
		// Include flag to indicate if DECAF evaluation framework is installed on this server
		out.printf(STR_ELEMENT, "hasDECAF",     Boolean.toString(hasDECAF));
		out.printf(STR_ELEMENT, "hasAMOSTools", Boolean.toString(hasAMOStools));
		// Close the server tag
		out.printf("%s</Server>\n", Indent);
	}
	
	@Override
	public String toString() {
		// Display non-default port number for remote server entry.
		final String portStr = (isRemote() && (port != 22) ? ":" + port : "");
		return name + portStr;
	}
	
	/**
	 * Helper method to return a fully qualified path to a given executable
	 * on this server.
	 * 
	 * This is a convenience helper method that can be used to obtain the
	 * fully qualified path to a given executable on this server.
	 * 
	 * @return A fully qualified path to a requested executable.
	 */
	public String getExecutable(EXEKind kind) {
		String peacePath     = (osType.equals(OSType.WINDOWS) ? "/" : "/peace/src/");
		String eastPath      = (osType.equals(OSType.WINDOWS) ? "/" : "/peace/EAST/C++/");
		String exeSuffix     = (osType.equals(OSType.WINDOWS) ? ".exe" : "");
		String scriptSuffix  = (osType.equals(OSType.WINDOWS) ? ".bat" : ".sh");
		String jobRunnerPath = "installFiles/";
		jobRunnerPath       += (osType.equals(OSType.WINDOWS) ? "windows" : "linux");
		
		switch(kind) {
		case PEACE_CLUSTER_MAKER:
			return installPath + peacePath + "peace" + exeSuffix;
		case EAST_EXE:
			return installPath + eastPath + "Main" + exeSuffix;
		case WIN_LAUNCHER:
			return installPath + "\\launcher" + exeSuffix;
		case JOB_RUNNER:
			return installPath + jobRunnerPath + "jobRunner" + scriptSuffix;
		default:
			return null;	
		}
	}
	
	/**
	 * Helper method to return a fully qualified path to a given executable
	 * on this server for a given type of job.
	 * 
	 * This is a convenience helper method that can be used to obtain the
	 * fully qualified path to a given executable on this server.
	 * 
	 * @param type The type of the job for which the executable is to be 
	 * returned.
	 * 
	 * @return A fully qualified path to a requested executable.
	 */
	public String getExecutable(JobType type) {
		String peacePath     = (osType.equals(OSType.WINDOWS) ? "/" : "/peace/src/");
		String eastPath      = (osType.equals(OSType.WINDOWS) ? "/" : "/peace/EAST/C++/");
		String exeSuffix     = (osType.equals(OSType.WINDOWS) ? ".exe" : "");
		
		switch(type) {
		case CLUSTERING:
			return installPath + peacePath + "peace" + exeSuffix;
		case EAST:
			return installPath + eastPath + "east.sh" + exeSuffix;
		default:
			return null;	
		}
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
	 * <p>The port number over which remote servers are to be contacted.
	 * For most traditional SSH installations, the default port is 22.
	 * However, for non-traditional hosts, the port number can vary.
	 * Varying the port number permits creation of tunnels etc. which
	 * makes it convenient to work around fire walls or with multiple
	 * clients.</p>
	 * 
	 * <p>Note that the port number is meaningful only for remote servers
	 * whose {@link #remote} flag is set to true.</p>
	 * 
	 * @see #isRemote()
	 */
	private int port;
	
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
	 * Flag to indicate if this server entry has EAST (the MST-based
	 * assembler) installed. This value is set when a Server entry is
	 * added to an workspace. This value is used to decide if EAST
	 * can be run on a server to perform assembly (after clustering).
	 */
	private boolean hasEAST;
	
	/**
	 * Flag to indicate if this server entry has DECAF (PEACE's
	 * Distributed Empirical Comparison and Analysis Framework)
	 * installed on this server.
	 */
	private boolean hasDECAF;
	
	/**
	 * Flag to indicate if this server has AMOS tools installed on it.
	 * 
	 * <p>AMOS (see http://sourceforge.net/apps/mediawiki/amos/) is a third-party 
	 * tool set that is used to convert ACE output file format generated by 
	 * EAST assembler to SAM file format. Similar to ACE, SAM is another widely 
	 * used text file format for storing assembly
	 * results (see http://samtools.sourceforge.net/). The process of converting
	 * an ACE file format to SAM is documented at: 
	 * http://sourceforge.net/apps/mediawiki/amos/index.php?title=Bank2contig#SAM_Conversion
	 * </p> 
	 * 
	 * <p><b>Note</b>We do not distribute AMOS tools with PEACE and we expect
	 * the user (or their system administrator) to install and setup
	 * AMOS tools for us. However, in this method we simply detect to see if
	 * the following three AMOS-tools that we use are available by running
	 * the following commands to simply return their version information:
	 * 
	 * <ul>
	 * <li>toAmos -v</li>
	 * <li>bank-transact -v</li>
	 * <li>bank2contig -v</li>
	 * </ul>
	 * 
	 * </p> 
	 * 
	 * If the above tools are present when this server entry was created,
	 * this this flag is set to true. The detection of tool is
	 * performed by the ServerInfoWizardPage class.
	 */
	private boolean hasAMOStools;
	
	/**
	 * The operating system type for this server. This value
	 * is determined when a new OS entry is added to a workspace.
	 * This value is persisted in the workspace and restored when
	 * a workspace is loaded.
	 */
	private OSType osType;
	
	/**
	 * This is a transient field that is never persisted (for security
	 * purposes). Typically, it is set once (each time the GUI is run)
	 * when an attempt is made to access the server.
	 */
	private transient String password;	
}