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
import java.util.ArrayList;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * A class to encapsulate information about a list of servers that have already
 * been configured in this work space. This class is instantiated from the Work
 * space class. This class is relatively straightforward in that it merely
 * contains a list of Servers. In addition, it facilitates marshalling and
 * unmarshalling of server data.
 * 
 */
public class ServerList {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * ServerList entry. This method is typically  used to create a suitable
	 * ServerList entry using data from a given DOM element. For each 
	 * Server entry this method uses the Server.createServer method to
	 * create suitable Server entries.
	 * 
	 * @param serverListNode The DOM element to be used for creating the 
	 * server list and to be used for creating the Server entries.
	 * 
	 * @return The newly created ServerList entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static ServerList create(Element serverNode) throws Exception {
		// Create the ServerList we are going to populate below.
		ServerList srvrList = new ServerList();
		// First extract the sequence counter information
		String seqCounter = DOMHelper.getStringValue(serverNode, "SeqCounter");
		srvrList.seqCounter = Integer.parseInt(seqCounter);
		// Now obtain all the Server elements under the serverList
		NodeList servers = serverNode.getElementsByTagName("Server");
		for(int idx = 0; (idx < servers.getLength()); idx++) {
			Node   tmpNode = servers.item(idx);
			// Ensure that this node is actually a server node
			if ((tmpNode.getNodeType() == Node.ELEMENT_NODE) &&
				("Server".equals(tmpNode.getNodeName())) && 
				(tmpNode instanceof Element)) {
				// OK, safe to type cast and try to parse..
				Element srvrNode = (Element) tmpNode;
				Server srvr = Server.create(srvrNode);
				// Add the valid server node to the list.
				srvrList.servers.add(srvr);
			}
		}
		// Return the newly created list.
		return srvrList;
	}
	
	/**
	 * The default constructor. It merely initializes all the instance
	 * variables to their default initial value.
	 */
	public ServerList() {
		servers    = new ArrayList<Server>();
		seqCounter = 1;
	}

	/**
	 * Reserves the next server ID for use in a server object. This method is
	 * used only when a new server object is added by the user. Server objects
	 * loaded from a file have their own server IDs.
	 * 
	 * @return A new, unique (within the work space) server ID for use with
	 * a new Server object instance.
	 */
	public String reserveServerID() {
		// Use current sequence counter value to create ID
		String id = "server" + seqCounter;
		// Setup sequence counter for next server.
		seqCounter++;
		// Return the ID for use by the caller.
		return id;
	}
	
	/**
	 * Method to add a new server entry to this server list.
	 * 
	 * This method sets a new server ID for the entry and 
	 * fires a WorkspaceEvent indicating the addition of the new
	 * server entry to all listeners.
	 * 
	 * @param server The new server entry to be added.
	 */
	public synchronized void add(Server server) {
		// First setup a new ID for this entry.
		server.setID(reserveServerID());
		// Add the entry
		servers.add(server);
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(server, WorkspaceEvent.Operation.INSERT);
		Workspace.get().fireWorkspaceChanged(we);
	}
	
	
	/**
	 * Method to remove a server entry to this server list.
	 * 
	 * This method must be used to remove a server entry from the
	 * server list. Once the entry has been removed, this method 
	 * fires a WorkspaceEvent indicating the removal of the
	 * server entry to all listeners.
	 * 
	 * @param server The new server entry to be removed.
	 */
	public synchronized void remove(Server server) {
		if (servers.remove(server)) {
			// Fire notification to listeners to update GUIs
			WorkspaceEvent we = new WorkspaceEvent(server, WorkspaceEvent.Operation.DELETE);
			Workspace.get().fireWorkspaceChanged(we);
		}
	}
	
	
	/**
	 * Obtain reference to a Server object, given its work space wide unique
	 * server ID. This method does a linear search on the serverList to
	 * locate the server entry (if you have 1000s of server entries in one
	 * workspace then you are missing the concept of a workspace). 
	 * 
	 * @param serverID The generated work space wide unique server ID. 
	 * Note that checks are case sensitive. 
	 * 
	 * @return Reference to a Server object whose serverID is the same as the
	 * supplied serverID. If a matching object was not found then this method
	 * returns null.
	 */
	public Server getServer(String serverID) {
		for(Server srvr: servers) {
			if (srvr.getID().equals(serverID)) {
				return srvr;
			}
		}
		// No matching entry found.
		return null;
	}
	
	/**
	 * Obtain the list of Server objects that are currently available
	 * in a work space.
	 * 
	 * @return The list of Server objects that are currently available
	 * in this work space.
	 */
	public ArrayList<Server> getServers() { return servers; }
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent Workspace node in the DOM tree.
	 * 
	 * @param workspace The DOM element corresponding to the Workspace
	 * node that contains this entry.
	 */
	public final void marshall(Element workspace) {
		// Create a top-level server list entry for this class
		Element srvrList = DOMHelper.addElement(workspace, "ServerList", null);
		// Add the sequence counter information to server list
		DOMHelper.addElement(srvrList, "SeqCounter", "" + seqCounter);
		// Add sub-elements for each server in our server list
		for (Server srvr : servers) {
			srvr.marshall(srvrList);
		}
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragement. The XML fragement is guaranteed to be compatible
	 * with the PEACE work space configuration data. This method provides
	 * better control on the XML formatting to generate a more readable
	 * XML output when compared to the DOM tree.
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t";
		// Create a top-level server list entry for this class
		out.printf("%s<ServerList>\n", Indent);
		// Add the sequence counter information to server list
		out.printf("%s\t<%2$s>%3$s</%2$s>\n", Indent, "SeqCounter", seqCounter);
		// Add sub-elements for each server in our server list
		for (Server srvr : servers) {
			srvr.marshall(out);
		}
		// Close the server list element.
		out.printf("%s</ServerList>\n\n", Indent);
	}
	
	/**
	 * The list of Server objects that have been configured and added
	 * to this list. This list is populated either via the create 
	 * method or new entries are added via the addServer() method
	 * in this class. 
	 */
	private ArrayList<Server> servers;
	
	/**
	 * Sequence counter that is maintained on a per-work space basis to
	 * generate unique/valid IDs for each new server entry added to this
	 * list.
	 */
	private long seqCounter;
}
