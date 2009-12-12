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

package org.peace_tools.core.session;

import java.awt.Component;

import org.peace_tools.workspace.Server;

/**
 * A helper factory class to create a suitable session depending on server.
 * 
 * This is a helper factory class that is used to create a local or a 
 * remote server session to interact with the server. The objective of
 * this class is to provide a single point for instantiating suitable
 * server session. 
 * 
 * @note This class is not meant to be instantiated. Instead directly use
 * the static methods in this class.
 */
public class SessionFactory {
	/**
	 * Factory method to create suitable server session.
	 * 
	 * This method must be used to create a suitable server session to
	 * interact with a remote server.   This method performs no special
	 * operations (such as connecting etc.) on the session.
	 * 
	 * @param parent The parent component that should be used to create 
	 * GUI elements that may be needed for any interactive operations.
	 * 
	 * @param server The server data-object that provides the necessary
     * information to connect to the server.
     * 
	 * @return Returns a suitable (local or remote) session for interacting
	 * with the server.
	 */
	public static ServerSession createSession(Component parent, Server server) {
		if (server.isRemote()) {
			return new RemoteServerSession(parent, server);
		} else {
			return new LocalServerSession(parent, server);
		}
	}
	
	/**
	 * Default constructor. The default constructor has been made private
	 * to ensure this class is never directly instantaited.
	 */
	private SessionFactory() {}
}
