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

package org.peace_tools.decagon;

/**
 * This class encapsulates the shared property keys associated with various
 * DECAGON configuration properties.
 * 
 * This class is a non-instantiable class that contains the keys that
 * are used by various classes constituting the DECAGON-GUI. The property
 * keys are used to obtain or update properties in the global properties
 * managed by PEACE.
 */
public class PropertyKeys {
	/**
	 * A convenience enumeration to make server configuration information
	 * and associated code more readable and maintainable. 
	 * 
	 * <p>This information is used to store the type of server entry being
	 * used for a specific purpose in DECAGON. This information is 
	 * stored as strings in the global properties managed via the 
	 * PEACEProperties class. For example, this value is used for
	 * the following configuration variables: 
	 * {@link PropertyKeys#STORAGE_SERVER_TYPE}, 
	 * {@link PropertyKeys#MYSQL_TUNNEL_SERVER_TYPE}, etc.</p>
	 *
	 * The enumeration values are interpreted in the following manner:
	 * 
	 * <ul>
	 * <li>If server type is {@value #LOCAL_HOST} then it indicates the
	 * local host on which PEACE is currently running.</li> 
	 * 
	 * <li>If the server type is {@value #KNOWN_HOST} then it indicates
	 * the server is an entry in the global ServerList in this
	 * workspace.</li>
	 * 
	 * <li>If the server type is {@link #CUSTOM_HOST} then this indicates
	 * a custom information about the machine.</li>
	 * 
	 * </ul>
	 */
	public enum ServerEntryType { LOCAL_HOST, KNOWN_HOST, CUSTOM_HOST};
	
	/**
	 * The property key to determine if Metasim is enabled or disabled.
	 * The value for this key is a boolean literal (either TRUE or FALSE).
	 */
	public static final String METASIM_ENABLED = "DECAGON.MetasimEnabled";

	/**
	 * The property key to obtain the path to the Metasim executable
	 * path.
	 */
	public static final String METASIM_PATH = "DECAGON.MetasimPath";

	/**
	 * The property key to determine if shared storage has been enabled.
	 * The value for this key is a boolean literal (either TRUE or FALSE).
	 */
	public static final String SHARED_SOTRAGE_ENABLED = "DECAGON.SharedStorageEnabled";

	/**
	 * The property key to obtain the type of server setting that was selected
	 * by the user when the storage was configured. The value for this string
	 * must be from the ServerEntryType enumeration.
	 */
	public static final String STORAGE_SERVER_TYPE = "DECAGON.StorageServerType";
	
	/**
	 * The property key to obtain the unique identifier for a known host to be
	 * used for the storage server. This value is meaningful only if the
	 * {@link #STORAGE_SERVER_TYPE} is set to {@value ServerEntryType#KNOWN_HOST}.
	 */
	public static final String STORAGE_SERVER_KNOWN_HOST_KEY = "Decagon.StorageServerKnownHostKey";
	
	/**
	 * The property key to obtain the host name/IP address/ID of the
	 * server that contains the shared storage space to be used by
	 * DECAGON. This value can be one of the following depending
	 * on the type of server set in {@value #STORAGE_SERVER_TYPE}:
	 * 
	 * <ul>
	 * <li>If server type is {@value ServerEntryType#LOCAL_HOST} then this 
	 * value is simply <code>localhost</code>.
	 * <li>If the server type is {@link ServerEntryType#CUSTOM_HOST} or
	 * {@value ServerEntryType#KNOWN_HOST} then this value
	 * is the actual hostName/IP address.</li>
	 * </ul>
	 * 
	 */
	public static final String STORAGE_SERVER_NAME = "DECAGON.StorageServerName";

	/**
	 * The property key to obtain the port number of the server that 
	 * contains the shared storage space to be used by DECAGON. This value is
	 * meaningful only if {@link #STORAGE_SERVER_TYPE} is set to
	 * {@value ServerEntryType#CUSTOM_HOST} or {@value ServerEntryType#KNOWN_HOST}.
	 */
	public static final String STORAGE_SERVER_PORT = "DECAGON.StorageServerPort";

	/**
	 * The property key to obtain the user ID to be used to log onto the
	 * server that contains the shared storage space to be used by
	 * DECAGON. This value is
	 * meaningful only if {@link #STORAGE_SERVER_TYPE} is set to
	 * {@value ServerEntryType#CUSTOM_HOST} or {@value ServerEntryType#KNOWN_HOST}. 
	 */
	public static final String STORAGE_SERVER_USERID = "DECAGON.StorageServerUserID";

	/**
	 * The property key to obtain the actual shared directory which serves
	 * as the shared storage space. This value is used for all types of
	 * storage server entries. This value is the absolute path on the 
	 * appropriate storage server.
	 */
	public static final String STORAGE_SERVER_PATH = "DECAGON.StorageServerPath";

	/**
	 * The property key to determine if a SSH-tunnel to access MySQL 
	 * database is enabled or disabled. The value for this key is a 
	 * boolean literal (either TRUE or FALSE).
	 */
	public static final String MYSQL_TUNNEL_ENABLED = "DECAGON.MySqlTunnelEnabled";

	/**
	 * The property key to obtain the type of server setting that was selected
	 * by the user when the storage was configured.  The value for this string
	 * must be from the ServerEntryType enumeration.
	 */
	public static final String MYSQL_TUNNEL_SERVER_TYPE = "DECAGON.MySqlTunnelServerType";
	
	/**
	 * The property key to obtain the unique identifier for a known host to be
	 * used for the storage server. This value is meaningful only if the
	 * {@link #STORAGE_SERVER_TYPE} is set to {@value ServerEntryType#KNOWN_HOST}.
	 */
	public static final String MYSQL_TUNNEL_SERVER_KNOWN_HOST_KEY = "Decagon.MySqlTunnelServerKnownHostKey";

	/**
	 * The property key to obtain the host name/IP address/ID of the
	 * server to be used to establish a SSH-tunnel to MySQL database.
	 * This value can be one of the following depending  on the type
	 * of server:
	 * 
	 * <ul>
	 * <li>If server type is {@value ServerEntryType#LOCAL_HOST} then this 
	 * value is simply <code>localhost</code>.
	 * <li>If the server type is {@link ServerEntryType#CUSTOM_HOST} or
	 * {@value ServerEntryType#KNOWN_HOST} then this value
	 * is the actual hostName/IP address.</li>
	 * </ul>
	 */
	public static final String MYSQL_TUNNEL_SERVER_NAME = "DECAGON.MySqlTunnelServerName";

	/**
	 * The property key to obtain the port number on the server to which 
	 * initial SSH-connection is to be established. This value is
	 * meaningful only if the {@link #MYSQL_TUNNEL_SERVER_TYPE} server 
	 * type is {@value ServerEntryType#KNOWN_HOST} or 
	 * {@value ServerEntryType#CUSTOM_HOST}.
	 */
	public static final String MYSQL_TUNNEL_SERVER_PORT = "DECAGON.MySqlTunnelServerPort";

	/**
	 * The property key to obtain the user ID to be used to log onto the
	 * server that contains the shared storage space to be used by
	 * DECAGON. This value is meaningful only if the 
	 * {@link #MYSQL_TUNNEL_SERVER_TYPE} server type is 
	 * {@value ServerEntryType#KNOWN_HOST} or 
	 * {@value ServerEntryType#CUSTOM_HOST}.
	 */
	public static final String MYSQL_TUNNEL_SERVER_USERID = "DECAGON.MySqlTunnelUserID";

	/**
	 * The property key to determine if MySQL database has been 
	 * enabled. The value for this key is a boolean literal (either
	 * TRUE or FALSE).
	 */
	public static final String MYSQL_ENABLED = "DECAGON.MySqlEnabled";

	/**
	 * The property key to obtain the host name/IP address/ID of the
	 * database host. This value is the actual hostName/IP address.
	 * 
	 */
	public static final String MYSQL_SERVER_NAME = "DECAGON.MySqlServerName";

	/**
	 * The property key to obtain the port number on the server through
	 * which connection to MySQL is to be established.
	 */
	public static final String MYSQL_SERVER_PORT = "DECAGON.MySqlServerPort";

	/**
	 * The property key to obtain the database user ID to be used to 
	 * log onto the MySQL server.
	 */
	public static final String MYSQL_SERVER_USERID = "DECAGON.MySqlDBUserID";

	/**
	 * The property key to obtain the actual name of the database that is
	 * being used to store data.
	 */
	public static final String MYSQL_DB_NAME = "DECAGON.MySqlDBName";

	/**
	 * The default constructor.
	 * 
	 * The constructor has been made private to ensure that this class
	 * is never instantiated. Instead the static property keys in this
	 * class can be directly accessed and used.
	 */
	private PropertyKeys() {
		// Nothing to be done in the constructor.
	}
}
