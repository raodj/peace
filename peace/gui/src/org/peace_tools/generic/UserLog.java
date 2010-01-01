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

package org.peace_tools.generic;

import java.util.Date;
import javax.swing.table.DefaultTableModel;

/**
 * Top-level logging class for creating user logs.
 * 
 * <p>
 * The UserLog provides a convenient mechanism the GUI system to log relevant
 * information at various levels of severity. This log is essentially used to
 * provide the user with useful information. Note that this class is a
 * centralized class for merely recording the logs (and not viewing/seeing) the
 * logs. The logs are viewed/shown by the UserLogPane which constitutes the
 * "View" as in the Model-View-Controller (MVC) pattern. This class represents
 * the "Model" in the MVC pattern while the top-level frame that contains the
 * UserLogPane serves as the controller.
 * </p>
 * 
 * Since the UserLog serves as an centralized location for cutting logs
 * relevant/useful to the user, there is only one unique instance of this class
 * in the GUI system. In order to enforce this property, this class has been
 * designed as a Singleton pattern -- that is, this class cannot be directly
 * instantiated. Instead the getLog() static method must be used to refer to the
 * globally unique instance of this class. On the other the convenient, static
 * log() method can be used to directly cut log entries.
 * 
 * <p>
 * <b>Note:</b>Logs for programming aids and debugging (that would be pertinent
 * for a programmer to view and understand) must be generated in the free form
 * ProgrammerLog and not in the UserLog.
 * </p>
 */
public class UserLog extends Log {
	/** Obtain reference to the globally unique instance of the UserLog.
	 * 
	 * This method must be used to obtain a reference to the globally unique
	 * instance of the UserLog. There is only one singleton instance of the
	 * UserLog that must be used to generate all user-log entries.
	 * 
	 * @return The globally unique instance of the UserLog.
	 */
	public static UserLog getLog() {
		return userLog;
	}
	
	/**
	 * This method logs the specified message at the given level. If the 
	 * level is valid (higher than the current setting) then the message
	 * is logged and true returned. Otherwise, false is returned.
	 * 
	 * @param level 
	 *            The level of the log to be cut.
	 * @param component
	 *            The component or subsystem that is generating the log message.
	 * @param text
	 *            The text or actual log message.
	 *            
	 * @return True if logged, false otherwise
	 */
	public static synchronized boolean log(LogLevel level, String component,
			String text) {
		if (userLog.currentLogLevel.compareTo(level) > 0) {
			// This log is not at a sufficient level of severity to log
			return false;
		}
		
		// Obtain the time when this log is being generated.
		Date logTime = new Date();
		// Create the encoded level+timestamp string.
		final String timestamp = "" + Log.encode(level) + logTime;
		// Ensure we have a valid component to work with.
		if (component == null) {
			component = "N/A";
		}
		// Create the data to be logged.
		final String[] toLog = { timestamp, component, text };
		userLog.logEntries.addRow(toLog); // Cut the log.
		// Write data to a log file if set.
		if (userLog.logFile != null) {
			String logLine = toLog[0] + '\t' + toLog[1] + '\t' 
					+ toLog[2] + '\t' + toLog[3];
				userLog.logFile.println(logLine);
		}
		// Let the listeners know that a new log entry has been cut.
		userLog.notifyViews(true, false, false);
		// Report that the log was cut
		return true;
	}
	
	/** Obtain the set of log entries contained in the user logs.
	 * 
	 * @return The current set of log entries contained in this class.
	 */
	public DefaultTableModel getLogEntries() {
		return logEntries;
	}
	
	/**
	 * The default (and only) constructor.
	 * 
	 * The UserLog represents a centralized facility for cutting logs 
	 * relevant/useful to the user. In other words, there must be only
	 * one, globally unique instance of this class that contains all the
	 * user logs in one centralized location. In order to ensure that
	 * only one unique instance is ever created, the constructor has
	 * been made private and a globally unique static variable is used
	 * to hold a singleton reference. The reference object can be 
	 * accessed via the getLog() method.
	 */
	private UserLog() {
		// Setup the default log level and log file name.
		super(Log.LogLevel.NOTICE, null);
		// Create the in-memory list of user logs
		final String[] Titles = { "Timestamp", "Component", "Message" };
		logEntries = new DefaultTableModel(Titles, 0);
	}
	
	/**
	 * This instance variable holds all the user logs that have been 
	 * generated thus far via the log() method in this class. The 
	 * entries added are displayed by various UserLogPane classes. 
	 * The entries are stored in a compact fashion to minimize the 
	 * in-memory size. The UserLogPane class uses a suitable
	 * renderer to appropriately display the log entries. 
	 */
	private DefaultTableModel logEntries;
	
	/**
	 * The globally unique instance of the UserLog. This singleton 
	 * instance is used to generate and store all of the user logs
	 * displayed by the various UserLogPane objects in various 
	 * windows.
	 */
	private static final UserLog userLog;
	
	/**
	 * Static initializer block to create the globally unique 
	 * user log and cut some initial logs in it.
	 */
	static {
		userLog = new UserLog();
		// Cut some initial logs
		// log(Log.LogLevel.NOTICE, "PEACE", Version.GUI_VERSION);
		// log(Log.LogLevel.NOTICE, "PEACE", Version.COPYRIGHT);
		// log(Log.LogLevel.NOTICE, "PEACE", Version.DISCLAIMER);
	}
}
