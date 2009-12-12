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

import java.util.Map;

/**
 * Top-level logging class for creating programmer logs.
 * 
 * <p>The ProgrammerLog provides a convenient mechanism the GUI system 
 * to log relevant information to help programmers debug and troubleshoot
 * problems. This log is essentially used to record information that may
 * be useful for a programmer but not useful for the user. <p>
 * 
 * <p>This class is a centralized class for merely recording the logs (and 
 * not viewing/seeing) the logs. The logs are viewed/shown by the 
 * Programmer LogPane which constitutes the "View" as in the Model-
 * View-Controller (MVC) pattern. This class represents the "Model" in 
 * the MVC pattern while the top-level frame that contains the 
 * ProgrammerLogPane serves as the controller. </p>
 * 
 * Since the ProgrammerLog serves as an centralized location for cutting 
 * logs relevant/useful to the programmer, there is only one unique 
 * instance of this class in the GUI system. In order to enforce this 
 * property, this class has been designed as a Singleton pattern -- 
 * that is, this class cannot be directly instantiated. Instead the 
 * getLog() static method must be used to refer to the globally unique
 * instance of this class. On the other the convenient, static log() 
 * method can be used to directly cut log entries. 
 * 
 * @note Logs for user (that would be pertinent for the user to view and
 * take action) must be generated in the UserLog.
 */
public class ProgrammerLog extends Log {
	/** Obtain reference to the globally unique instance of the
	 * ProgrammerLog.
	 * 
	 * This method must be used to obtain a reference to the globally 
	 * unique instance of the UserLog. There is only one singleton 
	 * instance of the ProgrammerLog that must be used to generate 
	 * all programmer-log entries.
	 * 
	 * @return The globally unique instance of the UserLog.
	 */
	public static ProgrammerLog getLog() {
		return programmerLog;
	}
	
	/**
	 * This method logs the specified message in the programmer log.
	 * After the log 
	 * 
	 * @param text
	 *            The text or actual log message.
	 *            
	 * @return Logs the message and returns true.
	 */
	public static synchronized boolean log(String text) {
		programmerLog.logEntries.append(text);
		// Write data to a log file if set.
		if (programmerLog.logFile != null) {
			programmerLog.logFile.println(text);
		}
		// Let the listeners know that a new log entry has been cut.
		programmerLog.notifyViews(true, false, false);
		// Report that the log was cut
		return true;
	}
	
	/**
	 * This method logs the specified message in the programmer log.
	 * After the log 
	 * 
	 * @param e Logs the exception in the programmer log.
	 *            
	 * @return Logs the message and returns true.
	 */
	public static synchronized boolean log(Exception e) {
		String text = e.toString();
		programmerLog.logEntries.append(text);
		// Write data to a log file if set.
		if (programmerLog.logFile != null) {
			programmerLog.logFile.println(text);
		}
		// Let the listeners know that a new log entry has been cut.
		programmerLog.notifyViews(true, false, false);
		// Report that the log was cut
		return true;
	}
	
	/** Obtain the set of log entries contained in the user logs.
	 * 
	 * @return The current set of log entries contained in this class.
	 */
	public StringBuffer getLogEntries() {
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
	private ProgrammerLog() {
		// Setup the default log level and log file name.
		super(Log.LogLevel.NOTICE, null);
		// Create the in-memory buffer for programmer logs
		logEntries = new StringBuffer();
	}
	
	/**
	 * This instance variable holds all the user logs that have been 
	 * generated thus far via the log() method in this class. The 
	 * entries added are displayed by various UserLogPane classes. 
	 * The entries are stored in a compact fashion to minimize the 
	 * in-memory size. The UserLogPane class uses a suitable
	 * renderer to appropriately display the log entries. 
	 */
	private StringBuffer logEntries;
	
	/**
	 * The globally unique instance of the UserLog. This singleton 
	 * instance is used to generate and store all of the user logs
	 * displayed by the various UserLogPane objects in various 
	 * windows.
	 */
	private static final ProgrammerLog programmerLog;
	
	/**
	 * Static initializer block to create the globally unique 
	 * user log and cut some initial logs in it.
	 */
	static {
		programmerLog = new ProgrammerLog();
		// Cut some initial logs in the programmer log.
		log("PEACE GUI Infrastructure started up.\n" + 
			"Programmer logs are buffered.\n");
		// Cut the environment in which the GUI is being run.
		StringBuffer envData = new StringBuffer();
		envData.append("PEACE is running in enviroment: \n");
		Map<String, String> env = System.getenv();
		for(Map.Entry<String, String> entry: env.entrySet()) {
			envData.append(entry.getKey());
			envData.append("='");
			envData.append(entry.getValue());
			envData.append("'\n");
		}
		log(envData.toString());
	}
}
