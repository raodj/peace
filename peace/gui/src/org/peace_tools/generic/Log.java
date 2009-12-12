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

import java.io.FileNotFoundException;
import java.io.PrintStream;
import javax.swing.event.EventListenerList;

/**
 * <p>
 * This is an abstract base class that provides the core features for managing
 * log entries. The logging system has been developed using the
 * Model-View-Controller (MVC) pattern. This class represents the core component
 * of the "Model" that provides the data needed by the "Views" for display
 * purposes. The corresponding "View" class is the LogPane. The top-level Frame
 * that contains the views acts as the "Controller".</p>
 * 
 * <p> This base class provides the following features:
 * 
 * <ul>
 * 
 * <li>It provides a generic place holder that is used by UserLog and
 * ProgrammerLog classes to store log information.</li>
 * 
 * <li>It provides a mechanism for enabling and disabling saving the logs to a
 * given log file.</li>
 * 
 * </ul>
 *
 * </p>
 */
public abstract class Log {
	/**
	 * This enumeration provides a convenient mechanism for referring to
	 * different error levels when generating logs. The error levels are
	 * interpreted in different manners as documented with each enumeration. The
	 * user typically has the option to control the level of logging that the
	 * user desires. Typically {@code DEBUG} and {@code INFO} levels are not
	 * logged by default.
	 */
	public enum LogLevel {
		/**
		 * This is an unimportant log entry that is typically used to provide
		 * additional information to the user.
		 */
		INFO, // Informational level
		/**
		 * This is a reasonably important log entry that is typically used to
		 * provide additional information to the user. This level is used to
		 * provide some useful information to the user (such as timings or Job
		 * completion notices) that does not impact the overall operation of the
		 * system.
		 */
		NOTICE, // Some useful log
		/**
		 * This log level value must be used to generate warning messages that
		 * indicate impeding problems to the user. Messages such as loss of
		 * network connectivity or reduced availability of memory fall into this
		 * category.
		 */
		WARNING, // Foreboding/Impeding problems
		/**
		 * This log level value must be used to generate error messages that
		 * report failure of operations and provide additional information
		 * regarding failure status.
		 */
		ERROR // Failures
	}

	/**
	 * The current log level at which log messages are being generated.
	 * Typically this value is set to LogLevel.NOTICE.
	 */
	protected LogLevel currentLogLevel = LogLevel.NOTICE;

	/**
	 * The file we are writing to (if any). The log file is typically set each
	 * time the program is launched.
	 */
	protected String logFileName = null;

	/**
	 * The log file stream to which the data is currently being written as logs
	 * are being generated. The file is opened and closed depending on whether
	 * the log is to be persisted or not.
	 */
	protected PrintStream logFile = null;

	/**
	 * The list of log listeners that are interested in receiving call back
	 * information whenever log information changes. 
	 */
	protected EventListenerList listeners;
	
	/**
	 * A constructor for the Log. It merely initializes all the instance
	 * variables to their default initial value. 
	 * 
	 * @param logLevel
	 *            The initial/default Log level for the log messages.
	 * @param fileName
	 *            An optional log file name to log to. This can be null.
	 */
	protected Log(LogLevel logLevel, String fileName) {
		// Setup the list to maintain the listeners
		listeners = new EventListenerList();
		// Set the preferred log level
		setLevel(logLevel);
		// Finally set the file name specified.
		setFileName(fileName);
	}

	/**
	 * Sets the Logging level to the specified level.
	 * 
	 * @param level
	 *            The new level at which log entries must be recorded.
	 */
	public void setLevel(LogLevel level) {
		if (currentLogLevel != level) {
			currentLogLevel = level;
			// Let views know about change is log level
			notifyViews(false, false, true);
		}
	}

	/**
	 * This is a helper method that is used to dispatch a call back to
	 * all interested listeners about the change in log information. This
	 * method is invoked by other methods in this class whenever the log
	 * information changes. This method walks the list of LogListener
	 * objects in the listeners list and invokes the logChanged() API
	 * method to dispatch notifications. 
	 * 
	 * @param logsChanged This parameter is true if log entries were added
	 * or removed from the logs. Typically this is the common scenario for
	 * the call back to occur.
	 * 
	 * @param fileNameChanged This parameter is true if the file name into
	 * which logs are written was changed by the user (via some view).
	 * 
	 * @param levelChanged This parameter is true if the current logging
	 * level set for the log was changed by the user (via some view).
	 */
	protected synchronized void notifyViews(final boolean logChanged, 
			final boolean fileNameChanged, final boolean levelChanged) {
		// Obtain the list of compatible listeners/views
		LogListener views[] = listeners.getListeners(LogListener.class);
		// Let each view know that the log status has changed.
		for(LogListener view : views) {
			view.logChanged(logChanged, fileNameChanged, levelChanged);
		}
	}
	
	/**
	 * Adds a log listener (typically a view derived off of LogPane) to
	 * the set of listeners to be notified whenever log information
	 * changes.
	 * 
	 * @param listener The listener to be notified on changes to the
	 * log.
	 */
	public synchronized void addLogListener(LogListener listener) {
		listeners.add(LogListener.class, listener);
	}
	
	/**
	 * Removes a log listener (typically a view derived off of LogPane)
	 * from the set of listeners to be notified whenever log information
	 * changes.
	 * 
	 * @param listener The listener to be removed from the list of 
	 * listeners interested in receiving log updates.
	 */
	public synchronized void removeLogListener(LogListener listener) {
		listeners.remove(LogListener.class, listener);
	}
	
	/**
	 * Retrieves the current logging level.
	 * 
	 * @return The current logging level at which log entries are being
	 *         recorded. Log entries below this specified level of severity are
	 *         not recorded.
	 */
	public LogLevel getLevel() {
		return currentLogLevel;
	}

	/**
	 * Parses a string containing a valid logging level and sets the level
	 * appropriately.
	 * 
	 * @param level
	 *            The String to parse and convert to a suitable Enum value. This
	 *            method throws an IllegalArgumentException if the level string
	 *            is an invalid log level.
	 */
	public void parseLevel(String level) {
		setLevel(LogLevel.valueOf(level));
	}

	/**
	 * This method clears the log and writes a special message to the
	 * PrintStream, if one exists.
	 */
	public void clearLog() {
		if (logFile != null) {
			logFile.println("Log cleared.");
		}
	}

	/**
	 * Set a file to which user logs are to be written for persistent storage.
	 * 
	 * This method can be used to set a file to which the user logs are to be
	 * written for future reference.
	 * 
	 * @note This method closes the existing log if one is open already.
	 * 
	 * @param file
	 *            The name of the file to which the user logs are to be written
	 *            for permanent storage. If this parameter is null then any open
	 *            log stream is closed and a new file is not opened.
	 *            
	 * @return This  method returns true if the operation was sucessful.
	 * On errors this method returns false.
	 */
	public boolean setFileName(String file) {
		if (logFile != null) {
			logFile.close();
			logFile = null;
		}
		if (file == null) {
			// Nothing to be done.
			logFileName = null;
			notifyViews(false, true, false);
			return true;
		}
		try {
			logFileName = file;
			logFile     = new PrintStream(file);
			// Now that we have opened the log file for writing
			// let the views know the file name has changed.
			notifyViews(false, true, false);
		} catch (FileNotFoundException e) {
			// Dang! That is not good
			logFileName = null;
			logFile     = null;
		}
		// When control drops here either the file was successfully
		// opened, or some error occurred. Let the views know.
		notifyViews(false, true, false);
		// Return the resulting status.
		return (logFile != null);
	}

	/**
	 * Obtain the name of the file to which the user logs are being written. If
	 * a valid file name is not set then this method returns null.
	 * 
	 * @return The name of the file to which user logs are being written.
	 */
	public String getFileName() {
		return logFileName;
	}
	
	/**
	 * Determine if logs are actually being written to a log file.
	 * 
	 * This method can be used to determine if the logs are currently being
	 * written to a log file.
	 * 
	 * @return This method returns true if the logs are being written to
	 * the file indicated by logFileName.
	 */
	public boolean isWritingLogs() {
		return (logFile != null);
	}

	/**
	 * method to start/stop writing logs to the given log file.
	 * 
	 * This method can be used to enable/disable writing data to a log
	 * file.
	 */
	public void toggleLogWriting() {
		if (logFile == null) {
			logFile.close();
			logFile = null;
		} else if ((logFile == null) && (logFileName != null)) {
			setFileName(logFileName);
		}
	}
	/**
	 * Utility method to convert a log level to a single character value.
	 * This method is used by the UserLog to combine log entries together
	 * to store information in a more compact form.
	 * 
	 * @param level The level that must be encoded to a single character.
	 * @return A single character representing the log level.
	 */
	public static char encode(LogLevel level) {
		return level.toString().charAt(0);
	}

	/**
	 * Utility method to convert an encoded log level character (obtained
	 * via a call to the encode() method) to corresponding LogLevel
	 * enumeration value. This method is used by the UserLog to combine log
	 * entries together to store information in a more compact form.
	 * 
	 * @param level A single character level code to be translated to 
	 * corresponding log level.
	 * 
	 * @return The log level enumeration representing the single character
	 * encoding for the log level.
	 */
	public static LogLevel decode(char code) {
		switch (code) {
		case 'I': return LogLevel.INFO;
		case 'N': return LogLevel.NOTICE;
		case 'W': return LogLevel.WARNING;
		case 'E':
		default: return LogLevel.ERROR;
		}
	}
}
