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

package org.peace_tools.core;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Properties;

import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.peace_tools.generic.Log;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;

/**
 * Helper class to manage process-wide, persist properties for PEACE GUI.
 * 
 * <p>This class provides a convenient wrapper around Java's default 
 * java.util.Properties class. This class provides a singleton instance 
 * for the properties and is used by various classes in PEACE GUI to 
 * manage their persistent properties. The {@link #get()} method in this
 * class should be used to obtain the process-wide unique instance of
 * this class for use.<p>
 * 
 * <p>All properties are stored as <name, value> pairs. The name associated
 * with all properties must be unique. Therefore it is required that each c
 * lass that stores properties in this file should include its class 
 * name as part of the property key (example: 
 * <code>WorkspaceChooer.WorkspaceList</code> would be a key in the 
 * properties list rather than simply <code>WorkspaceList</code>).
 * Note that the name and value must be standard Java Strings.
 * </p>
 */
public class PEACEProperties {
	/**
	 * The global set of persistent properties.
	 * 
	 * This instance variable maintains the process-wide, persistent
	 * set of properties for PEACE. The properties are shared across
	 * workspaces. The properties are stored as an XML file in
	 * the User's default home (obtained by {@link Utilities#getDefaultDirectory()}
	 * method) in file whose name is given by {@link #PEACE_PROPERTIES_FILE}  
	 */
	private static final PEACEProperties peaceProperties = new PEACEProperties();
	
	/**
	 * Static method to obtain access to the shared, process-wide unique 
	 * properties object.
	 * 
	 * This method must be used to obtain a reference to the globally
	 * shared singleton instance of this class. The object returned
	 * by this method must be used to get/set properties. In addition,
	 * the returned object can be used to save and load the properties.
	 * 
	 * @return The process-wide unique instance of this class that contains
	 * the current set of properties used by PEACE.
	 */
	public static PEACEProperties get() {
		return peaceProperties;
	}
	
	/**
	 * Obtain the value for a given property name (if present in the properties)
	 * 
	 * Searches for the property with the specified key in this property list. 
	 * If the key is not found in this property list, the default property list, 
	 * (if any) is then checked. The method returns null if the property is not 
	 * found.
	 *  
	 * @param name The name (or key) for the property. This value must be unique.
	 * 
	 * @return The value associated with the given property name if any. If
	 * a property with the given name was not found, then this method returns
	 * null.
	 */
	public synchronized String getProperty(final String name) {
		return properties.getProperty(name);
	}
	
	/**
	 * Obtain the value for a given property name (if present in the properties)
	 * 
	 * Searches for the property with the specified key in this property list. 
	 * If the key is not found in this property list, the default property list, 
	 * (if any) is then checked. This method returns the defValue if the property
	 * is not found.
	 *  
	 * @param name The name (or key) for the property. This value must be unique.
	 * 
	 * @param defValue The default value to be returned by this method if
	 * the property was not found in this list.
	 * 
	 * @return The value associated with the given property name if any. If
	 * a property with the given name was not found, then this method returns
	 * the default value.
	 */
	public synchronized String getProperty(final String name, final String defValue) {
		return properties.getProperty(name, defValue);
	}

	/**
	 * Obtain the boolean value for a given property name (if present in the properties)
	 * 
	 * Searches for the property with the specified key in this property list. 
	 * If the key is not found in this property list, the default property list, 
	 * (if any) is then checked. This method returns the defValue if the property
	 * is not found.
	 *  
	 * @param name The name (or key) for the property. This value must be unique.
	 * 
	 * @param defValue The default value to be returned by this method if
	 * the property was not found in this list.
	 * 
	 * @return The value associated with the given property name if any. If
	 * a property with the given name was not found, then this method returns
	 * the default value.
	 */
	public synchronized boolean getBooleanProperty(final String name, final boolean defValue) {
		final String strValue = properties.getProperty(name);
		return (strValue == null) ? defValue : Boolean.parseBoolean(strValue);
	}

	/**
	 * Remove a property from the properties list.
	 * 
	 * This method must be used to remove a property from the properties
	 * list. If the property does not exist then this method has no
	 * side effects.
	 * 
	 * @param name The name of the property to be removed. This value
	 * cannot be null and must uniquely identify the property to be 
	 * removed.
	 */
	public synchronized void removeProperty(final String name) {
		properties.remove(name);
	}
	
	/**
	 * Set/change the value for a given property name.
	 * 
	 * This method must be used to set (initially) or change the value
	 * associated with a given property. The name of properties must be
	 * unique. No checks are made to ensure property names are unique.
	 * Therefore care should be taken to ensure that properties are
	 * uniquely named by using the caller's class name in the name
	 * of the property.
	 *  
	 * @param name The name of the property whose value is to be set/changed.
	 * @param value The value to be set for this property. This value cannot
	 * be null.
	 * @return The previous value for this property (if it had a value)
	 */
	public synchronized String  setProperty(String name, String value) {
		return (String) properties.setProperty(name, value);
	}

	/**
	 * Set/change the value for a given property name.
	 * 
	 * This method must be used to set (initially) or change the value
	 * associated with a given property. The name of properties must be
	 * unique. No checks are made to ensure property names are unique.
	 * Therefore care should be taken to ensure that properties are
	 * uniquely named by using the caller's class name in the name
	 * of the property.
	 *  
	 * @param name The name of the property whose value is to be set/changed.
	 * @param value The value to be set for this property. 
	 * @return The previous value for this property (if it had a value) or false.
	 */
	public synchronized boolean  setProperty(String name, boolean value) {
		final String strValue = (String) properties.setProperty(name, Boolean.toString(value));
		return (strValue == null) ? false : Boolean.getBoolean(strValue);
	}

	
	/**
	 * Helper method to save the properties to disk.
	 * 
	 * This is a convenient helper method that must be used to save the
	 * properties file to disk. This method saves the properties file
	 * named {@link #PEACE_PROPERTIES_FILE} in the default home directory
	 * obtained via call to {@link Utilities#getDefaultDirectory()}. 
	 * 
	 * @param parent The parent window to be used to display any error
	 * messages.
	 * 
	 * @param showErrors If this flag is true then errors encountered 
	 * while saving are displayed to the user. If this flag is false
	 * then a GUI window is not displayed. In any case, errors are
	 * always logged in the ProgrammerLog.
	 * 
	 * @return The exception (if any) that occurred during writing the
	 * data. If no exceptions occur then this method returns null.
	 */
	public synchronized Exception save(final JComponent parent, 
			final boolean showErrors) {
		Exception result = null;
		final String propFile = Utilities.getDefaultDirectory() +
			File.separator + PEACE_PROPERTIES_FILE;
		try {
			FileOutputStream fos = new FileOutputStream(propFile);
			final String comment =  "PEACE Properties (Saved at: " + new Date() + ")";
			properties.storeToXML(fos, comment, "UTF-8");
			fos.close();
			UserLog.log(Log.LogLevel.INFO, "PEACEProperties", 
					"Properties saved successfully to " + propFile);
		} catch (Exception exp) {
			UserLog.log(Log.LogLevel.WARNING, "PEACEProperties", 
			"Error saving properties to " + propFile);
			ProgrammerLog.log(exp);
			if (showErrors) {
				final JPanel msg = 
					Utilities.collapsedMessage("Unable to save properties to " + propFile,
						Utilities.toString(exp));
				JOptionPane.showMessageDialog(parent, msg, 
						"Unable to save properties", JOptionPane.ERROR_MESSAGE);
			}
			result = exp;
		}
		return result;
	}
	
	/**
	 * Helper method to load the properties from disk.
	 * 
	 * This is a convenient helper method that must be used to load the
	 * properties from disk.  Typically, this method is invoked only once
	 * when PEACE-GUI starts up.  This method loads the properties from file
	 * named {@link #PEACE_PROPERTIES_FILE} in the default home directory
	 * obtained via call to {@link Utilities#getDefaultDirectory()}. 
	 * 
	 * @param parent The parent window to be used to display any error
	 * messages.
	 * 
	 * @param showErrors If this flag is true then errors encountered 
	 * while loading are displayed to the user. If this flag is false
	 * then a GUI window is not displayed. In any case, errors are
	 * always logged in the ProgrammerLog.
	 */
	public synchronized void load(final JComponent parent, final boolean showErrors) {
		final String propFile = Utilities.getDefaultDirectory() +
			File.separator + PEACE_PROPERTIES_FILE;
		try {
			FileInputStream fis = new FileInputStream(propFile);
			properties.loadFromXML(fis);
			fis.close();
			UserLog.log(Log.LogLevel.INFO, "PEACEProperties", 
					"Properties loaded successfully from " + propFile);
			// Dump properties to the programmer log to ease troubleshooting
			// user-problems on the long run.
			ByteArrayOutputStream baos = new ByteArrayOutputStream(1024);
			PrintWriter pw             = new PrintWriter(baos);
			properties.list(pw);
			ProgrammerLog.log("Properties loaded from: " + propFile);
			ProgrammerLog.log(baos.toString());
		} catch (Exception exp) {
			UserLog.log(Log.LogLevel.WARNING, "PEACEProperties", 
			"Error loading properties from " + propFile);
			ProgrammerLog.log(exp);
			if (showErrors) {
				final JPanel msg = 
					Utilities.collapsedMessage("Unable to load properties from " + propFile,
						Utilities.toString(exp));
				JOptionPane.showMessageDialog(parent, msg, 
						"Unable to load properties", JOptionPane.ERROR_MESSAGE);
			}
		}
	}
	
	/**
	 * The constructor.
	 * 
	 * The constructor is private in order to ensure that this class is
	 * never directly instantiated. Instead the process-wide unique
	 * singleton instance of this class should be used by calling
	 * {@link #get()} method.
	 */
	private PEACEProperties() {
		properties = new Properties();
	}
	
	/**
	 * The actual set of properties begin managed by this class.
	 * 
	 *  This object contains the actual set of properties being encapsulated
	 *  by this class.
	 */
	private final Properties properties;
	
	/**
	 * The name of the file to which the properties are written.
	 * 
	 * Avoid changing this file unless you know exactly what you are
	 * doing and you have discussed this change with the team. 
	 */
	private static final String PEACE_PROPERTIES_FILE = "PEACE_properties.xml";
}
