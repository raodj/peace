package org.peace_tools.data;

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

import java.awt.Component;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.WeakHashMap;

import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.ProgressMonitorInputStream;

import org.peace_tools.generic.Utilities;

/**
 * A helper class to load and maintain data files.
 * 
 * This class is a helper class that is used to load and maintain
 * information about ESTs, MSTs, and Cluster data files. This class
 * maintains references to currently loaded entries so that 
 * various data files can be shared and reused. The information is
 * maintained in a weak hash map. Therefore, if an entry is
 * not used (at least for some time) then it will be automatically
 * garbage collected and needs to be reloaded.
 *
 * @note This class is not meant to be directly instantiated to
 * preserve singleton property. Use the static get() method to
 * obtain a reference to the globally unique instance of this class.
 */
public class DataStore {
	/**
	 * Obtain reference to the globally unique singleton instance.
	 * 
	 * @return The globally unique singleton instance of the data store.
	 */
	public static DataStore get() { return dataStore; }
	
	/**
	 * Helper method to ensure there is sufficient memory to load data.
	 * 
	 * This is a helper method that is called just before data files
	 * are loaded to ensure that sufficient memory is available to 
	 * load large files. If not, an warning message is displayed to
	 * the user along with information on how to increase memory.
	 * 
	 * @param fileToBeLoaded The file that is going to be loaded into
	 * memory.
	 * 
	 * @param parent The parent component based on which any warning
	 * dialog boxes are to be displayed.
	 */
	public void memoryCheck(File fileToBeLoaded, Component parent) throws Exception {
		Runtime rt = Runtime.getRuntime();
		float maxMemory      = rt.maxMemory();
		float reservedMemory = rt.totalMemory();
		float availMemory    = maxMemory - reservedMemory;
		
		if (availMemory >= fileToBeLoaded.length() * 2.0f) {
			// Sufficient memory is available (we think).
			return;
		}
		final float TO_MEGS = 1024 * 1024;
		// Warn the user that sufficient memory is not available.
		final String WarnMsg = String.format(LOW_MEMORY_MSG, 
				fileToBeLoaded.getName(), fileToBeLoaded.length() / TO_MEGS,
				fileToBeLoaded.length() * 2.0f / TO_MEGS);
		final String MemUse = String.format(MEM_USAGE, 
				availMemory / TO_MEGS, reservedMemory / TO_MEGS,
				maxMemory / TO_MEGS);
		// Create a suitable collapsed pane to display message
		JPanel message = Utilities.collapsedMessage(WarnMsg, MemUse);
		int choice = JOptionPane.showConfirmDialog(parent, message,
				"Low memory warning", JOptionPane.YES_NO_OPTION, 
				JOptionPane.WARNING_MESSAGE);
		if (choice == JOptionPane.NO_OPTION) {
			throw new Exception("User decided to stop file load due to " +
					"low memory situation.");
		}
		// Try to free up some memory if possible.
		System.gc();
	}
	
	/**
	 * Method to load/get a given FASTA file.
	 * 
	 * This method can be used to obtain the ESTs loaded from a given
	 * FASTA file. The full path to the FASTA file must be supplied
	 * as the parameter to this method. If the specified file is not
	 * found in cache, then this method loads the file from disk,
	 * places the entry in cache, and returns a reference to the 
	 * ESTList.
	 * 
	 * @param fileName The FASTA file containing ESTs to be loaded.
	 * 
	 * @param parent If this parameter is not null then this method
	 * displays a progress monitor.
	 * 
	 * @return This method returns the list of ESTs from the given
	 * file.
	 * 
	 * @exception Exception This method throws exceptions if the 
	 * specified file could not be loaded due to various reasons.
	 */
	public ESTList getFASTA(String fileName, Component parent) throws Exception {
		// First check if entry is in cache.
		Object entry = cache.get(fileName);
		// If entry is present then check and return it. 
		if (entry != null) {
			if (entry instanceof ESTList) {
				return (ESTList) entry;
			} else {
				throw new IllegalArgumentException("The specified file did not map to a FASTA file.");
			}
		}
		// The entry was not found. It must be loaded. So do it now.
		File file = new File(fileName);
		if (!file.exists()) {
			throw new IOException("The file " + fileName + " was not found.");
		}
		// Do a memory check to ensure we have sufficient memory
		memoryCheck(file, parent);
		InputStream fis = new FileInputStream(file);
		// Wrap the input stream in progress monitor if requested.
		if (parent != null) {
			fis = new ProgressMonitorInputStream(parent, 
					"Please wait while the file is loaded...\n" +
					"(file: " + fileName + ")", fis);
		}
		// Load the ESTs
		ESTList ests = ESTList.loadESTs(fileName, fis);
		// Put entry in cache for future reference.
		cache.put(fileName, ests);
		// Now the caller gets a reference to hold on to.
		return ests;
	}
	
	/**
	 * Method to load/get a given cluster file.
	 * 
	 * This method can be used to obtain the cluster data loaded from 
	 * a given cluster file. The full path to the data file must be 
	 * supplied as the parameter to this method. If the specified file
	 * is not found in cache, then this method loads the file from 
	 * disk, places the entry in cache, and returns a reference to the 
	 * cluster information.
	 * 
	 * @param fileName The cluster file containing the data to be loaded.
	 * 
	 * @param parent If this parameter is not null then this method
	 * displays a progress monitor.
	 * 
	 * @return This method returns the cluster data  from the given
	 * file.
	 * 
	 * @exception Exception This method throws exceptions if the 
	 * specified file could not be loaded due to various reasons.
	 */
	public ClusterFile getClusterData(String fileName, Component parent) throws Exception {
		// First check if entry is in cache.
		Object entry = cache.get(fileName);
		// If entry is present then check and return it. 
		if (entry != null) {
			if (entry instanceof ClusterFile) {
				return (ClusterFile) entry;
			} else {
				throw new IllegalArgumentException("The specified file did not map to a FASTA file.");
			}
		}
		// The entry was not found. It must be loaded. So do it now.
		InputStream fis = new FileInputStream(fileName);
		// Wrap the input stream in progress monitor if requested.
		if (parent != null) {
			fis = new ProgressMonitorInputStream(parent, 
					"Please wait while the file is loaded...\n" +
					"(file: " + fileName + ")", fis);
		}
		// Load the ESTs
		ClusterFile clusters = ClusterFile.loadCluster(fileName, fis);
		// Put entry in cache for future reference.
		cache.put(fileName, clusters);
		// Now the caller gets a reference to hold on to.
		return clusters;
	}

	/**
	 * Method to load/get a given cluster file.
	 * 
	 * This method can be used to obtain the cluster data loaded from 
	 * a given cluster file. The full path to the data file must be 
	 * supplied as the parameter to this method. If the specified file
	 * is not found in cache, then this method loads the file from 
	 * disk, places the entry in cache, and returns a reference to the 
	 * cluster information.
	 * 
	 * @param fileName The cluster file containing the data to be loaded.
	 * 
	 * @param parent If this parameter is not null then this method
	 * displays a progress monitor.
	 * 
	 * @return This method returns the cluster data  from the given
	 * file.
	 * 
	 * @exception Exception This method throws exceptions if the 
	 * specified file could not be loaded due to various reasons.
	 */
	public MST getMSTData(String fileName, Component parent) throws Exception {
		// First check if entry is in cache.
		Object entry = cache.get(fileName);
		// If entry is present then check and return it. 
		if (entry != null) {
			if (entry instanceof MST) {
				return (MST) entry;
			} else {
				throw new IllegalArgumentException("The specified file did not map to a FASTA file.");
			}
		}
		// The entry was not found. It must be loaded. So do it now.
		InputStream fis = new FileInputStream(fileName);
		// Wrap the input stream in progress monitor if requested.
		if (parent != null) {
			fis = new ProgressMonitorInputStream(parent, 
					"Please wait while the file is loaded...\n" +
					"(file: " + fileName + ")", fis);
		}
		// Load the ESTs
		MST mst = MST.loadMST(fileName, fis);
		// Put entry in cache for future reference.
		cache.put(fileName, mst);
		// Now the caller gets a reference to hold on to.
		return mst;
	}

	/**
	 * The constructor. This class is not meant to be instantiated.
	 * Therefore the constructor is private. Use the get() method
	 * to obtain a reference to the global instance of this class.
	 */
	private DataStore() {
		// Initialize the cache.
		cache = new WeakHashMap<String, Object>();
	}
	
	/**
	 * A weak hash map that serves as a cache to hold MST, EST, 
	 * and cluster data file entries in it. Entries are added
	 * to this cache on demand, whenever new files need to be 
	 * loaded.
	 */
	private final WeakHashMap<String, Object> cache;
	
	/**
	 * Message to be displayed to the user to indicate that
	 * there is a low memory situation when loading a file.
	 * This message is formatted (to fill in missing
	 * information) by the memoryCheck() method.
	 */
	private static final String LOW_MEMORY_MSG = "<html>" +
		"Your Java VM is running low on memory and may not be able to open the<br>" +
		"file: %s (size=%.1f MB). Memory safety threshold: %.1f MB<br>" +
		"See details below for current memory usage statistics.<br>" +
		"You can proceed with opening the file (<i>but may experience problems</i>)<br>" +
		"or try closing open files to freeup some memory and then open this file.<br><br>" +
		"The best solution would be to exit PEACE and restart it with a larger heap<br>" +
		"as shown in the command line below (set suitable heap value instead of 2G):<br>" +
		"<b>java -Xmx2G -jar peace.jar<br>" +
		"</html>";

	/**
	 * Memory statistics to be displayed to the user after 
	 * formatting (to fill in the necessary information). This
	 * string is used in the memoryCheck() method.
	 */
	private static final String MEM_USAGE =  
		"Free memory available: %.1f MB\n" +
		"Current memory used: %.1f MB\n" +
		"Maximum memory available to Java: %.1f MB";

	/**
	 * The globally unique singleton instance of this class that is
	 * shared by all the other classes that require data files to
	 * be loaded.
	 */
	private static final DataStore dataStore = new DataStore();
}
