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

package org.peace_tools.data;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.ProgressMonitor;

import org.peace_tools.generic.Pair;
import org.peace_tools.workspace.ClassifierList;
import org.peace_tools.workspace.DBClassifier;
import org.peace_tools.workspace.Workspace;

/**
 * The top-level class that encapsulates all the pertinent information
 * regarding a cluster data file. This class deserializes the information
 * in a cluster data file generated by PEACE and stores it in memory.
 * The in-memory storage format for the core cluster information is
 * achieved using an hierarchically nested set of ClusterNode objects.
 * In addition, this class also maintains any generated information
 * that is placed in the file by PEACE. <br/><br/>
 * 
 * Note that the in-memory format represented by this class has been
 * primarily designed to provide more convenient access to the related
 * information and for display in a GUI. However, this class does
 * not directly perform any GUI related task. Instead, the GUI display
 * is organized using the MVC (Model-View-Controller) design pattern.
 * This class constitutes the "model" as in the MVC terminology.
 * 
 * @note In order to create a valid ClusterFile use the loadCluster() 
 * static method in this class.
 */
public class ClusterFile {
	/**
	 * The absolute path to the file name from where the cluster data was
	 * originally loaded. 
	 * 
	 * @return The absolute path to the file that uniquely identifies
	 * the contents of the cluster.
	 */
	public String getFileName() { return this.fileName; }
	
	/**
	 * Method to print the cluster in a simple text-based format.
	 * 
	 * This method is primarily used for validating the cluster data to
	 * ensure that the data was loaded correctly.
	 * 
	 * @param out The output stream to which the cluster data is to be 
	 * written.
	 */
	public void print(PrintStream out) {
		// First print the meta data out.
		for(Pair entry : metadata) {
			out.println(entry.getName() +
					(entry.getValue() != null ? ":" + entry.getValue() : ""));			
		}
		// Now recursively print the cluster using helper method.
		root.print(out, root, "");
	}
	
	/**
	 * Obtain the top-level root cluster for this cluster file.
	 * 
	 * @return This method must return the top-level root cluster for this
	 * cluster file.
	 */
	public ClusterNode getRoot() { return root; }
	
	/**
	 * Obtain the meta data associated with this cluster file.
	 * 
	 * @return The meta data associated with this cluster file.
	 */
	public ArrayList<Pair> getMetadata() { return metadata; }
	
	/**
	 * This method loads cluster data into an in-memory format.
	 * 
	 * This method must be used to load cluster data from a PEACE generated
	 * cluster data file and deserialize the information into the in-memory
	 * format. The in-memory format provides an hierarchical organization
	 * that is a bit more streamlined and easier to display in the GUI.
	 * 
	 * @param clusterFile The cluster file (generated by PEACE) from where the
	 * data is to be loaded in the in-memory format.
	 * 
	 * @return On success this method returns a valid cluster data structure
	 * loaded from the file.
	 * 
	 * @exception Exception This method throws an exception on errors.
	 */
	public static ClusterFile loadCluster(File clusterFile) throws Exception {
		FileInputStream fis = new FileInputStream(clusterFile);
		return loadCluster(clusterFile.getAbsolutePath(), fis);
	}
	
	/**
	 * This method loads cluster data into an in-memory format.
	 * 
	 * This method must be used to load cluster data from a PEACE generated
	 * cluster data file and deserialize the information into the in-memory
	 * format. The in-memory format provides an hierarchical organization
	 * that is a bit more streamlined and easier to display in the GUI.
	 * 
	 * @param fileName The absolute path to the file from where the data
	 * is being read.
	 * 
	 * @param is The input stream from where the data is to be read.
	 * 
	 * @return On success this method returns a valid cluster data structure
	 * loaded from the file.
	 * 
	 * @exception Exception This method throws an exception on errors.
	 */
	public static ClusterFile loadCluster(String fileName, InputStream is) throws Exception {
		// Wrap the cluster file into a scanner to make line-by-line 
		// processing easier.
		Scanner input = new Scanner(is);
		// Create a temporary cluster object to hold the data
		ClusterFile cluster = new ClusterFile(fileName);
		// Use a array list to maintain clusters to look up index
		ArrayList<ClusterNode> nodeList = new ArrayList<ClusterNode>();
		// Process line-by-line from the cluster file.
		while (input.hasNextLine()) {
			// Read the line in and ignore blank lines
			String line = input.nextLine().trim();
			if (line.length() < 1) {
				// Ignore blank lines.
				continue;
			}
			// The line is a meta data line if it starts with a '#'
			if (line.charAt(0) == '#') {
				// Process the meta data line.
				cluster.metadata.add(makeMetadataEntry(line));
			} else {
				// Process the line with cluster node information.
				makeClusterNode(nodeList, line);
			}
		}
		// OK, now that we have read all the data update the cluster
		cluster.root = nodeList.get(0);
		// return the newly loaded/created cluster.
		return cluster;
	}
	
	/**
	 * Determine if the clusters are current with work space classifiers.
	 * 
	 * This method can be used to determine if the clusters have already
	 * been classified using the current set of classifiers configured for
	 * this work space. 
	 * 
	 * @return This method returns true if the clusters have already been
	 * classified using the current set of classifiers.
	 */
	public synchronized boolean isClassified() {
		ClassifierList clasList = Workspace.get().getClassifierList();
		ArrayList<DBClassifier> classifiers = clasList.getClassifiers();
		// Intentional use of == to check references!
		return (classifiers == prevClusterList);
	}
	
	/**
	 * Method to recompute (as needed) the classification of ESTs in clusters.
	 * 
	 * This method must be used to recompute the classification of ESTs in a
	 * cluster whenever the data base classifiers in the EST change. This 
	 * method recomputes classifications only if the classifier list in the 
	 * work space has actually changed. Consequently, repeatedly calling this
	 * method does not have side effects. However, if classification is performed
	 * then it may be a long running process particularly, for large EST sets.
	 * Consequently, it is best to invoke this method from a separate thread
	 * so that the GUI does not appear to be hanging while classifications are
	 * computed.
	 * 
	 * @note This method assumes that the classifications for the ESTs
	 * in the estList have already been computed.
	 * 
	 * @param estList The list of ESTs associated with this cluster file to
	 * be used to classify the ESTs.
	 * 
	 * @param pm An optional progress monitor to be updated to indicate
	 * progress. This parameter can be null.
	 */
	public synchronized void classify(final ESTList estList, ProgressMonitor pm) {
		ClassifierList clasList = Workspace.get().getClassifierList();
		ArrayList<DBClassifier> classifiers = clasList.getClassifiers();
		if (this.prevClusterList != classifiers) {
			root.classify(estList, pm);
			this.prevClusterList = classifiers;
		}
	}
	
	/**
	 * Helper method to process a comma separated set of values
	 * representing a cluster node.
	 * 
	 *  This method is a helper method that is used to parse a data
	 *  line from the cluster file and convert it to a ClusterNode. This
	 *  method reads and validates the data. It then creates a clusterNode
	 *  and adds it to its parent (if one is present) and to the nodeList.
	 * 
	 * @param nodeList The list of nodes that have been read so far.
	 * 
	 * @param line The line containing node data to be processed and
	 * converted to a clusterNode.
	 * 
	 * @throws IOException This method throws an exception if the data
	 * was invalid or not read.
	 */
	protected static void makeClusterNode(ArrayList<ClusterNode> nodeList, 
			String line) throws IOException {
		// Extract the ',' separated values. The values are in the form:
		// C, <clstrID>, <parentClstrID>   OR 
		// E, <estIdx, clstrId>
		String[] entries = line.split(",");
		// Ensure we have exactly 3 entries.
		if (entries.length != 3) {
			throw new IOException("Invalid data line in cluster file " +
					"(line: " + line + ")");
		}
		// Get the clsterID or estIdx.
		int nodeID = Integer.parseInt(entries[1]);
		ClusterNode node = null;
		switch("CE".indexOf(entries[0].charAt(0))) {
		case 0: // Cluster node
			node = new ClusterNode(null, false, nodeID);
			// Add cluster node to the list of nodes for future reference
			nodeList.add(node);
			break;
		case 1: // EST node
			node = new ClusterNode(null, true, nodeID);
			break;
		default: // Bad node!
			throw new IOException("Unknown data line in cluster file " +
					"(line: " + line + ")");
		}
		// Add the node to its parent, if we have a valid parent.
		int parentID = Integer.parseInt(entries[2]);
		if (parentID != -1) {
			if ((parentID < 0) || (parentID >= nodeList.size())) {
				// This is an invalid parent reference!
				throw new IOException("Invalid parent ID in cluster file " +
						"(line: " + line + ")");
			}
			nodeList.get(parentID).addChild(node);
		}
	}
	
	/**
	 * This is a helper method that is used to parse a line of
	 * meta data entry (line starts with a '#' character) and convert 
	 * it to a a Pair containing a name, value pair and returns the 
	 * meta data as a pair. It assumes that the name and value are
	 * separated by ':'. 
	 * 
	 * @param line The line from the cluster file to be processed.
	 * 
	 * @return A pair object containing the name, value pair.
	 */
	protected static Pair makeMetadataEntry(String line) {
		// Locate the ":" in the line and use that as the separator
		// between name:value.
		int colonPos = line.indexOf(':');
		String name  = line;
		String value = null;
		// Extract value only if ':' was found.
		if (colonPos != -1) {
			// Extract name and value.
			name  = line.substring(1, colonPos).trim();
			value = line.substring(colonPos + 1).trim();
		}
		// Create and store the information
		return new Pair(name, value);
	}
	
	/**
	 * The constructor creates an empty cluster object. This method is
	 * called from the loadCluster method. This object is filled in
	 * later on as data is read from the cluster file. 
	 * 
	 * @param fileName The absolute path to the file from where the
	 * cluster data was loaded. This file name is used as an identifier
	 * to locate the files.
	 */
	private ClusterFile(String fileName) {
		this.fileName        = fileName;
		this.root            = null;
		this.metadata        = new ArrayList<Pair>();
		this.prevClusterList = null;
	}
	
	/**
	 * The file name from where the data has been read. This is typically
	 * the absolute path of the file from where the data was read.
	 */
	private String fileName;
	
	/**
	 * The root of the cluster node. Typically there is at least one root
	 * node.
	 */
	private ClusterNode root;
	
	/**
	 * The set of meta data that was loaded from the cluster file. The
	 * meta data is stored as a list of name, value pairs.
	 */
	private ArrayList<Pair> metadata;
	
	/**
	 * This reference is used to track the object used to hold the list
	 * of DBClassifiers in the work space. This information is used to
	 * decide if the classifications need to be recomputed. This works
	 * because the classifiers are modified as a single complete batch
	 * of entries. Each time a new object is set and if the objects have
	 * changed then the classifications need to be recomputed for this
	 * cluster file.
	 */
	private Object prevClusterList;
}