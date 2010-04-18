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

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * The top-level data set that encapsulates the meta data regarding data 
 * files that are all related to a single EST file. The objective of this
 * organization is to minimize data redundancy. Typically a single EST 
 * data file is analyzed multiple times in order to determine a suitable
 * clustering. This scheme (that is centered around the primary EST data
 * file) is meant to reflect the practical use of these files and work
 * with them.
 */
public class DataSet {
	/**
	 * Enumeration of the type of data file associated with this data set.
	 * 
	 * The following enumerations are used to indicate the file format of
	 * the data file associated with this data set.
	 */
	public enum DataFileType {
		/**
		 * This enumeration is used to indicate a FASTA (text) file. This
		 * is the default data type that is associated with a data file.
		 */
		FASTA,
		/**
		 * This enumeration is used to indicate a Standard Flowgram Format
		 * (SFF) file.
		 */
		SFF };
	
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * DataSet entry. This method is typically  used to create a suitable
	 * entry when loading a work space into the GUI.
	 * 
	 * @param data The DOM element to be used for creating the DataSet
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created data set object based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static DataSet create(Element data) throws Exception {
		// First extract the necessary information from the DOM tree.
		Element estData= DOMHelper.getElement(data, "ESTData");
		String id      = DOMHelper.getStringValue(estData, "ID"); 
		String path    = DOMHelper.getStringValue(estData, "Path");
		String desc    = DOMHelper.getStringValue(estData, "Description", true);
		String type    = estData.getAttribute("fileType");
		desc           = (desc != null) ? desc : "";
		// Convert file type to a suitable enumeration for further use
		DataFileType fileType= ((type != null) && (type.length() > 0)) ? 
				DataFileType.valueOf(type.toUpperCase()) : DataFileType.FASTA; 
		// Create the data set entry.
		DataSet dataSet= new DataSet(id, path, desc, fileType);
		// Now parse in any MSTData elements for this DataSet
		NodeList mstNodes = data.getElementsByTagName("MSTData");
		for(int idx = 0; (idx < mstNodes.getLength()); idx++) {
			Element node = (Element) mstNodes.item(idx);
			// Create a MST data set using the node.
			MSTData mst = MSTData.create(node);
			mst.setDataSet(dataSet);
			dataSet.mstList.add(mst);
		}
		// Now parse in any Cluster Data elements for this DataSet
		NodeList clstrNodes = data.getElementsByTagName("MSTClusterData");
		for(int idx = 0; (idx < mstNodes.getLength()); idx++) {
			Element node = (Element) clstrNodes.item(idx);
			// Create a MST data set using the node.
			MSTClusterData clstr = MSTClusterData.create(node);
			clstr.setDataSet(dataSet);
			dataSet.clusterList.add(clstr);
		}
		// All the data was parsed in successfully.
		return dataSet;
	}
	
	/**
	 * Constructor to create a Data Set entry.
	 * 
	 * This constructor must be used to create a new data set entry to
	 * be added to the work space. 
	 * 
	 * @param id The workspace-wide unique ID to be set for this data
	 * set. This ID must be obtained via a call to Workspace.reserveID().
	 * 
	 * @param path
	 *            The path to the actual EST file (on the local machine) that is
	 *            referred by this entry.
	 * @param description
	 *            A user defined description for the EST data file.
	 *            
	 * @param fileType The physical file format of the data file associated with
	 * this data set.
	 */
	public DataSet(String id, String path, String description, DataFileType fileType) {
		this.id          = id;
		this.path        = path;
		this.description = description;
		this.fileType    = fileType;
		this.mstList     = new ArrayList<MSTData>();
		this.clusterList = new ArrayList<MSTClusterData>();
	}
	
	/**
	 * Obtain the complete file name and path where the actual EST file is
	 * located. This value is set when a new data set is added to the work
	 * space. The path for the EST file is persisted in the work space
	 * configuration file and loaded when a work space is opened in the GUI.
	 * 
	 * @return This method returns the path to the EST FASTA file.
	 */
	public String getPath() { return path; }

	/**
	 * Set the complete file name and path where the actual EST file is
	 * located. This value is set when a new data set is added to the work
	 * space. The path for the EST file is persisted in the work space
	 * configuration file and loaded when a work space is opened in the GUI.
	 * 
	 * @param path The complete file name and path to the EST FASTA file.
	 */
	public void setPath(String path) { this.path = path; }

	/**
	 * Obtain the user supplied description for this this data. The value is set
	 * when a new data set is added to the work space. The description is
	 * persisted in the work space configuration file and loaded when a work
	 * space is opened in the GUI.
	 * 
	 * @return This method returns the user supplied description set for this
	 *         EST data file.
	 */
	public String getDescription() { return description; }

	/**
	 * Obtain the work space wide unique identifier set for 
	 * this data set. The ID value is created when data sets are added.
	 * The ID is persisted in the work space configuration file and 
	 * loaded when a work space is opened in the GUI.
	 * 
	 * @return This method returns the unique identifier set for this
	 * data set.
	 */
	public String getID() { return id; }
	
	/**
	 * Set the workspace-wide ID for this DataSet.
	 * 
	 * This method can be used to reset the ID associated with this data 
	 * set. This method can be used as long as the data set has not been
	 * added to the workspace. After it has been added, it is unwise to 
	 * change the ID.
	 * 
	 * @param id The new ID to be set for this entry. This value is obtained
	 * via a call to Workspace.reserveID() method.
	 */
	public void setID(String id) {
		this.id = id;
	}
	
	/**
	 * Set an user supplied description for this this data. The value is set
	 * when a new data set is added to the work space. The description is
	 * persisted in the work space configuration file and loaded when a work
	 * space is opened in the GUI.
	 * 
	 * @param desc The user supplied description set for this EST file.
	 */
	public void setDescription(String desc) { description = desc; }

	/**
	 * Obtain information about the physical file format of the data file.
	 * 
	 * This method must be used to determine the file format of the physical
	 * data file.
	 * 
	 * @return One of the pre-defined (and supported) file formats.
	 */
	public DataFileType getFileType() { return fileType; }
	
	/**
	 * Set the information about the physical file format of the data file.
	 * 
	 * This method must be used to set the format of the physical data file
	 * associated with this data set.
	 * 
	 * @param fileType A pre-defined enumerated value indicating the physical
	 * format of the data file associated with this data set.
	 */
	public void setFileType(DataFileType fileType) {
		this.fileType = fileType;
	}
	
	/**
	 * Convenience method to detect if the file type is a FASTA file.
	 * 
	 * This method must provides a convenient mechanism to determine if the
	 * data file associated with this data set is a FASTA file. 
	 * 
	 * @return This method returns true if the data file is in FASTA file 
	 * format. Otherwise it returns false.
	 */
	public boolean isFASTAFile() {
		return fileType.equals(DataFileType.FASTA);
	}
	
	/**
	 * Obtain the list of MST data files encapsulated by this data set.
	 * 
	 * @return The list (zero or more) of MST data files that are associated
	 * with this data set.
	 */
	public ArrayList<MSTData> getMSTList() { return mstList; }
	
	/**
	 * Obtain the list of clustering data files encapsulated by this data set.
	 * 
	 * @return The list (zero or more) of clustering data files that are
	 *         associated with this data set.
	 */
	public ArrayList<MSTClusterData> getClusterList() { return clusterList; }

	/**
	 * Obtain the MSTData entry for a given job ID.
	 * 
	 * @param jobID The ID of the job for which the MST data is to be
	 * searched and retrieved.
	 * 
	 * @return The MSTData corresponding to the given job ID. If the
	 * entry was not found this method returns null.
	 */
	public MSTData getMSTData(String jobID) {
		for(MSTData entry: mstList) {
			if (jobID.equals(entry.getJobSummary().getJobID())) {
				return entry;
			}
		}
		// Entry not found.
		return null;
	}
	
	/**
	 * Obtain the MSTClusterData entry for a given job ID.
	 * 
	 * @param jobID The ID of the job for which the cluster entry
	 * is to be searched and retrieved.
	 * 
	 * @return The MSTClusterData corresponding to the given job ID.
	 * If the entry was not found this method returns null.
	 */
	public MSTClusterData getClusterData(String jobID) {
		for(MSTClusterData entry: clusterList) {
			if (jobID.equals(entry.getJobSummary().getJobID())) {
				return entry;
			}
		}
		// Entry not found.
		return null;
	}

	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent Workspace node in the DOM tree.
	 * 
	 * @param workspace The DOM element corresponding to the "Workspace"
	 * node that contains this entry.
	 */
	public final void marshall(Element workspace) {
		// Create a top-level entry for this "DataSet"
		Element dataset = DOMHelper.addElement(workspace, "DataSet", null);
		// Add new sub-element for the ESTData node
		Element estData = DOMHelper.addElement(dataset, "ESTData", null);
		estData.setAttribute("fileType", fileType.toString().toLowerCase());
		DOMHelper.addElement(estData, "ID", id);
		DOMHelper.addElement(estData, "Path", path);
		DOMHelper.addElement(estData, "Description", description);
		// Add new sub-elements for each MSTData entries.
		for(MSTData mst : mstList) {
			mst.marshall(dataset);
		}
		// Add new sub-elements for each MSTClusterData entries.
		for(MSTClusterData clstr : clusterList) {
			clstr.marshall(dataset);
		}
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a XML
	 * fragment. The XML fragment is guaranteed to be compatible with the PEACE
	 * work space configuration data.
	 * 
	 * @param out
	 *            The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t";
		final String STR_ELEMENT = Indent + "\t\t" + "<%1$s>%2$s</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<DataSet>\n", Indent); 
		// Add new sub-elements for the ESTData element
		out.printf("%s\t<ESTData fileType=\"%s\">\n", Indent, fileType.toString().toLowerCase());
		out.printf(STR_ELEMENT, "ID", id);
		out.printf(STR_ELEMENT, "Path", path);
		out.printf(STR_ELEMENT, "Description", DOMHelper.xmlEncode(description));
		out.printf("%s\t</ESTData>\n", Indent); 
		// Add new sub-elements for each MSTData entries.
		for(MSTData mst : mstList) {
			mst.marshall(out);
		}
		// Add new sub-elements for each MSTClusterData entries.
		for(MSTClusterData clstr : clusterList) {
			clstr.marshall(out);
		}
		// Close the DataSet tag
		out.printf("%s</DataSet>\n", Indent);
	}
	
	/**
	 * Add a new Minimum Spanning Tree (MST) data file to this
	 * data set. The MST data must have been generated for the EST
	 * file associated with this data set.
	 *  
	 * <p><b>Note:</b> This method reports the newly added MSTData entry to
	 * all workspace listeners by firing a suitable event.</p>
	 * 
	 * @param mstData The new MSTData entry to be added to this 
	 * data set.
	 */
	public synchronized void add(MSTData mstData) {
		Workspace ws = Workspace.get();
		if (mstData != null) {
			mstList.add(mstData);
			mstData.setDataSet(this);
			WorkspaceEvent wse = new WorkspaceEvent(mstData,
					WorkspaceEvent.Operation.INSERT);
			ws.fireWorkspaceChanged(wse);
		}
	}
	
	/**
	 * Remove an existing Minimum Spanning Tree (MST) data file from 
	 * this data set.
	 *  
	 * <p><b>Note:</b> This method reports the removed MSTData entry to
	 * all workspace listeners by firing a suitable event.</p>
	 * 
	 * @param mstData The MSTData entry to be removed from this data set. 
	 */
	public synchronized void remove(MSTData mstData) {
		Workspace ws = Workspace.get();
		if (mstData != null) {
			mstList.remove(mstData);
			WorkspaceEvent wse = new WorkspaceEvent(mstData,
					WorkspaceEvent.Operation.DELETE);
			ws.fireWorkspaceChanged(wse);
		}
	}

	/**
	 * Add a new MST-based clustering data file to this
	 * data set. The MST data must have been generated for the EST
	 * file associated with this data set.
	 *  
	 * <p><b>Note:</b> This method reports the newly added entry to
	 * all workspace listeners by firing a suitable event.</p>
	 * 
	 * @param cluster The new clustering entry to be added to this 
	 * data set.
	 */
	public synchronized void add(MSTClusterData cluster) {
		Workspace ws = Workspace.get();
		if (cluster != null) {
			clusterList.add(cluster);
			cluster.setDataSet(this);
			WorkspaceEvent wse = new WorkspaceEvent(cluster,
					WorkspaceEvent.Operation.INSERT);
			ws.fireWorkspaceChanged(wse);
		}
	}

	/**
	 * Remove an existing MST-based clustering data file from 
	 * this data set.
	 *  
	 * <p><b>Note:</b> This method reports the removed cluster entry 
	 * to all work space listeners by firing a suitable event.</p>
	 * 
	 * @param cluster The cluster entry to be removed from this 
	 * data set. 
	 */
	public synchronized void remove(MSTClusterData cluster) {
		Workspace ws = Workspace.get();
		if (cluster != null) {
			clusterList.remove(cluster);
			WorkspaceEvent wse = new WorkspaceEvent(cluster,
					WorkspaceEvent.Operation.DELETE);
			ws.fireWorkspaceChanged(wse);
		}
	}
	
	/**
	 * Overrides the default implementation in the base class to 
	 * simply return the last part of the EST file associated
	 * with this data set.
	 * 
	 * @return A short string representation that is easy to display
	 * primarily in a tree view of the work space.
	 */
	@Override
	public String toString() {
		File tmpData = new File(path);
		return "Data Set [" + fileType + " file: " + tmpData.getName() + "]";
	}
	
	/**
	 * The optional list of MST data files that are currently associated
	 * with this data set. The mst list contains the MST data that was
	 * derived by analyzing the EST file associated with this data set.
	 */
	private ArrayList<MSTData> mstList;
	
	/**
	 * The optional list of clustering data files that are currently associated
	 * with this data set. The list contains the clustering data that was
	 * derived by analyzing the MST file(s) associated with this data set.
	 */
	private ArrayList<MSTClusterData> clusterList;
	
	/**
	 * The complete path to the actual EST file on the local machine. The
	 * EST file is expected to be in FASTA file format. 
	 */
	private String path;

	/**
	 * A user defined description for this EST data file. This description 
	 * is set when a new job is scheduled and persisted in the work space 
	 * configuration for future references.
	 */
	private String description;
	
	/**
	 * The unique generated data set ID value for this data. This value is
	 * created when a new data set is added. This value is persisted in
	 * the work space configuration. It provides a convenient mechanism
	 * to refer to a specific entry within the workspace. The IDs are 
	 * generated via  a call to Workspace.reserveID() method.
	 */
	private String id;

	/**
	 * Enumeration to indicate the physical file format of the data file
	 * associated with this data set. This value is persisted in the
	 * work space configuration. It enables loading the necessary
	 * sequences for viewing and analysis.
	 */
	private DataFileType fileType;
}
