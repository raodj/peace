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
import org.w3c.dom.Element;

/**
 * <p>
 * This class encapsulates all the data necessary to use a given clustering
 * information generated from a Spanning Tree (MST) data file. This element
 * contains meta data about the MST file. The MSTClusterData class is created
 * each time a new Job is run to cluster based on a given MST data set. The MST
 * data is the primary input from which the cluster data is derived.
 * </p>
 * 
 * <p>
 * In addition to encapsulating the data this class also provides convenient
 * interfaces for marshalling and unmarshalling XML data compatible with the
 * PEACE GUI configuration XML file format.
 * </p>
 */
public class MSTClusterData {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * MSTClusterData entry. This method is typically  used to create a 
	 * suitable entry when loading a Work space into the GUI.
	 * 
	 * @param clstrData The DOM element to be used for creating the cluster
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created data entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static MSTClusterData create(Element clstrData) throws Exception {
		// First extract the necessary information from the DOM tree.
		String id      = DOMHelper.getStringValue(clstrData, "ID");
		String mstID   = DOMHelper.getStringValue(clstrData, "MSTRef");
		String path    = DOMHelper.getStringValue(clstrData, "Path");
		String desc    = DOMHelper.getStringValue(clstrData, "Description", true);
		int    thresh  = DOMHelper.getIntValue(clstrData, "Threshold");
		// load up the job summary information.
		Element jobData    = DOMHelper.getElement(clstrData, "JobSummary");
		JobSummary summary = JobSummary.create(jobData);
		// Now that we have all the information create the actual 
		// cluster data node.
		return new MSTClusterData(id, mstID, path, desc, thresh, summary);
	}
	
	/**
	 * The constructor to create a fully populated MSTClusterData object that 
	 * contains all the meta data regarding a MST file that has been generated
	 * from a given EST file.
	 * 
	 * @param id The workspace wide unique ID associated with this
	 * entry.
	 * 
	 * @param mstID
	 *            The unique MST data set ID value for this data.
	 * @param path
	 *            The path to the actual cluster file (on the local machine) that is
	 *            referred by this entry.
	 * @param description
	 *            A user defined description for this file.
	 * @param threshold
	 * 			  The threshold value used to partition a MST.
	 * @param summary
	 *            The core/useful information about the job that was run to
	 *            compute the MST.
	 */
	public MSTClusterData(String id, String mstID, String path, 
			String description, int threshold, JobSummary summary) {
		this.id          = id;
		this.mstID       = mstID;
		this.path        = path;
		this.description = description;
		this.threshold   = threshold;
		this.jobSummary  = summary;
	}
	
	/**
	 * Obtain the work space wide unique identifier set for 
	 * this entry. 
	 * 
	 * The ID value is created when a job is scheduled to compute
	 * clustering from a data set. The ID is persisted in the work
	 * space configuration file and loaded when a work space is 
	 * opened in the GUI.
	 * 
	 * @return This method returns the unique identifier set for this
	 * data set.
	 */
	public String getID() { return id; }
	
	/**
	 * Obtain the work space wide unique identifier set for this MST data set.
	 * The ID value is created when new jobs are scheduled to create cluster
	 * data. These IDs are persisted in the work space configuration file and
	 * loaded when a work space is opened in the GUI.
	 * 
	 * @return This method returns the unique MSTData identifier set for 
	 * this data set.
	 */
	public String getMSTID() { return mstID; }
		
	/**
	 * Obtain the complete file name and path where the actual cluster file is
	 * located. This value is set when a new job to generate this cluster data
	 * is scheduled. The path for the file is persisted in the work space
	 * configuration file and loaded when a work space is opened in the GUI.
	 * 
	 * @return This method returns the path to the cluster file.
	 */
	public String getPath() { return path; }
	
	/**
	 * Obtain the user supplied description for this this data. The value is set
	 * when new jobs are scheduled. The description is persisted in the work
	 * space configuration file and loaded when a work space is opened in the
	 * GUI.
	 * 
	 * @return This method returns the user supplied description set for this
	 *         job.
	 */
	public String getDescription() { return description; }

	/**
	 * Obtain meta data about the job that was run to generate the cluster. The
	 * JobID in the job summary can be used to look up additional information
	 * about the job. Similarly, the serverID value in the job summary can be
	 * used to look up information about the server on which the job was run.
	 * 
	 * @return This method returns information about the job that was run to
	 *         generate the cluster data.
	 */
	public JobSummary getJobSummary() { return jobSummary; }

	/**
	 * The threshold value used for cluster generating. The threshold defines
	 * the metric on a given edge of a MST that identifies transition to a 
	 * new cluster.
	 * 
	 * @return The threshold value that was used to generate the clusters.
	 */
	public int getThreshold() { return threshold; }
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent DataSet node in the DOM tree.
	 * 
	 * @param dataset The DOM element corresponding to the "DataSet"
	 * node that contains this entry.
	 */
	public final void marshall(Element dataset) {
		// Create a top-level entry for this "Job"
		Element clstrData = DOMHelper.addElement(dataset, "MSTClusterData", null);
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(clstrData, "ID", id);
		DOMHelper.addElement(clstrData, "MSTRef", mstID);
		DOMHelper.addElement(clstrData, "Path", path);
		DOMHelper.addElement(clstrData, "Description", 
				(description != null) ? description : "");
		DOMHelper.addElement(clstrData, "Threshold", "" + threshold);
		// Finally marshall the job summary information.
		jobSummary.marshall(clstrData);
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t";
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		final String NUM_ELEMENT = Indent + "\t" + "<%1$s>%2$d</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<MSTClusterData>\n", Indent); 
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "ID", id);
		out.printf(STR_ELEMENT, "MSTRef", mstID);
		out.printf(STR_ELEMENT, "Path", path);
		out.printf(STR_ELEMENT, "Description", 
				(description != null) ? DOMHelper.xmlEncode(description) : "");
		out.printf(NUM_ELEMENT, "Threshold", threshold);
		// Marshall the job summary out.
		jobSummary.marshall(out);
		// Close the MSTData tag
		out.printf("%s</MSTClusterData>\n", Indent);
	}
	
	/**
	 * Overrides the default implementation in the base class to 
	 * simply return the last part of the cluster file associated
	 * with this entry.
	 * 
	 * @return A short string representation that is easy to display
	 * primarily in a tree view of the work space.
	 */
	@Override
	public String toString() {
		File tmpData = new File(path);
		return "Cluster Data [File: " + tmpData.getName() + "]";
	}
	
	/**
	 * Return the information in the form of a partial PEACE command line 
	 * parameter.
	 * 
	 * This method can be used to obtain the information needed to
	 * generate the MST based on the supplied information in the form
	 * of a command line parameter. 
	 * 
	 * @return Return the information as a command line parameter.
	 */
	public String toCmdLine() {
		String cmdLine = "";
		File file = new File(path);
		cmdLine += " --gui-print --output " + file.getName(); 
		return cmdLine;
	}
	
    /**
     * Set revised updated information about the job associated with
     * this MSTData entry.
     * 
     * @param job The job from where the necessary information is
     * to be copied.
     */
    public void updateJobSummary(Job job) {
    	jobSummary = new JobSummary(job);
   		// Notify all listeners about the change.
    	Workspace ws = Workspace.get();
		ws.fireWorkspaceChanged(new WorkspaceEvent(this, 
				WorkspaceEvent.Operation.INSERT));
    }
    
    /**
     * Set the data set that this object logically belongs to. 
     * 
     * <p><b>Note:</b>  This method is used by the DataSet. Use 
     * {@link DataSet#add(MSTClusterData)} method to add a 
     * {@link MSTClusterData} object to a data set.<p>
     * 
     * @param dataSet The data set to which this object has been added.
     */
    protected void setDataSet(DataSet dataSet) {
            this.dataSet = dataSet;
    }
    
    /**
     * Obtain the data set that logically contains this MST Data.
     * 
     * @return The data set that logically contains this object. This value
     * is set only after this object has been added to a data set.
     */
    public DataSet getDataSet() { return dataSet; }
    
    /**
	 * Provides a textual multi-line summary for this element.
	 * 
	 * This method is typically used by GUI to display a simple 
	 * representation of the data stored in this element.
	 * 
	 * @param indent The leading spaces to be used for indenting each
	 * line. This parameter cannot be null.
	 * 
	 * @return A textual multi-line summary of the data in this element.
	 */
	public String getSummary(String indent) {
		StringBuilder sb = new StringBuilder(128);
		sb.append(indent + "Path: " + path + "\n");
		sb.append(indent + "Description: " + description + "\n");
		sb.append(indent + "Threshold: " + threshold + "\n");
		sb.append(indent + "MST ID: " + mstID + "\n");
		return sb.toString();
	}
	
    /**
     * Reference to the data set that logically contains this MSTData 
     * element. This element is not serialized to XML format 
     * (and consequently is flagged transient).
     */
    private transient DataSet dataSet;

	/**
	 * The unique generated ID value for this entry. This value is
	 * created when a new entry is added. This value is persisted in
	 * the work space configuration. It provides a convenient mechanism
	 * to refer to a specific entry within the workspace. The IDs are 
	 * generated via  a call to Workspace.reserveID() method.
	 */
	private String id;
	
	/**
	 * The unique value for the MSTData entry based on which this cluster data
	 * was generated. This value value is persisted in the work space
	 * configuration.
	 */
	protected final String mstID;
	
	/**
	 * The path to the actual MST file (on the local machine) that is 
	 * referred by this entry.
	 */
	private final String path;

	/**
	 * A user defined description for this file. This description is set when
	 * a new job is scheduled and persisted in the work space configuration for
	 * future references.
	 */
	private final String description;
	
	/**
	 * The threshold value used for cluster generating. The threshold defines
	 * the metric on a given edge of a MST that identifies transition to a 
	 * new cluster.
	 */
	private final int threshold;
	
	/**
	 * The core/useful information about the job that was run to compute
	 * the MST.
	 */
	private JobSummary jobSummary;
}
