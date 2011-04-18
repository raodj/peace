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

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.generic.Utilities;
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
		 * (SFF) file. Further details on the SFF file format can be found at:
		 * http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats#sff
		 */
		SFF,
		/**
		 * This enumeration is used to indicate an ACE file format used in
		 * genomics. See http://en.wikipedia.org/wiki/ACE_file_format.
		 */
		ACE,
		/**
		 * This enumeration indicates a SAM (Sequence Alignment/Map) text file 
		 * format that is used for storing large nucleotide sequence alignments. 
		 * See http://samtools.sourceforge.net/ for details.
		 */
		SAM,
		/**
		 * This is the binary version of SAM. 
		 */
		BAM,
		/**
		 * A general Tab Separated Value (TSV) text file. 
		 */
		TSV,
		/**
		 * A regular ASCII text file.
		 */
		TXT
	};
	
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
		// Parse in any statistics element about this data set
		Element estStats = DOMHelper.getElement(data, "Stats");
		if (estStats != null) {
			dataSet.stats = DataFileStats.create(path, estStats); 
		}
		// Now parse in any MSTData elements for this DataSet
		NodeList mstNodes = data.getElementsByTagName("GeneratedFileList");
		for(int idx = 0; (idx < mstNodes.getLength()); idx++) {
			Element node = (Element) mstNodes.item(idx);
			// Create a generated file list using the node.
			GeneratedFileList gfl = GeneratedFileList.create(node);
			dataSet.gflList.add(gfl);
			gfl.setDataSet(dataSet);
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
		this.gflList     = new ArrayList<GeneratedFileList>();
		this.stats       = null;
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
	 * Obtain the list of GeneratedFileList (GFL) data files 
	 * encapsulated by this data set.
	 * 
	 * @return The list (zero or more) of GeneratedFileList objects
	 * data files that are associated with this data set. This
	 * return value is never null (but can be an empty list).
	 */
	public ArrayList<GeneratedFileList> getGflList() { return gflList; }

	/**
	 * Obtain the Generated File List (GFL) entry for a given job ID.
	 * 
	 * @param jobID The ID of the job for which the MST data is to be
	 * searched and retrieved.
	 * 
	 * @return The GFL object corresponding to the given job ID. If the
	 * entry was not found this method returns null.
	 */
	public GeneratedFileList getGFL(String jobID) {
		for(GeneratedFileList gfl: gflList) {
			if (jobID.equals(gfl.getJobSummary().getJobID())) {
				return gfl;
			}
		}
		// Entry not found.
		return null;
	}
	

	/**
	 * Method to marshal the data stored in this object to become part of
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
		// Add stats element if we have one.
		if (stats != null) {
			stats.marshall(estData);
		}
		// Add new sub-elements for each GFL entry.
		for(GeneratedFileList gfl: gflList) {
			gfl.marshall(dataset);
		}
	}
	
	/**
	 * Method to marshal the data stored in this object directly to a XML
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
		// Marshal out any statistics object we may have
		if (stats != null) {
			stats.marshall(out);
		}
		out.printf("%s\t</ESTData>\n", Indent); 
		// Add new sub-elements for each GeneratedFileList entries.
		for(GeneratedFileList gfl : gflList) {
			gfl.marshall(out);
		}
		// Close the DataSet tag
		out.printf("%s</DataSet>\n", Indent);
	}
	
	/**
	 * Add a new GeneratedFileList (GFL) object to this
	 * data set. The GeneratedFileList data must have been generated
	 * for the EST file associated with this data set.
	 *  
	 * <p><b>Note:</b> This method reports the newly added entry to
	 * all workspace listeners by firing a suitable event.</p>
	 * 
	 * @param gfl The new GeneratedFileList entry to be added to this 
	 * data set.
	 */
	public synchronized void add(GeneratedFileList gfl) {
		Workspace ws = Workspace.get();
		if (gfl != null) {
			gflList.add(gfl);
			gfl.setDataSet(this);
			WorkspaceEvent wse = new WorkspaceEvent(gfl,
					WorkspaceEvent.Operation.INSERT);
			ws.fireWorkspaceChanged(wse);
		}
	}
	
	/**
	 * Remove an existing GeneratedFileList (GFL) data file from 
	 * this data set.
	 *  
	 * <p><b>Note:</b> This method reports the removed entry to
	 * all workspace listeners by firing a suitable event.</p>
	 * 
	 * @param gfl The MSTData entry to be removed from this data set. 
	 */
	public synchronized void remove(GeneratedFileList gfl) {
		Workspace ws = Workspace.get();
		if (gfl != null) {
			gflList.remove(gfl);
			WorkspaceEvent wse = new WorkspaceEvent(gfl,
					WorkspaceEvent.Operation.DELETE);
			ws.fireWorkspaceChanged(wse);
		}
	}


	/**
	 * Set aggregate statistics information about the cDNA data file.
	 * 
	 * This method can be used to setup the aggregate statistics about
	 * the cDNA fragments associated with this data set. This information
	 * is persisted in the workspace XML providing rapid access to the
	 * aggregate statistics without having to recompute them.
	 * 
	 * @param stats The aggregate statistics about this data set.
	 */
	public void setStats(DataFileStats stats) {
		this.stats = stats;
	}
	
	/**
	 * Set aggregate statistics information about the cDNA data file.
	 * 
	 * This method can be used to obtain the aggregate statistics about
	 * the cDNA fragments associated with this data set. This information
	 * is typically obtained from its persisted copy in the workspace XML 
	 * providing rapid access to the aggregate statistics without having 
	 * to recompute them.
	 * 
	 * @param stats The aggregate statistics about this data set.
	 */
	DataFileStats getStats() {
		return stats;
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
	 * Method to write summary information about the data set.
	 * 
	 * This method is a convenience method that is used by various 
	 * wizards to display summary information about the data set.
	 * The summary information about the data set include information
	 * about the source cDNA data file along with statistics about
	 * the data set.
	 * 
	 * @param sw The summary writer to which the data is to be written.
	 */
	public void summarize(SummaryWriter sw) {
		sw.addSection("Data Set Summary");
		final File tempFile = new File(path);
		sw.addSummary("File Name", tempFile.getName(), description);
		sw.addSummary("Path",      tempFile.getPath(), null);
		sw.addSummary("File type", fileType.toString(), null);
		// Add statistics about the data set if available.
		if (stats != null) {
			stats.summarize(sw);
		}
	}
	
	/**
	 * Helper method to check if the MST file exists and is readable.
	 * 
	 * This is a convenience method that can be used to verify that
	 * the MST file (associated with this entry) exists and is readable.
	 * 
	 * @return This method returns true if the MST file exists and is
	 * readable. Otherwise this method returns false.
	 */
	public boolean isGood() {
		File f = new File(path);
		return (f.exists() && f.canRead());
	}
	
	/**
	 * Helper method to get a tool-tip text for GUI components to use.
	 * 
	 * This is a helper method that provides an HTML formatted tool-tip
	 * text that can be readily displayed by the GUI. The tool-tip
	 * is a multi-line HTML fragment that includes: file path,
	 * file type, description, and statistics (if available).
	 * 
	 * @return A HTML document that contains a HTML-formatted tool-tip
	 * text. This string is never null.
	 */
	public String getToolTipText() {
		final File   tempFile = new File(path);
		final String shortPath= Utilities.trim(tempFile.getPath(), 50);
		final String htmlDesc = Utilities.wrapStringToHTML(description, 45);
		final String formatStr= ((stats != null) ? TOOL_TIP_WITH_STATS_TEMPLATE :
			                                          TOOL_TIP_NO_STATS_TEMPLATE);
		final String toolTip  =
			String.format(formatStr, shortPath, fileType.toString(), 
					htmlDesc, stats.getCount(), stats.getAvgLength(), 
					stats.getLengthSD());
		return toolTip;
	}
	
	/**
	 * The optional list of generated data files that are currently 
	 * associated with this data set. Each entry in this list contains
	 * a list of generated output files/artifacts from jobs run using
	 * this data set (or one of the underlying artifacts).
	 */
	private ArrayList<GeneratedFileList> gflList;
	
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
	
	/**
	 * A statistics object that contains the meta data about the cDNA file
	 * associated with this data set. This element is an optional element
	 * and may not be present in older data set entries created by PEACE.
	 */
	private DataFileStats stats;
	
	/**
	 * A fixed string constant to ease generation of tool tip text
	 * for use/display by GUI components. This text string is
	 * suitably formatted (by the {@link #getToolTipText()} method)
	 * via printf to fill-in values for various parameters.
	 */
	private static final String TOOL_TIP_WITH_STATS_TEMPLATE = "<html>" +
		"<b>Path:</b> %s<br/>" +
		"<b>File Type:</b> %s<br/>" +
		"<b>Description:</b> %s<br/>" + 
		"<b>cDNA/EST Entries:</b> %d<br/>" +
		"<b>Avg. Size per Entry:</b> %.2f (SD: %.2f)" +
		"</html>";
	
	/**
	 * A fixed string constant to ease generation of tool tip text
	 * for use/display by GUI components. This text string is
	 * suitably formatted (by the {@link #getToolTipText()} method)
	 * via printf to fill-in values for various parameters. This
	 * string is used in cases where the data set does not have
	 * pre-computed aggregate statistics associated with it.
	 */
	private static final String TOOL_TIP_NO_STATS_TEMPLATE = "<html>" +
		"<b>Path:</b> %s<br/>" +
		"<b>File Type:</b> %s<br/>" +
		"<b>Description:</b> %s<br/>" + 
		"</html>";
}
