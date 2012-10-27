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

import java.io.PrintWriter;
import java.util.ArrayList;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.w3c.dom.Element;


/**
 * <p>This class corresponds to a "BatonJob" element in a 
 * PEACE work space XML data. This class encapsulates the core 
 * information associated with a Job. This class serves merely
 * as a read-only type class that is created when a new job is run.
 * The data is persisted in the XML work space for future reference
 * so that users can determine jobs that were scheduled and check on
 * status of long-running jobs.</p>
 * 
 * <p>In addition to encapsulating the data this class also provides
 * convenient interfaces for marshaling and un-marshaling XML data compatible
 * with the PEACE GUI configuration XML file format. </p>
 */
public class BatonJob extends Job {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * Job entry. This method is typically  used to create a suitable
	 * Job entry when loading a Work space into the GUI.
	 * 
	 * @param jobNode The DOM element to be used for creating the job
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created job entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static BatonJob create(Element jobNode) throws Exception {
		// First extract the necessary information to create a minimally
		// populated clustering job object from the DOM tree.
		String jobID    = DOMHelper.getStringValue(jobNode, "JobID");
		String srvrID   = DOMHelper.getStringValue(jobNode, "ServerID");
		// Now we have the information to create a minimal baton
		// job class. So do it and use it for further operations.
		BatonJob job= new BatonJob(jobID, srvrID); 
		// Now let the job class of job marshal the non-final values.
		job.unmarshal(jobNode);
		// Return the fully populated entry back to the caller
		return job;
	}
	
	@Override
	protected void unmarshal(Element jobNode) throws Exception {
		// Let base class un-marshal common elements
		super.unmarshal(jobNode);
		// TODO: implement parameters
		//parameters = parseParameters(jobNode);
	}
	
	/**
	 * Constructor to create a minimally populated clustering job
	 * object with just the fixed value fields initialized to 
	 * specific values.
	 * 
	 * This constructor is typically used to create a temporary job
	 * object into which additional values are going to be loaded
	 * from a XML work space. This constructor is currently used
	 * only by the {@link #create(Element)} method.
	 * 
	 * @param jobID The workspace-wide unique, generated job ID for this job.
	 * The jobID is generated via a call to JobList.reserveJobID() method.
	 * This is a valid string (typically in the form job####)
	 * 
	 * @param serverID The ID of the server on which this job is running.
	 * This is a cross reference ID of a Server configured in this workspace.
	 * 
	 * @param analyzer Information about the distance/similarity analyzer
	 * that was used to create this MST data set.
	 *            
	 * @param threshold The threshold value used to partition a MST.
	 */
	private BatonJob(String jobID, String serverID) {
		super(JobType.CLUSTERING, jobID, serverID);
	}
	
	/**
	 * Set number of nucleotides constituting baton head.
	 * This could be used in the future when filters and parameters are fully
	 * implemented.
	 * @param kmers The number of nucleotides constituting baton head.
	 */
	protected void setKmers(int kmers){
		this.kmers = kmers;
	}
	/**
	 * Constructor to create a Job object with the fixed value fields
	 * initialized to specific values.
	 * 
	 * @param jobID The workspace-wide unique, generated job ID for this job.
	 * The jobID is generated via a call to JobList.reserveJobID() method.
	 * This is a valid string (typically in the form job####)
	 * 
	 * @param description A user-supplied description for this job entry. This
	 * maybe an empty string (but cannot be null).
	 * 
	 * @param serverID The ID of the server on which this job is running.
	 * This is a cross reference ID of a Server configured in this workspace.
	 * 
	 * @param path The directory on the server where the data for this job is stored.
	 * 
	 * @param nodes The number of compute nodes that were requested on a cluster for
	 * running this job. This value must be at least 1.
	 * 
	 * @param cpusPerNode The number of CPUs on each node that were requested for this
	 * job. This value must be at least 1.
	 * 
	 * @param memory The total memory (sum of all memory used by all processes, 
	 * (in MB) that was allocated for this job.
	 * 
	 * @param maxRunTime The maximum runtime (in hours) that was assigned for 
	 * this job.
	 * 
	 * @param analyzer Information about the distance/similarity analyzer
	 * that was used to create this MST data set.
	 *            
	 * @param threshold The threshold value used to partition a MST.
	 * 
	 * @param heuristics The list of heuristics that were used to accelerate 
	 * the MST generation algorithm. This list may be null for no heuristics at all
	 * and the defaults will kick in. If this list not-null but empty then 
	 * all heuristics are disabled.
	 * 
	 * @param filters The list of filters that were used to filter out short ESTs or
	 * ESTs with low complexity sections to ensure that the resultant clustering is
	 * high quality. This list may be null.
	 * 
	 * @param parameters The list of general purpose parameters to be set
	 * for this job. These parameters are used to further customize the
	 * operations of the clustering engine.
	 */
	public BatonJob(String jobID, String description, String serverID,
			String path, int nodes, int cpusPerNode,
			int memory, int maxRunTime, ArrayList<Param> parameters) {
		super(JobType.CLUSTERING, jobID, serverID, description, path,
			nodes, cpusPerNode, memory, maxRunTime);
		this.parameters = parameters;
		this.dataFileName = path;
	}

	
	/**
	 * Helper method to build the set of general purpose command-line
	 * parameters associated with this job.
	 * 
	 * This is a helper method creates the list of general purpose
	 * command-line parameters associated with this job. These
	 * parameters are not associated with filters or heuristics. However,
	 * they are used to further customize the operations of the 
	 * clustering engine.
	 */
	public void setupParameters(DataSet ds, GeneratedFileList gfl,
			boolean maskBases) {
		// Add command-line for input data set file type. Note that
		// here we map the path from local machine to remote machine
		String fileType = (ds.getFileType().equals(DataFileType.SFF) ?
				"--sffFile" : "--fastaFile");
		parameters.add(new Param(fileType, getServerESTFile(ds)));
		// Check and suppress base masking.
		if (!maskBases) {
			// Disable base masking. That is atcg is
			// equal to ATCG.
			parameters.add(new Param("--no-mask-bases", null));
		}
		// Cut the file names with appropriate parameter entries.
		FileEntry mstFE = gfl.findEntry(FileEntry.FileEntryType.MST);
		parameters.add(new Param("--output-mst-file", mstFE.getName()));
		// Add entry for clustering output in GUI compatible format
		parameters.add(new Param("--gui-print", null));
		FileEntry clsFE = gfl.findEntry(FileEntry.FileEntryType.CLS);
		parameters.add(new Param("--output-cls-file", clsFE.getName()));
		// Add command line option to generate progress information
		parameters.add(new Param("--progress", "progress.dat"));
	}
	
	
	/**
	 * Return the information in the form of a partial PEACE command line 
	 * parameter.
	 * 
	 * This method can be used to obtain the information needed to
	 * generate the MST and clusters based on the supplied information in the form
	 * of a command line parameter. 
	 * 
	 * @return Return the information as a command line to be passed to the
	 * PEACE clustering engine.
	 */
	public String toCmdLine() {
	// Use helper method to build the command line.
	//TODO: Fix command line arguments: 
	// When filters are implemented, they may be added to the cmd string	
		return " peace --assembler baton --analyzer null --estFiles " + dataFileName;
	}
	
	@Override
	public String toString() {
		return jobID;
	}
	
	/**
	 * Method to marshal the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent JobList node in the DOM tree.
	 * 
	 * @param jobList The DOM element corresponding to the "JobList"
	 * node that contains this entry.
	 */

	public final void marshall(Element jobList) {
		// Create a top-level entry for this "Job"
		Element job = DOMHelper.addElement(jobList, "BatonJob", null);
		// Let base class add common sub-elements
		super.marshall(job);
		// Now marshal the information regarding k-mers.
		if (kmers == 2 || kmers == 3) {
		  DOMHelper.addElement(job, "K-mers", "" + kmers);
		}
		// Now marshal the information regarding general parameters
		//marshallParameters(parameters, job);
	}
	
	/**
	 * Method to marshal the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 * 
	 * @param indentPrefix The extra indentation to be done to make the
	 * output look nice. If no additional indentation is needed then
	 * an empty string ("") must be passed in.
	 */
	public final void marshall(PrintWriter out, final String indentPrefix) {
		final String Indent = indentPrefix + "\t\t";
		final String NUM_ELEMENT = Indent + "\t" + "<%1$s>%2$d</%1$s>\n";
		// Create a top-level server entry for this server
		out.printf("%s<BatonJob>\n", Indent);
		// Let base class write out common elements first
		super.marshall(out, indentPrefix);
		// Marshal out the k-mer information
		out.printf(NUM_ELEMENT, "K-mers", kmers );
		// Close the job tag
		out.printf("%s</BatonJob>\n", Indent);
	}
	
	/**
	 * The list of filters configured by the user for this job. These
	 * filters are used to weed out fragments that may interfere with
	 * clustering and reduce overall quality of clustering.
	 * This array list can be null (depending on the job type).
	 */
	private int kmers;
		
	/**
	 * The list of parameters configured by the user for this job. These
	 * parameters are passed to corresponding assembly tools.
	 * This array list cannot be null (but can be empty).
	 */
	private ArrayList<Param> parameters;
	
	/**
	 * This string is used to pass the filename of the data set for a 
	 * baton job to the getCmd() function. This is a temporary fix,
	 * and will disappear when parameters have been set up correctly.  
	 */
	private String dataFileName;
}
