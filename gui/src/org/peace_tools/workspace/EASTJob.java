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

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.workspace.FileEntry.FileEntryType;
import org.w3c.dom.Element;

/**
 * <p>This class corresponds to a "EASTJob" element in a 
 * PEACE work space XML data. This class encapsulates the core 
 * information associated with an assembly job performed via EAST.
 * This class serves merely as a read-only type class that is 
 * created when a new job is run. The data is persisted in the
 * XML work space for future reference so that users can determine
 * jobs that were scheduled and check on status of long-running jobs.
 * A lot of the data associated with this type of job is provided/
 * managed by the base class.</p>
 * 
 * <p>In addition to encapsulating the data this class also provides
 * convenient interfaces for marshaling and un-marshaling XML data compatible
 * with the PEACE GUI configuration XML file format. </p>
 */
public class EASTJob extends Job {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * job entry. This method is typically  used to create a suitable
	 * job entry when loading a Work space into the GUI.
	 * 
	 * @param jobNode The DOM element to be used for creating the job
	 * entry and populating with the needed data. This node must
	 * correspond to "EASTJob" element.
	 * 
	 * @return The newly created job entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static EASTJob create(Element jobNode) throws Exception {
		// First extract the necessary information to create a minimally
		// populated clustering job object from the DOM tree.
		String jobID    = DOMHelper.getStringValue(jobNode, "JobID");
		String srvrID   = DOMHelper.getStringValue(jobNode, "ServerID");
		
		// Now we have the information to create a minimal east
		// job class. So do it and use it for further operations.
		EASTJob job= new EASTJob(jobID, srvrID);
		
		// Now let the job class of job marshal the non-final values.
		job.unmarshal(jobNode);
		// Return the fully populated entry back to the caller
		return job;
	}
		
	@Override
	protected void unmarshal(Element jobNode) throws Exception {
		// Let base class un-marshal common elements
		super.unmarshal(jobNode);
		// Extract the parameters using helper method
		parameters = parseParameters(jobNode);
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
	 */
	private EASTJob(String jobID, String serverID) {
		super(JobType.EAST, jobID, serverID);
		// Save information about heuristics and filters
		this.parameters = null;
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
	 * @param memory The total memory (sum of all memory used by all processes, 
	 * (in MB) that was allocated for this job.
	 * 
	 * @param maxRunTime The maximum runtime (in hours) that was assigned for 
	 * this job.
	 * 
	 * @param parameters The list of parameters that were used to customize
	 * operation of EAST. This list cannot be null (but can be empty). 
	 */
	public EASTJob(String jobID, String description, String serverID,
			String path, int memory, int maxRunTime, ArrayList<Param> parameters) {
		super(JobType.EAST, jobID, serverID, description, path, 
			1, 1, memory, maxRunTime);
		// Save information about parameters
		this.parameters = parameters;
	}
	
	/**
	 * Command line arguments for EAST.
	 * 
	 * This method can be used to obtain the parameter information
	 * in the form of command line parameters that can be readily 
	 * passed to PEACE/EAST.  The command line parameters are used
	 * to configure the filters used by the clustering engine to 
	 * mirror the configuration setup by the user for this job.
	 * 
	 * @return The command line parameters to be passed to the PEACE/EAST. 
	 * If no parameters are defined then this method returns an empty
	 * string 
	 */
	@Override
	public String getParametersCmdLine() {
		return toCmdLine(parameters);
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
		return " " + getParametersCmdLine();
	}
	
	@Override
	public String toString() {
		return jobID;
	}

	/**
	/**
	 * Helper method to add input and output files as parameters to this job.
	 * 
	 * EAST requires the input and output file to be specified as command-line
	 * arguments in a specific order. This method is used by the 
	 * EastJobWizard to add additional input/output files as parameters to
	 * this job.
	 * 
	 * @param ds The data set that contains the source FASTA file to be assembled. 
	 * This parameter cannot be null.
	 * 
	 * @param mstFile The file that contains the MST data to be used for assembly.
	 * This parameter cannot be null.
	 *  
	 * @param qualityFile An optional quality file (in FASTAQ format) to be used
	 * during assembly. This parameter can be null.
	 *  
	 * @param gfl The list of generated output files to be specified as target
	 * output files to EAST.
	 */
	public void addParameters(DataSet ds, FileEntry mstFile, FileEntry qualityFile, 
			GeneratedFileList gfl) {
		parameters.add(new Param(getServerESTFile(ds), null));
		parameters.add(new Param(mstFile.getName(), null));
		
		parameters.add(new Param(gfl.findEntry(FileEntryType.ASM).getName(), null));
		parameters.add(new Param(gfl.findEntry(FileEntryType.SINGLETONS).getName(), null));
		parameters.add(new Param(gfl.findEntry(FileEntryType.STATS).getName(), null));
		
		if (qualityFile != null) {
			parameters.add(new Param(qualityFile.getName(), null));
		}
	}
	
	/**
	 * Obtain the list of parameters that were used to customize this
	 * job. These command line parameters were passed to the final
	 * executable to customize its operations for this job.
	 * 
	 * @return The list of parameters associated with this job. This 
	 * value cannot be null (but can be an empty list)
	 */
	public ArrayList<Param> getParameters() { return parameters; }
	
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
		Element job = DOMHelper.addElement(jobList, "EASTJob", null);
		// Let base class add common sub-elements
		super.marshall(job);
		// Now marshal the information regarding parameters.
		marshallParameters(parameters, job);
	}
	
	/**
	 * Method to marshal the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t";
		// Create a top-level server entry for this server
		out.printf("%s<EASTJob>\n", Indent);
		// Let base class write out common elements first
		super.marshall(out);
		// Now marshal the information regarding parameters.
		marshallParameters(parameters, out);
		// Close the job tag
		out.printf("%s</EASTJob>\n", Indent);
	}

	/**
	 * Method to write summary information about the job.
	 * 
	 * This method is a convenience method that is used by various 
	 * wizards to display summary information about the job.
	 * The summary information about the job include information
	 * about the filter, heuristics, server-configuration etc.
	 * 
	 * @param sw The summary writer to which the data is to be written.
	 * 
	 */
	@Override
	public void summarize(SummaryWriter sw) {
		// Cut summary information about parameters
		if (parameters != null) {
			// Display list of parameters
			sw.addSection("Command line parameter list");
			for(Param p: parameters) {
				sw.addSummary(p.getName(), p.getValue(), null);
			}
		}
		// Let the base class summarize common information
		super.summarize(sw);
	}
		
	/**
	 * The list of parameters configured by the user for this job. These
	 * parameters are passed to corresponding assembly tools.
	 * This array list cannot be null (but can be empty).
	 */
	private ArrayList<Param> parameters;
}
