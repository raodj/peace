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

import javax.xml.datatype.DatatypeFactory;

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.generic.Utilities;
import org.w3c.dom.Element;


/**
 * Class to encapsulate some basic information about a Job. This class is used
 * within the MSTData and MSTClusterData classes. This class serves merely as a
 * read-only type class that is created when a new job is run. The data is
 * persisted in the XML work space for future reference so that users can
 * determine some core characteristics of the Job that was run to generate a
 * given data set.
 */
public class JobSummary extends JobBase {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * JobSummary entry. This method is typically  used to create a suitable
	 * JobSummary entry when loading a Work space into the GUI.
	 * 
	 * @param jobNode The DOM element to be used for creating the summary
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created job summary entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static JobSummary create(Element jobNode) throws Exception {
		// First extract the necessary information from the DOM tree.
		String jobID      = jobNode.getAttribute("jobID");
		String jTypeStr   = DOMHelper.getStringValue(jobNode, "Type");
		String statusStr  = DOMHelper.getStringValue(jobNode, "Status");
		int    cpus       = DOMHelper.getIntValue(jobNode, "CPUs");
		// Obtain the server name and serverID attribute
		Element server    = DOMHelper.getElement(jobNode, "ServerName");
		String srvrName   = DOMHelper.getStringValue(jobNode, "ServerName");
		srvrName          = srvrName.trim();
		String srvrID     = server.getAttribute("serverID");
		// Extract the heuristics and filter summary strings.
		String heuristics = DOMHelper.getStringValue(jobNode, "HeuristicsSummary", "");
		String filters    = DOMHelper.getStringValue(jobNode, "FiltersSummary", "");
		String parameters = DOMHelper.getStringValue(jobNode, "ParametersSummary", "");
		String prevJobID  = null;
		if (DOMHelper.hasElement(jobNode, "PrevJobID")) {
			prevJobID  = DOMHelper.getStringValue(jobNode, "PrevJobID");
		}
		// Now that we have sufficient information create the job summary
		final JobBase.JobType type = JobBase.JobType.valueOf(jTypeStr.toUpperCase());
		final JobBase.JobStatusType status = JobStatusType.valueOf(statusStr.toUpperCase());
		JobSummary job = new JobSummary(type, jobID, srvrID, cpus, 
				srvrName, heuristics, filters, parameters, prevJobID);
		job.setStatus(status);
		// Now update its various other if available.
		if (DOMHelper.hasElement(jobNode, "RunTime")) { 
			String runTime  = DOMHelper.getStringValue(jobNode, "RunTime");
			DatatypeFactory codec = DatatypeFactory.newInstance();
			job.runtime = codec.newDuration(runTime);
		}

		// Return the newly created object.
		return job;
	}
	
	/**
	 * Constructor to create a common job object with the fixed value fields
	 * initialized to specific values.
	 * 
	 * @param type The type of job entry for which this summary is being
	 * created.
	 * 
	 * @param jobID The work space-wide unique, generated job ID for this job.
	 * The jobID is generated via a call to JobList.reserveJobID() method.
	 * This is a valid string (typically in the form job####)
	 * 
	 * @param serverID The ID of the server on which this job is running.
	 * This is a cross reference ID of a Server configured in this work space.
	 * 
	 * @param cpus The total number of CPUs that were used to run this job.
	 * 
	 * @param serverName The name (or IP address) of the server on which this job
	 * was run. 
	 * 
	 * @param heuristicsSummary The command line parameter(s) passed to the PEACE
	 * clustering engine (C++ side) to configure heuristics to accelerate clustering.
	 * 
	 * @param filtersSummary The command line parameter(s) passed to the PEACE
	 * clustering engine (C++ side) to configure filters to improve clustering quality.
	 * 
	 * @param parameters The command line parameter(s) passed to PEACE/EAST. This
	 * value may be an empty string but not null.
	 * 
	 * @param prevJobID The unique generated jobID value for another/previous
	 * job that this job is dependent on. This value is typically set when a 
	 * clustering + assembly type jobs are scheduled. The assembly job
	 * is flagged as being dependent on the clustering job. This value can be
	 * null.
	 */
	public JobSummary(JobBase.JobType type, String jobID, String serverID, 
			int cpus, String serverName, String heuristicsSummary, 
			String filtersSummary, String parameters, String prevJobID) {
		super(type, jobID, serverID);
		this.cpus              = cpus;
		this.serverName        = serverName;
		this.heuristicsSummary = heuristicsSummary;
		this.filtersSummary    = filtersSummary;
		this.parameters        = parameters;
		this.prevJobID         = prevJobID;
	}
	
	/** A convenience constructor.
	 * 
	 * The constructor uses the complete Job object to populate the various 
	 * summary information that is contained in this class.
	 * 
	 * @param job The job object to be used for summarizing the information 
	 * associated with this job.
	 */
	public JobSummary(Job job) {
		super(job);
		// Set up the non-common data members next.
		this.cpus           = job.getNodes() * job.getCPUsPerNode();
		// Get the server name indirectly using the serverID
		String srvrID       = job.getServerID();
		Workspace workspace = Workspace.get();
		Server server       = workspace.getServerList().getServer(srvrID);
		this.serverName     = server.getName();
		// Save heuristics and filter summary information separately.
		heuristicsSummary   = job.getHeuristicsCmdLine();
		filtersSummary      = job.getFiltersCmdLine();
		parameters          = job.getParametersCmdLine();
	}
	
	/** A convenience constructor.
	 * 
	 * The constructor uses the complete Job object to populate the various 
	 * summary information that is contained in this class.
	 * 
	 * @param job The job object to be used for summarizing the information 
	 * associated with this job.
	 * 
	 * @param analyzerType The type of analyzer used for the job. If the analyzer's type
	 * is TwoPassD2 that uses automatic heuristics, then the heuristic summary
	 * is set to an empty string.
	 */
	public JobSummary(Job job, FWAnalyzer.FWAnalyzerType analyzerType) {
		super(job);
		// Set up the non-common data members next.
		this.cpus           = job.getNodes() * job.getCPUsPerNode();
		// Get the server name indirectly using the serverID
		String srvrID       = job.getServerID();
		Workspace workspace = Workspace.get();
		Server server       = workspace.getServerList().getServer(srvrID);
		this.serverName     = server.getName();
		// Save heuristics and filter summary information separately.
		heuristicsSummary   = job.getHeuristicsCmdLine();
		filtersSummary      = job.getFiltersCmdLine();
		parameters          = job.getParametersCmdLine();
	}
	
	/**
	 * The server name that was obtained from the Server entry and stored
	 * in this summary when the JobSummary object was created.
	 * 
	 * @return The domain name (or IP address) of the server on which the
	 * Job was run.
	 */
	public String getServerName() { return serverName; }
	
	/**
	 * The total number of CPUS that was computed and stored
	 * in this summary when the JobSummary object was created.
	 * 
	 * @return The total number of CPUs that were used for this job.
	 */
	public int getCPUs() { return cpus; }
	
	/**
	 * Obtain a summary of the heuristics run as a part of the data.
	 * 
	 * This method can be used to obtain a summary representation of the
	 * heuristics used in the job to generate the associated data.
	 * 
	 * @return A summary information about the heuristics used. This 
	 * information is stored as the command line parameters passed to the
	 * PEACE clustering engine. This value can be an empty string.
	 */
	public String getHeuristicsSummary() { return heuristicsSummary; }

	/**
	 * Obtain a summary of the filters run as a part of the data.
	 * 
	 * This method can be used to obtain a summary representation of the
	 * filters used in the job to generate the associated data.
	 * 
	 * @return A summary information about the filters used. This 
	 * information is stored as the command line parameters passed to the
	 * PEACE clustering engine. This value can be an empty string.
	 */
	public String getFiltersSummary() { return filtersSummary; }

	/**
	 * Obtain a summary of the parameters run as a part of the job.
	 * 
	 * This method can be used to obtain a summary representation of the
	 * parameters used in the job to generate the associated data.
	 * 
	 * @return A summary information about the parameters used. This 
	 * information is stored as the command line parameters passed to the
	 * PEACE clustering engine. This value can be just an empty
	 * string.
	 */
	public String getParameterSummary() { return parameters; }
	
	/**
	 * Method to marshal the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent JobList node in the DOM tree.
	 * 
	 * @param jobList The DOM element corresponding to the "JobList"
	 * node that contains this entry.
	 */
	public final void marshall(Element jobList) {
		// Create a top-level entry for this "JobSummary"
		Element job = DOMHelper.addElement(jobList, "JobSummary", null);
		// Add the type attributes for this server.
		job.setAttribute("jobID", jobID);
		// Add the type of job to the node
		DOMHelper.addElement(job, "Type", type.toString().toLowerCase());
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(job, "Status", status.toString().toLowerCase());
		// Setup the serverName element with serverID as attribute
		Element server = DOMHelper.addElement(job, "ServerName", serverName);
		server.setAttribute("serverID", serverID);
		DOMHelper.addElement(job, "CPUs", "" + cpus);
		if (runtime != null) {
			DOMHelper.addElement(job, "RunTime", runtime.toString());
		}
		// Marshal the heuristics and filter summary information as well.
		DOMHelper.addElement(job, "HeuristicsSummary", heuristicsSummary);
		DOMHelper.addElement(job, "FiltersSummary",    filtersSummary);
		DOMHelper.addElement(job, "ParametersSummary", parameters);
		if (prevJobID != null) {
			DOMHelper.addElement(job, "PrevJobID", prevJobID);
		}
	}
	
	/**
	 * Method to marshal the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t\t";
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		final String NUM_ELEMENT = Indent + "\t" + "<%1$s>%2$d</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<JobSummary jobID=\"%s\">\n", Indent, jobID); 
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "Type", type.toString().toLowerCase());
		out.printf(STR_ELEMENT, "Status", status.toString().toLowerCase());
		out.printf("%s\t<%2$s serverID=\"%3$s\">%4$s</%2$s>\n", Indent, 
				"ServerName", serverID, serverName);
		out.printf(NUM_ELEMENT, "CPUs", cpus);
		if (runtime != null) {
			out.printf(STR_ELEMENT, "RunTime", runtime.toString());
		}
		// Display heuristic summary and filter summary strings.
		out.printf(STR_ELEMENT, "HeuristicsSummary", heuristicsSummary);
		out.printf(STR_ELEMENT, "FiltersSummary",    filtersSummary);
		out.printf(STR_ELEMENT, "ParametersSummary", parameters);
		if (prevJobID != null) {
			out.printf(STR_ELEMENT, "PrevJobID", prevJobID);
		}
		// Close the job tag
		out.printf("%s</JobSummary>\n", Indent);
	}
	
	/**
	 * Method to write summary information about the Job Summary data.
	 * 
	 * This method is a convenience method that is used by various 
	 * wizards to display summary information about the MST data.
	 * The summary information about the data set include information
	 * about the MST data file location and the analyzer used to 
	 * generate the MST.
	 * 
	 * @param sw The summary writer to which the data is to be written.
	 */
	public void summarize(SummaryWriter sw) {
		// See if we can get description from the full job entry.
		Job job = Workspace.get().getJobList().getjob(jobID);
		String description = (job != null ? job.getDescription() : "");
		// Summarize information about this job
		sw.addSubSection("Job Summary", null, description);
		sw.addSubSummary("Server name", serverName, null);
		sw.addSubSummary("#Processes", "" + cpus, 
				"This value is CPUs * Nodes/CPU");
	}
	
	/**
	 * Helper method to get a tool-tip text for GUI components to use.
	 * 
	 * This is a helper method that provides an HTML formatted tool-tip
	 * text that can be readily displayed by the GUI. The tool-tip
	 * is a multi-line HTML fragment that includes some summary 
	 * information about the job that generated the files.
	 * 
	 * @return A HTML document that contains a HTML-formatted tool-tip
	 * text. This string is never null.
	 */
	public String getToolTipText() {
		final String NA = "<i><font size=\"-2\">UNAVAILABLE</font></i>";
		// Get main job entry and job description if available.
		Job job = Workspace.get().getJobList().getjob(jobID);
		String description = (job != null ? job.getDescription() : NA);
		final String htmlDesc = Utilities.wrapStringToHTML(description, 45);
		String lastUpdate  = (job != null ? job.getLastUpdateTimestamp().toString() : NA);
		
		// Create the formatted tool-tip
		final String toolTip  =
			String.format(TOOL_TIP_TEMPLATE, jobID, type.toString(), 
					htmlDesc, status.toString(), serverID, serverName, lastUpdate);
		return toolTip;
	}
	
	/**
	 * The domain name (or IP address) of the server on which the job was actually
	 * run.
	 */
	private final String serverName;

	/**
	 * The total number of CPUs (Nodes * CPUsPerNode) that were used for
	 * running the job.
	 */
	private final int cpus;
	
	/**
	 * A string that contains summary of heuristics run as a part of the job
	 * used to generate the associated data. This string is simply stored as
	 * the command line parameter passed to the PEACE clustering tool.
	 */
	private final String heuristicsSummary;
	
	/**
	 * A string that contains summary of filters run as a part of the job
	 * used to generate the associated data. This string is simply stored as
	 * the command line parameter passed to the PEACE clustering tool.
	 */
	private final String filtersSummary;
	
	/**
	 * A string that contains summary of parameters set up as a part of the job
	 * used to generate the associated data. This string is simply stored as
	 * the command line parameter passed to the PEACE/EAST tools.
	 */
	private final String parameters;
	
	/**
	 * A fixed string constant to ease generation of tool tip text
	 * for use/display by GUI components. This text string is
	 * suitably formatted (by the {@link #getToolTipText()} method)
	 * via printf to fill-in values for various parameters.
	 */
	private static final String TOOL_TIP_TEMPLATE = "<html>" +
		"<b>Job ID:</b> %s<br/>" +
		"<b>Job Type:</b> %s<br/>" +
		"<b>Description:</b> %s<br/>" +
		"<b>Status:</b> %s<br/>" +
		"<b>Run on server [ID: %s]:</b> %s<br/>" +
		"<b>Status updated on:</b> %s" +
		"</html>";
}
