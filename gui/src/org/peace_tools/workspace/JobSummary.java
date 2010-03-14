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
		int    cpus       = DOMHelper.getIntValue(jobNode, "CPUs");
		// Obtain the server name and serverID attribute
		Element server    = DOMHelper.getElement(jobNode, "ServerName");
		String srvrName   = DOMHelper.getStringValue(jobNode, "ServerName");
		srvrName          = srvrName.trim();
		String srvrID     = server.getAttribute("serverID");
		// Extract the heuristics and filter summary strings.
		String heuristics = DOMHelper.getStringValue(jobNode, "HeuristicsSummary");
		String filters    = DOMHelper.getStringValue(jobNode, "FiltersSummary");
		// Now that we have sufficient information create the job summary
		JobSummary job = new JobSummary(jobID, srvrID, cpus, srvrName, heuristics, filters);
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
	 */
	public JobSummary(String jobID, String serverID, int cpus, String serverName,
					  String heuristicsSummary, String filtersSummary) {
		super(jobID, serverID);
		this.cpus              = cpus;
		this.serverName        = serverName;
		this.heuristicsSummary = heuristicsSummary;
		this.filtersSummary    = filtersSummary;
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
		this.cpus       = job.getNodes() * job.getCPUsPerNode();
		// Get the server name indirectly using the serverID
		String srvrID   = job.getServerID();
		Workspace workspace = Workspace.get();
		Server server       = workspace.getServerList().getServer(srvrID);
		this.serverName     = server.getName();
		// Save heuristics and filter summary information separately.
		heuristicsSummary   = job.getHeuristicsCmdLine();
		filtersSummary      = job.getFiltersCmdLine();
	}
	
	/** A convenience constructor.
	 * 
	 * The constructor uses the complete Job object to populate the various 
	 * summary information that is contained in this class.
	 * 
	 * @param job The job object to be used for summarizing the information 
	 * associated with this job.
	 * 
	 * @param analyzer The type of analyzer used for the job. If the analyzer's type
	 * is TwoPassD2 that uses automatic heuristics, then the heuristic summary
	 * is set to an empty string.
	 */
	public JobSummary(Job job, FWAnalyzer.FWAnalyzerType analyzerType) {
		super(job);
		// Set up the non-common data members next.
		this.cpus       = job.getNodes() * job.getCPUsPerNode();
		// Get the server name indirectly using the serverID
		String srvrID   = job.getServerID();
		Workspace workspace = Workspace.get();
		Server server       = workspace.getServerList().getServer(srvrID);
		this.serverName     = server.getName();
		// Save heuristics and filter summary information separately.
		boolean isTwoPassD2 = analyzerType.equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2); 
		heuristicsSummary   = isTwoPassD2 ? "" : job.getHeuristicsCmdLine();  
		filtersSummary      = job.getFiltersCmdLine();
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
	 * PEACE clustering engine.
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
	 * PEACE clustering engine.
	 */
	public String getFiltersSummary() { return filtersSummary; }
	
	/**
	 * Method to marshall the data stored in this object to become part of
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
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
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
		// Close the job tag
		out.printf("%s</JobSummary>\n", Indent);
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
}
