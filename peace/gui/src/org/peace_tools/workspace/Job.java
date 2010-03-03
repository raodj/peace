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
import java.util.GregorianCalendar;

import javax.xml.datatype.DatatypeConfigurationException;
import javax.xml.datatype.DatatypeFactory;
import javax.xml.datatype.XMLGregorianCalendar;

import org.peace_tools.generic.ProgrammerLog;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * This class corresponds to a "Job" element in a PEACE work space
 * XML data. This class encapsulates the core information associated
 * with a Job. This class serves merely as a read-only type class 
 * that is created when a new job is run. The data is persisted in
 * the XML work space for future reference so that users can determine
 * jobs that were scheduled and check on status of long-running jobs.
 */
public class Job extends JobBase {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * Job entry. This method is typically  used to create a suitable
	 * Job entry when loading a Work space into the GUI.
	 * 
	 * @param jobNode The DOM element to be used for creating the Job
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created job entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static Job create(Element jobNode) throws Exception {
		// First extract the necessary information from the DOM tree.
		String jobID    = DOMHelper.getStringValue(jobNode, "JobID");
		String desc     = DOMHelper.getStringValue(jobNode, "Description", true);
		desc            = (desc != null) ? desc : "";
		String serverID = DOMHelper.getStringValue(jobNode, "ServerID");
		String path     = DOMHelper.getStringValue(jobNode, "Path");
		int    nodes    = DOMHelper.getIntValue(jobNode, "Nodes");
		int    cpus     = DOMHelper.getIntValue(jobNode, "CPUsPerNode");
		int    memory   = DOMHelper.getIntValue(jobNode, "Memory");
		int    maxTime  = DOMHelper.getIntValue(jobNode, "MaxRunTime");
		String lastTime = DOMHelper.getStringValue(jobNode, "LastUpdateTimestamp");
		String statStr  = DOMHelper.getStringValue(jobNode, "Status");
		// Extract the heuristic and filter lists using helper methods.
		ArrayList<Heuristic> heuristics = parseHeuristicChain(jobNode);
		ArrayList<Filter>    filters    = parseFilterChain(jobNode);
		// Now that we have sufficient information create the core
		// job object.
		Job job = new Job(jobID, desc, serverID, path, nodes, cpus, memory, maxTime, 
				heuristics, filters);
		// Now update its various properties
		statStr = statStr.toUpperCase();
		job.status = JobBase.JobStatusType.valueOf(JobBase.JobStatusType.class, statStr);
		// Convert the various strings to suitable data types
		DatatypeFactory codec   = DatatypeFactory.newInstance();
		job.lastUpdateTimestamp = codec.newXMLGregorianCalendar(lastTime);
		if (DOMHelper.hasElement(jobNode, "StartTimestamp")) {
			String startTime= DOMHelper.getStringValue(jobNode, "StartTimestamp");
		    job.startTimestamp = codec.newXMLGregorianCalendar(startTime);
		}
		if (DOMHelper.hasElement(jobNode, "Runtime")) {
			String runTime  = DOMHelper.getStringValue(jobNode, "RunTime");
			job.runtime             = codec.newDuration(runTime);
		}
		// Return the newly created object.
		return job;
	}
	
	/**
	 * Helper method to unmarshall the heuristic list into an in-memory array list.
	 * 
	 * This is a helper method that is invoked from the the create() method to 
	 * unmarshall the XML DOM tree corresponding to the heuristic list node into
	 * suitable in-memory classes for further processing.
	 * 
	 * @param jobNode The top-level job node from where the heuristic chain element
	 * is to be extracted and processed. 
	 * 
	 * @return An array list containing a list of Heuristic objects. 
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	private static ArrayList<Heuristic> parseHeuristicChain(Element jobNode) throws Exception {
		// Parse out the heuristic chain using a helper method.
		Element chainNode= DOMHelper.getElement(jobNode, "HeuristicChain");
		NodeList chain   = chainNode.getElementsByTagName("Heuristic");
		ArrayList<Heuristic> heuristics = new ArrayList<Heuristic>(); 
		for(int idx = 0; (idx < chain.getLength()); idx++) {
			Element node = (Element) chain.item(idx);
			// Create a heuristic entry using the DOM data
			Heuristic heuristic = Heuristic.create(node);
			heuristics.add(heuristic);
		}
		// Return the parsed-in heuristic chain for further use.
		return heuristics;
	}
	
	/**
	 * Helper method to unmarshall the filter list into an in-memory array list.
	 * 
	 * This is a helper method that is invoked from the the create() method to 
	 * unmarshall the XML DOM tree corresponding to the filter list node into
	 * suitable in-memory classes for further processing.
	 * 
	 * @param jobNode The top-level job node from where the filter chain element
	 * is to be extracted and processed. 
	 * 
	 * @return An array list containing the list of filter objects in the 
	 * filter chain. 
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	private static ArrayList<Filter> parseFilterChain(Element jobNode) throws Exception {
		// Parse out the heuristic chain using a helper method.
		Element chainNode= DOMHelper.getElement(jobNode, "FilterChain");
		NodeList chain   = chainNode.getElementsByTagName("Filter");
		ArrayList<Filter> filters = new ArrayList<Filter>(); 
		for(int idx = 0; (idx < chain.getLength()); idx++) {
			Element node = (Element) chain.item(idx);
			// Create a heuristic entry using the DOM data
			Filter heuristic = Filter.create(node);
			filters.add(heuristic);
		}
		// Return the parsed-in heuristic chain for further use.
		return filters;
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
	 * @param memory The total memory (sum of all memory used by all processes, in MB) 
	 * that was allocated for this job.
	 * 
	 * @param maxRunTime The maximum runtime (in hours) that was assigned for this job.
	 * 
	 * @param heuristics The list of heuristics that were used to accelerate 
	 * the MST generation algorithm.
	 * 
	 * @param filters The list of filters that were used to filter out short ESTs or
	 * ESTs with low complexity sections to ensure that the resultant clustering is
	 * high quality.
	 */
	public Job(String jobID, String description, String serverID,
			String path, int nodes, int cpusPerNode,
			int memory, int maxRunTime, ArrayList<Heuristic> heuristics, 
			ArrayList<Filter> filters) {
		super(jobID, serverID);
		this.description         = description;
		this.path                = path;
		this.nodes               = nodes;
		this.cpusPerNode         = cpusPerNode;
		// Set non-final fields to null for now.
		this.startTimestamp      = null;
		this.lastUpdateTimestamp = null;
		this.runtime             = null;
		this.status              = JobStatusType.STARTING;
		this.memory              = memory;
		this.maxRunTime          = maxRunTime;
		// Set default last update timestamp
		setLastUpdateTime();
		// Save information about heuristics and filters
		this.heuristics          = heuristics;
		this.filters             = filters;
	}

	/**
	 * Obtain the user supplied description for this this job. The
	 * value is set when new jobs are scheduled. The description is
	 * persisted in the work space configuration file and loaded when a 
	 * work space is opened in the GUI.
	 * 
	 * @return This method returns the user supplied description
	 * set for this job.
	 */
	public String getDescription() { return description; }

	/**
	 * Set a description for this job. 
	 * 
	 * The description is persisted in the work space configuration
	 * file and loaded when a work space is opened in the GUI.
	 * @param desc
	 */
	public void setDescription(String desc) {
		description = desc;
	}
	
	/**
	 * Obtain the system generated directory where the output files for
	 * this job are stored on the server. This value is set when new 
	 * jobs are scheduled. The path for the job is persisted in the 
	 * work space configuration file and loaded when a 
	 * work space is opened in the GUI.
	 * 
	 * @return This method returns the directory where the files for this
	 * job are located.
	 */
	public String getPath() { return path; }

	/**
	 * Set the path where the files for this job are located.
	 * 
	 * @param path The path where the files for this job are 
	 * located.
	 */
	public void setPath(String path) {
		this.path = path;
	}
	
	/**
	 * Obtain the number of nodes that were requested for running this job.
	 * 
	 * @return The number of nodes that were requested for running this job.
	 */
	public int getNodes() { return nodes; }

	/**
	 * Obtain the number of CPUs on each node that were requested for running
	 * this job.
	 * 
	 * @return The number of CPUs per node that were requested for running this
	 *         job.
	 */
	public int getCPUsPerNode() { return cpusPerNode; }

	/**
	 * Obtain the memory to be requested for this job.
	 * 
	 * @return The memory to be requested for this job. This value is
	 * in Gigabytes.
	 */
	public int getMemory() { return memory; }
	
	/**
	 * Obtain the total time to be requested for this job.
	 * 
	 * This method must be used to determine the maximum run time
	 * to be requested when submitting this job via a job scheduling
	 * system.
	 * 
	 * @return The maximum wall clock time to be requested when 
	 * submitting this job via PBS. This value is in hours. 
	 */
	public int getMaxRunTime() { return maxRunTime; }
	
	/**
	 * This is a helper method that is used to set the lastUpdateTimestamp
	 * to the current system time.
	 */
	private void setLastUpdateTime() {
		try {
			DatatypeFactory codec = DatatypeFactory.newInstance();
			GregorianCalendar now = new GregorianCalendar();
			lastUpdateTimestamp = codec.newXMLGregorianCalendar(now);
		} catch (DatatypeConfigurationException e) {
			// Cut log entry in programmer log.
			ProgrammerLog.log(e);
		}
	}
	
	/**
	 * Change the status for this job. This method also sets the last
	 * update time stamp.
	 * 
	 * @param status The new status value to be set for this job.
	 */
	@Override
	public void setStatus(JobStatusType status) {
		super.setStatus(status);
		// Record the last time the status was updated
		setLastUpdateTime();
		// Notify all the listeners about the status change.
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(this, WorkspaceEvent.Operation.UPDATE);
		Workspace.get().fireWorkspaceChanged(we);
		// Update status of any job summaries in the data set
		MSTData mst = Workspace.get().getMSTData(jobID);
		if (mst != null) {
			mst.updateJobSummary(this);
			// Check and update the cluster file as needed
			MSTClusterData cluster = mst.getDataSet().getClusterData(jobID);
			if (cluster != null) {
				cluster.updateJobSummary(this);
			}
		}
	}
	
	/**
	 * Returns the last time the status of this job was set/changed.
	 * 
	 * @return The time stamp when the status of this job was last set.
	 */
	public XMLGregorianCalendar getLastUpdateTimestamp() { 
		return lastUpdateTimestamp; 
	}

	/**
	 * This method can be used to explicitly set the time stamp when this
	 * job actually started running on the server.
	 * 
	 * @param timestamp The timestamp when the job actually started running
	 * on the server. This information is typically determined from the output
	 * logs generated by the job.
	 */
	public void setStartTimestamp(GregorianCalendar timestamp) {
		try {
			DatatypeFactory codec = DatatypeFactory.newInstance();
			startTimestamp = codec.newXMLGregorianCalendar(timestamp);
		} catch (DatatypeConfigurationException e) {
			// Cut log entry in programmer log.
			e.printStackTrace();
		}
	}
	
	/**
	 * Obtain the timestamp indicating when the job actually started running on
	 * the server.
	 * 
	 * @return The timestamp when the job started running. If the job has not
	 * yet started running, then this method returns null.
	 */
	public XMLGregorianCalendar getStartTimestamp() {
		return startTimestamp;
	}

	/**
	 * Command line for PEACE clustering engine to configure heuristics.
	 * 
	 * <p>This method can be used to obtain the heuristic information in the form
	 * of command line parameters that can be readily passed to the PEACE
	 * clustering engine. The command line parameters are used to configure the
	 * clustering engine to suit the configuration setup by the user for this
	 * job.</p>
	 * 
	 * <p><b>Note:<b> The return value of this method must be ignored if this
	 * call is being made in conjunction with the two pass d2 analyzer.</p> 
	 * 
	 * @return The command line (to correspondingly setup the heuristics) to
	 * be passed to the PEACE clustering engine.
	 */
	public String getHeuristicsCmdLine() {
		String cmdLine = "";				
		// Convert heuristic information to command line parameters.
		String heuristicParams = "";
		cmdLine += "--heuristics ";
		if (heuristics.size() > 0) {
			for(int i = 0; (i < heuristics.size()); i++) {
				Heuristic heur = heuristics.get(i);
				cmdLine += (i > 0 ? "-" : "") + heur.getName();
				heuristicParams += heur.toCmdLine();
			}
		} else {
			// The command line parameter for no heuristics case is null
			cmdLine += "null";
		}
		cmdLine += heuristicParams;
		// Return the cmd line for the heuristics
		return cmdLine;
	}

	/**
	 * Command line for PEACE clustering engine to configure filters.
	 * 
	 * This method can be used to obtain the filters information in the form
	 * of command line parameters that can be readily passed to the PEACE
	 * clustering engine. The command line parameters are used to configure the
	 * filters used by the clustering engine to mirror the configuration setup 
	 * by the user for this job.
	 * 
	 * @return The command line (to correspondingly setup the filters) to
	 * be passed to the PEACE clustering engine.
	 */
	public String getFiltersCmdLine() {
		// Next covert filter information to command line parameters.
		String filterParams = "";
		String cmdLine = "--filters ";
		if (filters.size() > 0) {
			for(int i = 0; (i < filters.size()); i++) {
				Filter filter = filters.get(i);
				cmdLine += (i > 0 ? "-" : "") + filter.getName();
				filterParams += filter.toCmdLine();
			}
		} else {
			// The command line parameter for no heuristics case is null
			cmdLine += "null";
		}
		cmdLine += filterParams;
		// Return the cmd line.
		return cmdLine;
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
		return " " + getHeuristicsCmdLine() + " " + getFiltersCmdLine();
	}

	/**
	 * Set the progress status for this job. This method generates 
	 * notifications to all registered workspace listeners.
	 * 
	 * @param estsAnalyzed The number of ESTs that have been analyzed.
	 * @param totalEstCount The total number of ESTs to be analyzed.
	 */
	public void setProgress(int estsAnalyzed, int totalEstCount) {
		progressInfo[0] = estsAnalyzed;
		progressInfo[1] = totalEstCount;
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(this, WorkspaceEvent.Operation.UPDATE);
		Workspace.get().fireWorkspaceChanged(we);
	}
	
	/**
	 * Obtain the current progress information status for this job.
	 * 
	 * <p><b>Note:</b>  The progress information returned may be -1, -1 
	 * if the progress information is not known.</p>
	 * 
	 * @return The progress information for the job. The first 
	 * entry is the number of ESTs analyzed thus far and the second
	 * entry is the total number of ESTs to be analyzed.
	 */
	public int[] getProgressInfo() { return progressInfo; }
	
	@Override
	public String toString() {
		return jobID;
	}
	
	/**
	 * Obtain the list of heuristics that were used to accelerate the
	 * process of constructing the MST. Specifically many of these
	 * heuristics accelerate the frame-word analyzer that was used to
	 * build the MST.
	 * 
	 * @return The list of heuristics that were used to acclerate the
	 * MST construction process.
	 */
	public ArrayList<Heuristic> getHeuristicList() { return heuristics; }
	
	/**
	 * Obtain the list of filters that were used to improve quality of
	 * clustering. Specifically the filters weed out ESTs that are known
	 * to interfere with clustering and deteriorate overall quality of 
	 * results.
	 * 
	 * @return The list of filters hat were used to improve quality of
	 * clustering.
	 */
	public ArrayList<Filter> getFilterList() { return filters; }
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent JobList node in the DOM tree.
	 * 
	 * @param jobList The DOM element corresponding to the "JobList"
	 * node that contains this entry.
	 */
	public final void marshall(Element jobList) {
		// Create a top-level entry for this "Job"
		Element job = DOMHelper.addElement(jobList, "Job", null);
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(job, "JobID", jobID);
		DOMHelper.addElement(job, "Description", DOMHelper.xmlEncode(description));
		DOMHelper.addElement(job, "ServerID", serverID);
		DOMHelper.addElement(job, "Path", path);
		DOMHelper.addElement(job, "Nodes", "" + nodes);
		DOMHelper.addElement(job, "CPUsPerNode", "" + cpusPerNode);
		DOMHelper.addElement(job, "Memory", "" + memory);
		DOMHelper.addElement(job, "MaxRunTime", "" + maxRunTime);
		
		if (startTimestamp != null) {
			DOMHelper.addElement(job, "StartTimestamp", startTimestamp.toString());
		}
		DOMHelper.addElement(job, "Status", status.toString().toLowerCase());
		DOMHelper.addElement(job, "LastUpdateTimestamp", lastUpdateTimestamp.toString());
		if (runtime != null) {
			DOMHelper.addElement(job, "RunTime", runtime.toString());
		}
		// Now marshall the information regarding heuristics.
		Element chain = DOMHelper.addElement(job, "HeuristicChain", null);
		for(Heuristic heuristic : heuristics) {
			heuristic.marshall(chain);
		}
		// Now marshall the information regarding filters.
		chain = DOMHelper.addElement(job, "FilterChain", null);
		for(Filter filter : filters) {
			filter.marshall(chain);
		}
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
		out.printf("%s<Job>\n", Indent); 
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "JobID", jobID);
		out.printf(STR_ELEMENT, "Description", description);
		out.printf(STR_ELEMENT, "ServerID", serverID);
		out.printf(STR_ELEMENT, "Path", path);
		out.printf(NUM_ELEMENT, "Nodes", nodes);
		out.printf(NUM_ELEMENT, "CPUsPerNode", cpusPerNode);
		out.printf(NUM_ELEMENT, "Memory", memory);
		out.printf(NUM_ELEMENT, "MaxRunTime", maxRunTime);
		
		if (startTimestamp != null) {
			out.printf(STR_ELEMENT, "StartTimestamp", startTimestamp.toString());
		}
		out.printf(STR_ELEMENT, "Status", status.toString().toLowerCase());
		out.printf(STR_ELEMENT, "LastUpdateTimestamp", lastUpdateTimestamp.toString());
		if (runtime != null) {
			out.printf(STR_ELEMENT, "RunTime", runtime.toString());
		}
		// Marhsall out the heuristic chain.
		out.printf("%s\t<HeuristicChain>\n", Indent);
		for(Heuristic heuristic : heuristics) {
			heuristic.marshall(out);
		}
		out.printf("%s\t</HeuristicChain>\n", Indent);
		
		// marhsall out the filter chain.
		out.printf("%s\t<FilterChain>\n", Indent);
		for(Filter filter : filters) {
			filter.marshall(out);
		}
		out.printf("%s\t</FilterChain>\n", Indent);
		
		// Close the job tag
		out.printf("%s</Job>\n", Indent);
	}

	/**
	 * Set the monitoring thread for this Job.
	 * 
	 * This method sets up the monitoring thread for this job. It
	 * also broadcasts an update event to all workspace listeners.
	 * This enables any GUI components to update their display.
	 * 
	 * @param monitor The monitoring thread for this job. If this
	 * parameter is null, then the job monitor thread is cleared.
	 */
	public synchronized void setMonitor(Thread monitor) {
		this.monitor = monitor;
		// Notify all the listeners about the status change.
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(this, WorkspaceEvent.Operation.UPDATE);
		Workspace.get().fireWorkspaceChanged(we);
	}
	
	/**
	 * Obtain the job monitoring thread (if any)
	 * 
	 * @return The job monitoring thread associated with this job. If
	 * a monitor thread is not set, then this method returns null.
	 */
	public Thread getMonitor() {
		return monitor;
	}
	
	/**
	 * A user defined description for this job. This description is set when
	 * a new job is scheduled and persisted in the work space configuration for
	 * future references.
	 */
	private String description;
	
	/**
	 * The path to the directory on the server (may it be local or remote)
	 * where the files corresponding to this job are stored.
	 */
	private String path;
	
	/**
	 * The number of nodes on a cluster that were requested for running this
	 * job.
	 */
	private int nodes;

	/**
	 * The number of CPUS per Nodes (on a cluster) that were requested for
	 * running this job.
	 */
	private int cpusPerNode;

	/**
	 * The peak memory to be reserved for this job when submitting via
	 * a PBS system. This value is in giga bytes.
	 */
	private int memory;

	/**
	 * The total wall clock run time to be reserved for this job when
	 * submitting via a PBS system. The value is in hours.
	 */
	private int maxRunTime;

	/**
	 * The time when this job was actually created and queued to run on 
	 * a server. This value does not have a bearing on the runtime of
	 * a job. This value is persisted for future reference on long running
	 * jobs.
	 */
	private XMLGregorianCalendar startTimestamp;
	
	/**
	 * The time stamp value of the last time the status of this job was
	 * updated. This value is periodically updated until the job completes
	 * running.
	 */
	private XMLGregorianCalendar lastUpdateTimestamp;
	
	/**
	 * The list of heuristics configured by the user for this job. These
	 * heuristics are used to accelerate the MST generation algorithm.
	 */
	private ArrayList<Heuristic> heuristics;
	
	/**
	 * The list of filters configured by the user for this job. These
	 * filters are used to weed out fragments that may interfer with
	 * clustering and reduce overall quality of clustering.
	 */
	private ArrayList<Filter> filters;
	
	/**
	 * A transient (not persisted) progress information about the job.
	 * The first entry is the number of ESTs that have been analyzed
	 * and the second entry is the total number of ESTs being analyzed.
	 * This information is useful for GUI's to display some progress
	 * information is desired. This value is periodically updated in 
	 * the run() method if the job is in QUEUED or RUNNING status.
	 */
	private transient int progressInfo[] = {-1, -1};
	
	/**
	 * A transient (not persistent) reference to the monitor thread
	 * for this job (if any). This value is set whenever a new 
	 * monitor thread is created for this job. The reference is reset
	 * when the monitor thread is stopped. The monitoring thread
	 * is actually implemented by org.peace_tools.core.JobMonitor class.
	 */
	private transient Thread monitor;
}
