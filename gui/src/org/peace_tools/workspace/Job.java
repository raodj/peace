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
import java.util.GregorianCalendar;

import javax.xml.datatype.DatatypeConfigurationException;
import javax.xml.datatype.DatatypeFactory;
import javax.xml.datatype.XMLGregorianCalendar;

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.generic.ProgrammerLog;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * <p>This class is a common base class that includes some of the
 * data members and functionality that is shared between various
 * types of derived job classes. A similar extension strategy
 * is used to define the XML schema for jobs. However, unlike in
 * XML definition, this class cannot be directly instantiated.
 * Instead, one of the derived classes must be instantiated
 * and used. The instantiation and use of derived classes is
 * performed by the Workspace element when it processes an 
 * XML schema.</p>
 * 
 * <p>This class corresponds to a generic Job  element in a PEACE
 * work space XML data. This class encapsulates the core information
 * associated with a Job. This class serves merely as a read-only 
 * type class that is created when a new job is run. The data is 
 * persisted in  the XML work space for future reference so that
 * users can determine jobs that were scheduled and check on status
 * of long-running jobs.</p>
 * 
 * <p>In addition to encapsulating the data this class also provides
 * convenient interfaces for marshaling and un-marshaling XML data compatible
 * with the PEACE GUI configuration XML file format. </p>
 */
public abstract class Job extends JobBase {
	/**
	 * <p>Helper method to utilize data from a DOM tree to populate various
	 * common information corresponding to a given job entry.
	 * This method is typically  used to populate the job object with
	 * appropriate data when loading a Work space into the GUI.</p>
	 * 
	 * <p><b>Note:</b>This method only unmarshals non-final fields. Some
	 * of the elements are immutable when a job entry is created.</p>
	 * 
	 * @param jobNode The DOM element to be used for creating the Job
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created job entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	protected void unmarshal(Element jobNode) throws Exception {
		// First extract the generic job type information
		// Extract description (if available)
		String desc = DOMHelper.getStringValue(jobNode, "Description", true);
		description = (desc != null) ? desc : "";
		// Extract the various server and common job information.
		path        = DOMHelper.getStringValue(jobNode, "Path");
		nodes       = DOMHelper.getIntValue(jobNode, "Nodes");
		cpusPerNode = DOMHelper.getIntValue(jobNode, "CPUsPerNode");
		memory      = DOMHelper.getIntValue(jobNode, "Memory");
		maxRunTime  = DOMHelper.getIntValue(jobNode, "MaxRunTime");
		// Extract the time information and convert it to appropriate
		// usable objects.
		String lastTime       = DOMHelper.getStringValue(jobNode, "LastUpdateTimestamp");
		DatatypeFactory codec = DatatypeFactory.newInstance();
		lastUpdateTimestamp   = codec.newXMLGregorianCalendar(lastTime);
		// Extract status and convert it to an enumeration object.
		String statStr        = DOMHelper.getStringValue(jobNode, "Status");
		statStr               = statStr.toUpperCase();
		status                = JobBase.JobStatusType.valueOf(statStr);

		// Extract common optional information for a job
		if (DOMHelper.hasElement(jobNode, "StartTimestamp")) {
			String startTime = DOMHelper.getStringValue(jobNode, "StartTimestamp");
		    startTimestamp   = codec.newXMLGregorianCalendar(startTime);
		}
		if (DOMHelper.hasElement(jobNode, "Runtime")) {
			String runtimeStr = DOMHelper.getStringValue(jobNode, "RunTime");
			runtime           = codec.newDuration(runtimeStr);
		}
		// Load ID of a job that we are dependent upon (if any)
		if (DOMHelper.hasElement(jobNode, "PrevJobID")) {
			String prevJobID  = DOMHelper.getStringValue(jobNode, "PrevJobID");
			setPreviousJobID(prevJobID);
		}
	}

	/**
	 * Helper method to un-marshal the parameters (if any)
	 * into an in-memory array list.
	 * 
	 * This is a helper method that is invoked from the the create() method to 
	 * un-marshal the XML DOM tree corresponding to the parameters into
	 * suitable in-memory classes for further processing.
	 * 
	 * @param jobNode The top-level job node from where the parameter elements
	 * are to be extracted and processed. 
	 * 
	 * @return An array list containing the list of parameter objects in the 
	 * job entry. If no parameters were found, then this method returns
	 * null. 
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	protected static ArrayList<Param> parseParameters(Element jobNode) throws Exception {
		// Parse out the parameters (if any) for this job.
		NodeList paramList = jobNode.getElementsByTagName("Param");
		if ((paramList == null) || (paramList.getLength() == 0)) {
			return null;
		}
		ArrayList<Param> parameters = new ArrayList<Param>(); 
		for(int idx = 0; (idx < paramList.getLength()); idx++) {
			Element paramNode = (Element) paramList.item(idx);
			// Skip Param nodes that are not direct children of the jobNode
			if (paramNode.getParentNode() != jobNode) {
				// This is not a direct parameter node for this job. This could
				// be parameter element present inside filter chain or some other element.
				continue;
			}
			Param   entry     = Param.create(paramNode);
			// Add the parameter information to the parameters list.
			parameters.add(entry);
		}
		// Return the parsed-in heuristic chain for further use.
		return parameters;
	}
	
	/**
	 * Constructor to create a minimally populated Job object with just
	 * the fixed value fields initialized to specific values.
	 * 
	 * This constructor is typically used to create a temporary job
	 * object into which additional values are going to be loaded
	 * from a XML work space. This constructor is typically used by
	 * derived job classes.
	 * 
	 * @param type The type of job entry being created.
	 * 
	 * @param jobID The workspace-wide unique, generated job ID for this job.
	 * The jobID is generated via a call to JobList.reserveJobID() method.
	 * This is a valid string (typically in the form job####)
	 * 
	 * @param serverID The ID of the server on which this job is running.
	 * This is a cross reference ID of a Server configured in this workspace.
	*/
	protected Job(JobBase.JobType type, String jobID, String serverID) {
		super(type, jobID, serverID);
		this.description         = "";
		this.path                = null;
		this.nodes               = -1;
		this.cpusPerNode         = -1;
		// Set non-final fields to null for now.
		this.startTimestamp      = null;
		this.lastUpdateTimestamp = null;
		this.runtime             = null;
		this.status              = JobStatusType.STARTING;
		this.memory              = -1;
		this.maxRunTime          = -1;
		this.lastUpdateTimestamp = null;
	}
	
	/**
	 * Constructor to create a Job object with the fixed value fields
	 * initialized to specific values.
	 * 
	 * @param type The type of job entry being created.
	 * 
	 * @param jobID The workspace-wide unique, generated job ID for this job.
	 * The jobID is generated via a call to JobList.reserveJobID() method.
	 * This is a valid string (typically in the form job####)
	 * 
	 * @param serverID The ID of the server on which this job is running.
	 * This is a cross reference ID of a Server configured in this workspace.
	 * 
	 * @param description A user-supplied description for this job entry. This
	 * maybe an empty string (but cannot be null).
	 * 
	 * @param path The directory on the server where the data for this job 
	 * is stored on the local machine on which PEACE is running.
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
	 */
	public Job(JobBase.JobType type, String jobID, String serverID, 
			String description, String path, int nodes, int cpusPerNode,
			int memory, int maxRunTime) {
		super(type, jobID, serverID);
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
		// Set default last update time stamp
		setLastUpdateTime();
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
		// Fire notification to listeners to update GUIs. In addition
		// MainFrame will update dependent job status and start a monitoring
		// thread for the dependent job. 
		WorkspaceEvent we = new WorkspaceEvent(this, WorkspaceEvent.Operation.UPDATE);
		Workspace.get().fireWorkspaceChanged(we);
		// Update status of any job summaries in the data sets
		GeneratedFileList gfl = Workspace.get().getGFL(jobID);
		if (gfl != null) {
			gfl.updateJobSummary(this);
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
	 * Obtain the time stamp indicating when the job actually started running on
	 * the server.
	 * 
	 * @return The time stamp when the job started running. If the job has not
	 * yet started running, then this method returns null.
	 */
	public XMLGregorianCalendar getStartTimestamp() {
		return startTimestamp;
	}


	/**
	 * Set the progress status for this job. This method generates 
	 * notifications to all registered workspace listeners.
	 * 
	 * @param estsAnalyzed The number of ESTs that have been analyzed.
	 * @param totalEstCount The total number of ESTs to be analyzed.
	 */
	public final void setProgress(int estsAnalyzed, int totalEstCount) {
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
	public final int[] getProgressInfo() { return progressInfo; }
	
	@Override
	public String toString() {
		return jobID;
	}
	
	/**
	 * Return the information in the form of a complete general purpose
	 * command line.
	 * 
	 * This method can be used to obtain the information needed to
	 * run this job with all the necessary command line parameters.
	 * 
	 * <p><b>Note:</b>The command line must not include the executable
	 * path but just the executable name. The executable path is added
	 * at a later time when the job is submitted to the server based on
	 * the installation path on the specific server.</p>
	 * 
	 * @return Return the information as a command line parameter.
	 */
	public abstract String toCmdLine();
	
	/**
	 * Method to marshal the common data stored in this object 
	 * to become part of a DOM tree element passed in. This 
	 * method assumes that the element passed in corresponds to
	 * the parent JobList node in the DOM tree.
	 * 
	 * @param job The DOM element corresponding to the "Job"
	 * node that contains this entry.
	 */
	public void marshall(Element job) {
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(job, "Type", type.toString().toLowerCase());
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
	 * 
	 * @param indentPrefix The extra indentation to be done to make the
	 * output look nice. If no additional indentation is needed then
	 * an empty string ("") must be passed in.
	 */
	public void marshall(PrintWriter out, final String indentPrefix) {
		final String Indent = indentPrefix + "\t\t";
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		final String NUM_ELEMENT = Indent + "\t" + "<%1$s>%2$d</%1$s>\n";
		final String xmlDesc     = DOMHelper.xmlEncode(description);
		
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "Type", type.toString().toLowerCase());
		out.printf(STR_ELEMENT, "JobID", jobID);
		out.printf(STR_ELEMENT, "Description", xmlDesc);
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
		if (prevJobID != null) {
			out.printf(STR_ELEMENT, "PrevJobID", prevJobID);
		}
	}

	/**
	 * Helper method to marshal the parameters associated with a job
	 * to a given DOM element.
	 * 
	 * This method is used by derived classes to conveniently marshal
	 * any parameters they may have to a given DOM node.
	 * 
	 * @param parameters The list of parameters to be marshaled to a
	 * given DOM element. If this parameter is null, then this method
	 * does not perform any operation.
	 * 
	 * @param job The DOM element representing a job to which the 
	 * parameters are to be added.
	 */
	protected void marshallParameters(ArrayList<Param> parameters, Element job) {
		if (parameters != null) {
			for(Param p : parameters) {
				p.marshall(job);
			}
		}
	}
	
	/**
	 * Helper method to marshal the parameters associated with a job
	 * to a given output stream.
	 * 
	 * This method is used by derived classes to conveniently marshal
	 * any parameters they may have to a given output stream.
	 * 
	 * @param parameters The list of parameters to be marshaled to a
	 * given DOM element. If this parameter is null, then this method
	 * does not perform any operation.
	 * 
	 * @param out The stream to which the XML must be serialized.	 
	 */
	protected void marshallParameters(ArrayList<Param> parameters, PrintWriter out) {
		if (parameters != null) {
			for(Param p : parameters) {
				p.marshall(out);
			}
		}
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
	public void summarize(SummaryWriter sw) {
		// Summarize information about the server & config
		final Server srvr = Workspace.get().getServerList().getServer(serverID);
		sw.addSection("Server-Job Summary");
		sw.addSummary("Server", srvr.getName(), srvr.getDescription());
		sw.addSummary("Nodes", "" + nodes, null);
		sw.addSummary("CPUs per Node", "" + cpusPerNode, null);
		sw.addSummary("Time requested", "" + maxRunTime + " hours", 
			"The job should not take more CPU time than this");
		sw.addSummary("Memory per Node", "" + memory + " MB",
				"Peak memor per node");
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
	public final synchronized void setMonitor(Thread monitor) {
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
	public final Thread getMonitor() {
		return monitor;
	}
	
	/**
	 * Command line for PEACE clustering engine to configure heuristics,
	 * if applicable.
	 * 
	 * <p>This method can be used to obtain the heuristic information 
	 * in the form of command line parameters that can be readily passed
	 * to the PEACE clustering engine. The command line parameters are 
	 * used to configure the clustering engine to suit the configuration
	 * setup by the user for this job.</p>
	 * 
	 * <p><b>Note:<b> The return value of this method may be an empty
	 * string if heuristics are not applicable for this job. Derived
	 * classes override this method to return suitable heuristic
	 * arguments if applicable.</p> 
	 * 
	 * @return The command line (to correspondingly setup the heuristics) to
	 * be passed to the PEACE clustering engine. If no heuristics are 
	 * defined then this method returns an empty string.
	 */
	public String getHeuristicsCmdLine() {
		return "";
	}
	
	/**
	 * General purpose command line with generic set of parameters
	 * setup for this job.
	 *  
	 * <p>This method can be used to obtain the set of generic parameters
	 * set for this job in the form of command line parameters. </p>
	 * 
	 * <p><b>Note:<b> The return value of this method may be an empty
	 * string if parameters are not applicable for this job. Derived
	 * classes override this method to return suitable heuristic
	 * arguments if applicable.</p> 
	 * 
	 * @return The command line corresponding to the parameters set for
	 * this job. If parameters are not applicable, then this method 
	 * returns an empty string.
	 */
	public String getParametersCmdLine() {
		return "";
	}

	/**
	 * Helper method to convert parameters to a suitable command line.
	 * 
	 * This is a helper method that can be used by the derived classes
	 * in their implementation of the {@link #getParametersCmdLine()}
	 * method to convert parameters (if any) to a command line.
	 * 
	 * @param parameters The list of parameters to be converted to a
	 * command line. This parameter can be null (in which case an
	 * empty string is returned).
	 * 
	 * @return A partial command line representing the set of 
	 * parameters in the given list.
	 */
	protected String toCmdLine(ArrayList<Param> parameters) {
		if (parameters == null) {
			return ""; // no filters at all
		}
		// Next covert parameter information to command line parameters.
		String cmdLine = "";
		for(int i = 0; (i < parameters.size()); i++) {
			Param p = parameters.get(i);
			cmdLine += (i > 0 ? " " : "") + p;
		}
		// Return the cmd line.
		return cmdLine;
	}
	
	/**
	 * Command line for PEACE clustering engine to configure filters,
	 * if applicable for this type of job.
	 * 
	 * This method can be used to obtain the filters information in the form
	 * of command line parameters that can be readily passed to the PEACE
	 * clustering engine. The command line parameters are used to configure the
	 * filters used by the clustering engine to mirror the configuration setup 
	 * by the user for this job.
	 * 
	 * <p><b>Note:<b> The return value of this method may be an empty
	 * string if heuristics are not applicable for this job. Derived
	 * classes override this method to return suitable heuristic
	 * arguments if applicable.</p>
	 * 
	 * @return The command line (to correspondingly setup the filters) to
	 * be passed to the PEACE clustering engine. If no filters are 
	 * defined or if filters are not applicable then this method 
	 * returns an empty string.
	 */
	public String getFiltersCmdLine() {
		return "";
	}
	
    /**
     * Helper method to compute the full path to a given data set file
     * on the server.
     * 
     * PEACE may have been installed at different paths on different
     * servers/machines. This method helps to add the necessary
     * path-prefix to a given data set file so that it can be 
     * appropriately located on the remote machine.
     *
     * @param ds The data set containing the source cDNA fragment file
     * to be mapped to an absolute path on the remote server.
     * 
     * @return The full path (based on install path) to the EST data
     * file on the server set for this job.
     */
	protected String getServerESTFile(DataSet ds) {
        String localFileName = ds.getPath();
        File localFile = new File(localFileName);
        Server srvr = Workspace.get().getServerList().getServer(getServerID());
        // Construct server-specific EST file location.
        String serverESTFileName = srvr.getServerPath(localFile);
        return serverESTFileName;
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
