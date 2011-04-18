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
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * A class to encapsulate information about a list of jobs that have already
 * been configured in this work space. This class is instantiated from the work
 * space class. This class is relatively straightforward in that it merely
 * contains a list of Job objects. In addition, it facilitates marshaling and
 * unmarshalling of job data.
 * 
 */
public class JobList {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * JobList entry. This method is typically  used to create a suitable
	 * JobList entry using data from a given DOM element. For each 
	 * Server entry this method uses the Job.create method to
	 * create suitable Job objects.
	 * 
	 * @param jobListNode The DOM element to be used for creating the 
	 * job list and to be used for creating the Job entries.
	 * 
	 * @return The newly created JobList entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static JobList create(Element jobListNode) throws Exception {
		// Create the jobList we are going to populate below.
		JobList srvrList = new JobList();
		// First extract the sequence counter information
		String seqCounter = DOMHelper.getStringValue(jobListNode, "SeqCounter");
		srvrList.seqCounter = Integer.parseInt(seqCounter);
		// Now check and process all the job elements under the jobList
		NodeList jobs = jobListNode.getChildNodes();
		for(int idx = 0; (idx < jobs.getLength()); idx++) {
			Node tmpNode = jobs.item(idx);
			// Ensure that this node is actually a job node
			if ((tmpNode.getNodeType() == Node.ELEMENT_NODE) &&
				(tmpNode instanceof Element)) {
				// OK, safe to type cast and try to parse..
				Element jobNode = (Element) tmpNode;
				final String nodeName = jobNode.getNodeName();
				if ("ClusteringJob".equals(nodeName)) { 
					Job j = ClusteringJob.create(jobNode);
					// Add the valid job node to the list.
					srvrList.jobs.add(j);
				} else if ("EASTJob".equals(nodeName)) {
					Job j = EASTJob.create(jobNode);
					srvrList.jobs.add(j);
				}
			}
		}
		// Return the newly created list.
		return srvrList;
	}
	
	/**
	 * The default constructor. It merely initializes all the instance
	 * variables to their default initial value.
	 */
	public JobList() {
		jobs       = new ArrayList<Job>();
		seqCounter = 1;
	}

	/**
	 * Reserves the next job ID for use in a job object. This method is
	 * used only when a new job object is added by the user. job objects
	 * loaded from a file have their own job IDs.
	 * 
	 * @return A new, unique (within the work space) job ID for use with
	 * a new job object instance.
	 */
	public String reserveJobID() {
		// Use current sequence counter value to create ID
		String id = "job" + seqCounter;
		// Setup sequence counter for next job.
		seqCounter++;
		// Return the ID for use by the caller.
		return id;
	}
	
	/**
	 * Obtain reference to a job object, given its work space wide unique
	 * job ID. This method does a linear search on the jobList to
	 * locate the job entry (if you have 1000s of job entries in one
	 * workspace then you are missing the concept of a workspace). 
	 * 
	 * @param jobID The generated work space wide unique job ID. 
	 * Note that checks are case sensitive. 
	 * 
	 * @return Reference to a job object whose jobID is the same as the
	 * supplied jobID. If a matching object was not found then this method
	 * returns null.
	 */
	public Job getjob(String jobID) {
		for(Job j: jobs) {
			if (j.getJobID().equals(jobID)) {
				return j;
			}
		}
		// No matching entry found.
		return null;
	}
	
	/**
	 * Method to add a new job entry to the job list.
	 * 
	 * This method adds the specified job to the job list and 
	 * fires a WorkspaceEvent indicating the addition of the new
	 * entry to all listeners.
	 * 
	 * @param job The new job entry to be added to the job list.
	 */
	public synchronized void add(Job job) {
		// Add the entry
		jobs.add(job);
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(job, WorkspaceEvent.Operation.INSERT);
		Workspace.get().fireWorkspaceChanged(we);
	}

	/**
	 * Method to remove a job entry.
	 * 
	 * This method removes the given job entry from the job list
	 * and fires a WorkspaceEvent indicating the removal of the
	 * job entry to all listeners.
	 * 
	 * <p><b>Note:</b>  This method does not delete any of the files associated
	 * with the job. It simply removes the job entry from the work
	 * space.</p>
	 * 
	 * @param job The job entry to be removed from the job list.
	 */
	public synchronized void remove(Job job) {
		// Add the entry
		if (jobs.remove(job)) {
			// Fire notification to listeners to update GUIs
			WorkspaceEvent we = new WorkspaceEvent(job, WorkspaceEvent.Operation.DELETE);
			Workspace.get().fireWorkspaceChanged(we);
		}
	}
	/**
	 * Obtain the list of job objects that are currently available
	 * in a work space.
	 * 
	 * @return The list of job objects that are currently available
	 * in this work space.
	 */
	public ArrayList<Job> getJobs() { return jobs; }

	/**
	 * Obtain the sequence counter used to generate unique job IDs.
	 * 
	 * @return The sequence counter for this work space.
	 */
	public long getSeqCounter() { return seqCounter; }
	
	/**
	 * Method to marshal the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent Workspace node in the DOM tree.
	 * 
	 * @param workspace The DOM element corresponding to the Workspace
	 * node that contains this entry.
	 */
	public final void marshall(Element workspace) {
		// Create a top-level job list entry for this class
		Element srvrList = DOMHelper.addElement(workspace, "JobList", null);
		// Add the sequence counter information to job list
		DOMHelper.addElement(srvrList, "SeqCounter", "" + seqCounter);
		// Add sub-elements for each job in our job list
		for (Job j : jobs) {
			j.marshall(srvrList);
		}
	}
	
	/**
	 * Method to marshal the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. This method provides
	 * better control on the XML formatting to generate a more readable
	 * XML output when compared to the DOM tree.
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t";
		// Create a top-level job list entry for this class
		out.printf("%s<JobList>\n", Indent);
		// Add the sequence counter information to job list
		out.printf("%s\t<%2$s>%3$s</%2$s>\n", Indent, "SeqCounter", seqCounter);
		// Add sub-elements for each job in our job list
		for (Job j : jobs) {
			j.marshall(out);
		}
		// Close the job list element.
		out.printf("%s</JobList>\n\n", Indent);
	}
	
	/**
	 * The list of job objects that have been configured and added
	 * to this list. This list is populated either via the create 
	 * method or new entries are added via the addjob() method
	 * in this class. 
	 */
	private ArrayList<Job> jobs;
	
	/**
	 * Sequence counter that is maintained on a per-work space basis to
	 * generate unique/valid IDs for each new job entry added to this
	 * list.
	 */
	private long seqCounter;
}
