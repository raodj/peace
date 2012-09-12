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

import org.peace_tools.core.JobMonitor;
import org.peace_tools.core.MainFrame;
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
	 * This is an inner class that is used to a list that can
	 * consist of either Job entries or JobList entries (but not both
	 * in the same object). This class is used to maintain the
	 * entries in the same order that the user organizes the jobs.
	 */
	public class JobOrListWrapper {
		/**
		 * Convenience constructor to create a wrapper class to contain
		 * a job entry in it. 
		 * 
		 * @param job The job entry to be encapsulated in this object.
		 * This parameter cannot be null.
		 */
		public JobOrListWrapper(Job job) {
			this.job        = job;
			this.jobSubList = null;
		}

		/**
		 * Convenience constructor to create a wrapper class to contain
		 * a sub-list entry in it.
		 * 
		 * @param subList The job sub-list entry to be encapsulated 
		 * in this object. This parameter cannot be null.
		 */
		public JobOrListWrapper(JobList jobSubList) {
			this.job        = null;
			this.jobSubList = jobSubList;
		}

		/**
		 * Method to determine if the wrapper object contain a job or
		 * a sub-list.
		 * 
		 * @return This method returns true if the wrapper object 
		 * contains a sub-list. If the wrapper object contains a 
		 * job then this method returns false.
		 */
		public boolean isSubList() {
			return (jobSubList != null);
		}

		/**
		 * Obtain the job entry encapsulated by this wrapper object.
		 * This method will return a valid object when the
		 * {@link #isSubList()} method returns false.
		 * 
		 * @return The job object (if any) that is encapsulated by 
		 * this object. The return value will be null if the wrapper
		 * is object is encapsulating a sub-list instead of a job.
		 */
		public Job getJob() {
			return job;
		}

		/**
		 * Obtain the job sub-list entry encapsulated by this wrapper 
		 * object.
		 * 
		 * This method will return a valid object when the
		 * {@link #isSubList()} method returns true.
		 * 
		 * @return The job sub-list object (if any) that is encapsulated
		 * by this object. The return value will be null if the wrapper
		 * is object is encapsulating a job object instead of a 
		 * job sub-list.
		 */
		public JobList getSubList() {
			return jobSubList;
		}

		@Override
		public String toString() {
			return (isSubList() ? jobSubList.getName() : job.jobID);
		}

		/**
		 * The job object (if any) that is encapsulated by this object.
		 * This object can be null if the wrapper is object is
		 * encapsulating a sub-list instead of a job.
		 */
		private final Job job;

		/**
		 * The job object (if any) that is encapsulated by this object.
		 * This object can be null if the wrapper is object is
		 * encapsulating a sub-list instead of a job.
		 */
		private final JobList jobSubList;
	}

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
		JobList jobList = new JobList();
		// First extract the sequence counter information
		if (DOMHelper.hasElement(jobListNode, "SeqCounter")) {
			String seqCounter = DOMHelper.getStringValue(jobListNode, "SeqCounter");
			jobList.seqCounter = Integer.parseInt(seqCounter);
		} else {
			jobList.seqCounter = -1;
		}
		// Extract the name and description for the job list
		jobList.name        = DOMHelper.getStringValue(jobListNode, "Name");
		jobList.description = DOMHelper.getStringValue(jobListNode, "Description");
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
					jobList.jobs.add(jobList.new JobOrListWrapper(j));
				} else if ("EASTJob".equals(nodeName)) {
					Job j = EASTJob.create(jobNode);
					jobList.jobs.add(jobList.new JobOrListWrapper(j));
				} else if ("DECAGONJob".equals(nodeName)) {
					Job j = DECAGONJob.create(jobNode);
					jobList.jobs.add(jobList.new JobOrListWrapper(j));
				} else if ("JobList".equals(nodeName)) {
					JobList subList = JobList.create(jobNode);
					jobList.jobs.add(jobList.new JobOrListWrapper(subList));
				}
			}
		}
		// Return the newly created list.
		return jobList;
	}

	/**
	 * The default constructor. It merely initializes all the instance
	 * variables to their default initial value.
	 */
	public JobList() {
		jobs        = new ArrayList<JobOrListWrapper>();
		seqCounter  = 1;
		name        = "JobList";
		description = "The main job list";
	}

	/**
	 * The default constructor to create a job list with given name and
	 * description. 
	 * 
	 * @param name The name to be set for this job list. This value 
	 * cannot be null.
	 * 
	 * @param description The description to be set for this job list.
	 */
	public JobList(String name, String description) {
		this.jobs        = new ArrayList<JobOrListWrapper>();
		this.seqCounter  = 1;
		this.name        = name;
		this.description = description;
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
	 * job ID. This method does a recursive search on the jobList to
	 * locate the job entry (if you have 1000s of job  entries in 100s 
	 * of sub-lists in one workspace then you are missing the concept 
	 * of a workspace). 
	 * 
	 * @param jobID The generated work space wide unique job ID. 
	 * Note that checks are case sensitive. 
	 * 
	 * @return Reference to a job object whose jobID is the same as the
	 * supplied jobID. If a matching object was not found then this method
	 * returns null.
	 */
	public Job getjob(String jobID) {
		Job retVal = null;
		for(JobOrListWrapper wrapper: jobs) {
			if ((wrapper.isSubList()) &&
					((retVal = wrapper.getSubList().getjob(jobID)) != null)) {
				// Found matching entry in a sub-list.
				break;
			} else if ((!wrapper.isSubList()) &&
					(wrapper.getJob().getJobID().equals(jobID))) {
				retVal = wrapper.getJob();
				break;
			}
		}
		// Return job object back to caller.
		return retVal;
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
		// Add the entry suitably wrapped.
		JobOrListWrapper entry = new JobOrListWrapper(job);
		jobs.add(entry);
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(job, WorkspaceEvent.Operation.INSERT);
		Workspace.get().fireWorkspaceChanged(we);
	}

	/**
	 * Method to insert a new job entry to the job list while ensuring
	 * the job entry is prior to any other job entries that depend on it.
	 * 
	 * This method inserts the specified job to the job list and 
	 * fires a WorkspaceEvent indicating the addition of the new
	 * entry to all listeners. The insertion is done such that
	 * the order of dependence is reflected in the list to help with
	 * user understanding of dependence. To accomplish this objective,
	 * the job list performs a recursive depth-first search for the
	 * dependent job prior to inserting it.
	 * 
	 * @param job The new job entry to be added to the job list.
	 */
	public synchronized void insert(Job job) {
		insert(job, true);
	}

	/**
	 * Get the list of jobs that are marked as being dependent on a given 
	 * job ID.
	 * 
	 * This method provides a convenient mechanism to recursively search
	 * the job-hierarchy to locate various jobs that are dependent (typically,
	 * the dependent jobs are waiting for an earlier job to complete) on
	 * a given job ID. All the dependent jobs (if found) are added to the
	 * given array list.
	 * 
	 * @param jobID The job ID for which dependent jobs are to be determined.
	 * This value cannot be null.
	 * 
	 * @param dependentJobs An array list to which all dependent jobs are to 
	 * be added. Entries are only added to the array list and existing values
	 * are not modified.
	 */
	public void getDependentJobs(final String jobID, ArrayList<Job> dependentJobs) {
		for(JobOrListWrapper wrapper: jobs) {
			if (wrapper.isSubList()) {
				wrapper.getSubList().getDependentJobs(jobID, dependentJobs);
			} else {
				final String prevJobID = wrapper.getJob().getPreviousJobID();
				if ((prevJobID != null) && prevJobID.equals(jobID)) {
					dependentJobs.add(wrapper.getJob());
				}
			}
		}
	}

	/**
	 * Helper method to insert a new job entry to the job list 
	 * while ensuring the job entry is prior to any other job 
	 * entries that depend on it.
	 * 
	 * This method is invoked from the {@link #insert(Job)} method
	 * This method inserts the specified job to the job list and 
	 * fires a WorkspaceEvent indicating the addition of the new
	 * entry to all listeners. The insertion is done such that
	 * the order of dependence is reflected in the list to help with
	 * user understanding of dependence. To accomplish this objective,
	 * the job list performs a recursive depth-first search for the
	 * dependent job prior to inserting it.
	 * 
	 * @param job The new job entry to be added to the job list.
	 * 
	 * @param addAtEnd If this flag is true, then this method adds
	 * the job to the end of its job list. This flag is set to true
	 * for top-level calls, while recursive descends into sublists
	 * has this flag set to false so that entries are added only at 
	 * the level.
	 * 
	 * @return This method returns true if the job was successfully inserted.
	 * Otherwise this method returns false.
	 */
	public synchronized boolean insert(Job job, final boolean addAtEnd) {
		// Find first job entry (if any) that is dependent on this job.
		int index;
		for(index = 0; (index < jobs.size()); index++) {
			final JobOrListWrapper wrapper = jobs.get(index);
			if (wrapper.isSubList()) {
				// Recursively descend into this list to find dependent job
				if (wrapper.getSubList().insert(job, false)) {
					// The dependent job was found and the given job has
					// been successfully inserted
					return true;
				}
			} else {
				final String currDepJobID = wrapper.getJob().getPreviousJobID();
				if ((currDepJobID != null) && (currDepJobID.equals(job.getJobID()))) {
					// Found the first job that is dependent on this job. So insert
					// new job before the current one.
					break;
				}
			}
		}

		if ((index < jobs.size()) || addAtEnd) {
			JobOrListWrapper entry = new JobOrListWrapper(job);
			jobs.add(index, entry);
			// Fire notification to listeners to update GUIs
			WorkspaceEvent we = new WorkspaceEvent(job, WorkspaceEvent.Operation.INSERT);
			Workspace.get().fireWorkspaceChanged(we);
			// Entry added.
			return true;
		}
		// Entry not added
		return false;
	}

	/**
	 * Helper method to insert a new job entry to the job list at
	 * a given index position.
	 * 
	 * This method inserts the specified job to this job list and 
	 * fires a WorkspaceEvent indicating the addition of the new
	 * entry to all listeners. The job entry is inserted at the
	 * given index position.
	 * 
	 * @param job The new job entry to be added to the job list.
	 * 
	 * @param position The zero-based index position where the 
	 * job entry is to be inserted. If the position is -1, then the
	 * entry is appended at the end of the job list.
	 * 
	 * @return This method returns true if the job was successfully inserted.
	 * Otherwise this method returns false.
	 */
	public synchronized JobOrListWrapper insert(Job job, final int position) {
		JobOrListWrapper entry = new JobOrListWrapper(job);
		if ((position < 0) || (position > jobs.size())) {
			// Append at the end
			jobs.add(entry);
		} else {
			// Insert in the middle.
			jobs.add(position, entry);
		}
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(job, WorkspaceEvent.Operation.INSERT);
		Workspace.get().fireWorkspaceChanged(we);
		// Fire notifications regarding change this job list.
		we = new WorkspaceEvent(this, WorkspaceEvent.Operation.INSERT, position, entry);
		Workspace.get().fireWorkspaceChanged(we);
		return entry;
	}

	/**
	 * Helper method to insert a new sublist to the job list at
	 * a given index position.
	 * 
	 * This method inserts the specified job sublist to this job 
	 * list and fires a WorkspaceEvent indicating the addition of
	 * the new entry to all listeners. The job entry is inserted 
	 * at the given index position.
	 * 
	 * @param jobList The job sublist entry to be added to the job list.
	 * 
	 * @param position The zero-based index position where the 
	 * job entry is to be inserted. If the position is -1, then the
	 * entry is appended at the end of the job list.
	 * 
	 * @return This method returns true if the job was successfully inserted.
	 * Otherwise this method returns false.
	 */
	public synchronized JobOrListWrapper insert(JobList jobList, final int position) {
		JobOrListWrapper entry = new JobOrListWrapper(jobList);
		if ((position < 0) || (position > jobs.size())) {
			// Append at the end
			jobs.add(entry);
		} else {
			// Insert in the middle.
			jobs.add(position, entry);
		}		
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(this, WorkspaceEvent.Operation.INSERT, position, entry);
		Workspace.get().fireWorkspaceChanged(we);
		return entry;
	}

	/**
	 * Method to remove an entry (could be job or job-list entry).
	 * 
	 * This method removes the given entry from this job list
	 * and fires a WorkspaceEvent indicating the removal of the
	 * job entry to all listeners.
	 * 
	 * <p><b>Note:</b>  This method does not delete any of the files 
	 * associated with a job or a job list. It simply removes the 
	 * entry from this list. In addition, this method does not recursively
	 * descend sublists to locate the given entry. It only operates
	 * within the scope of this job list.</p>
	 * 
	 * @param entry The job entry to be removed from the job list.
	 * 
	 * @return This method returns true if the job entry was successfully
	 * removed. Otherwise this method returns false.
	 */
	public synchronized boolean remove(JobOrListWrapper entry) {
		final int position = jobs.indexOf(entry);
		// Remove the entry
		if (jobs.remove(entry)) {
			// Fire notification to listeners to update GUIs
			if (!entry.isSubList()) {
				// For jobs we fire two events.
				WorkspaceEvent we = new WorkspaceEvent(entry.getJob(), WorkspaceEvent.Operation.DELETE);
				Workspace.get().fireWorkspaceChanged(we);
			}
			// Fire change in this list
			WorkspaceEvent we = new WorkspaceEvent(this, WorkspaceEvent.Operation.DELETE, position, entry);
			Workspace.get().fireWorkspaceChanged(we);
			// Job entry successfully removed.
			return true;
		}
		// Job entry was not removed.
		return false;
	}

	/**
	 * Method to remove an entry (could be job or job-list entry).
	 * 
	 * This method removes the given entry from this job list
	 * and fires a WorkspaceEvent indicating the removal of the
	 * job entry to all listeners.
	 * 
	 * <p><b>Note:</b>  This method does not delete any of the files 
	 * associated with a job or a job list. It simply removes the 
	 * entry from this list. In addition, this method does not recursively
	 * descend sublists to locate the given entry. It only operates
	 * within the scope of this job list.</p>
	 * 
	 * @param entry The job entry to be removed from the job list.
	 * 
	 * @return This method returns true if the job entry was successfully
	 * removed. Otherwise this method returns false.
	 */
	public synchronized boolean remove(Job job) {
		// Recursively search the set of jobs until the matching job entry
		// is found.
		for(int index = 0; (index < jobs.size()); index++) {
			JobOrListWrapper wrapper = jobs.get(index);
			if (wrapper.isSubList()) {
				if (wrapper.getSubList().remove(job)) {
					// A sublist had the job entry and it was removed.
					// No need for any more operations
					return true;
				}
			} else if (wrapper.getJob() == job) {				
				if (jobs.remove(wrapper)) {
					// Fire notification to listeners to update GUIs
					WorkspaceEvent we =	new WorkspaceEvent(job, 
							WorkspaceEvent.Operation.DELETE);
					Workspace.get().fireWorkspaceChanged(we);
					// Fire change in this list
					we = new WorkspaceEvent(this, 
							WorkspaceEvent.Operation.DELETE, index, wrapper);
					Workspace.get().fireWorkspaceChanged(we);
					// Job entry successfully removed.
					return true;
				}					
			}

		}
		// Job entry was not removed.
		return false;
	}

	/**
	 * Recursively search this list and sublists to locate the JobList 
	 * that contains a given entry.
	 * 
	 * This method can be used to locate a job list that contains a given
	 * entry. This method searches in a depth-first manner to locate 
	 * the job list that contains the given entry. 
	 *  
	 * @param entry The entry to search for.
	 *  
	 * @return The job list containing the entry, if the entry is found.
	 * Otherwise this method returns null.
	 */
	public JobList find(JobOrListWrapper entry) {
		JobList retVal = null;
		for(JobOrListWrapper wrapper: jobs) {
			if (wrapper == entry) {
				retVal = this;
				break;
			} else if (wrapper.isSubList()) {
				if ((retVal = wrapper.getSubList().find(entry)) != null) {
					break;
				}
			}
		}
		return retVal;
	}

	/**
	 * Method to create threads for all jobs whose status is 
	 * to be monitored.
	 * 
	 * This method is typically invoked from {@link MainFrame#createJobThreads()}
	 * method after the main frame is created. This method performs the task of 
	 * checking job status and creating a job thread if a job needs to be
	 * monitored. This method recursively descends the job list hierarchy
	 * to create job monitors.
	 * 
	 * <p><b>Note:</b>  Invoking this method twice will cause unnecessary
	 * monitoring threads to start up. So avoid duplicate calls.</p>
	 */
	public void createJobThreads(final MainFrame mf) {
		for(JobOrListWrapper wrapper: jobs) {
			if (wrapper.isSubList()) {
				wrapper.getSubList().createJobThreads(mf);
			} else {
				final Job job = wrapper.getJob();
				if (!job.isDone() && !job.isWaiting()) {
					// This job requires some monitoring. Start thread for this job.
					JobMonitor.create(job, mf);
				}
			}
		}
	}

	/**
	 * Method to add a sub-list to this job list.
	 * 
	 * This method must be used to add a given job sub-list to this
	 * list of jobs. This method fires notification events
	 * to all the interested listeners.
	 * 
	 * @param subList The sub-list of jobs to be added to this job
	 * list. This object cannot be null.
	 */
	public void addSubList(JobList subList) {
		// Let the insert method do the actual work.
		insert(subList, jobs.size());
	}
	
	/**
	 * Obtain the list of job and job-sublist objects that are 
	 * currently encapsulated in this job list.
	 * 
	 * @return The list of job and job-sublist  objects that are 
	 * currently encapsulated by this job list.
	 */
	public ArrayList<JobOrListWrapper> getJobs() { return jobs; }

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
		if (seqCounter > 0) {
			DOMHelper.addElement(srvrList, "SeqCounter", "" + seqCounter);
		}
		// Add sub-elements for each job in our job list
		for (JobOrListWrapper wrapper: jobs) {
			if (wrapper.isSubList()) {
				wrapper.getSubList().marshall(srvrList);
			} else {
				wrapper.getJob().marshall(srvrList);
			}
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
	 * 
	 * @param indentPrefix The extra indentation to be done to make the
	 * output look nice. If no additional indentation is needed then
	 * an empty string ("") must be passed in.
	 */
	public final void marshall(PrintWriter out, final String indentPrefix) {
		final String Indent = indentPrefix + "\t";
		final String STR_FORMAT = Indent + "\t<%1$s>%2$s</%1$s>\n";
		// Create a top-level job list entry for this class
		out.printf("%s<JobList>\n", Indent);
		// Add the sequence counter information to job list (if applicable)
		if (seqCounter > 0) {
			out.printf("%s\t<%2$s>%3$s</%2$s>\n", Indent, "SeqCounter", seqCounter);
		}
		// Serialize the name and description for the job list
		out.printf(STR_FORMAT, "Name", name);
		out.printf(STR_FORMAT, "Description", description);
		// Add sub-elements for each job in our job list
		for (JobOrListWrapper wrapper: jobs) {
			if (wrapper.isSubList()) {
				wrapper.getSubList().marshall(out, Indent);
			} else {
				wrapper.getJob().marshall(out, Indent);
			}
		}
		// Close the job list element.
		out.printf("%s</JobList>\n\n", Indent);
	}

	/**
	 * Obtain the description associated with this job list. 
	 * 
	 * The description is brief information that is either entered
	 * by the user or filled-in by PEACE when a job list is created.
	 * The top-level job list has a default description that is fixed.
	 * 
	 * @return the description The description associated with the
	 * job list. This method never returns null. If description is
	 * unavailable, then this method returns an empty string.
	 */
	public String getDescription() {
		return description;
	}

	/**
	 * Set the description for this job list.
	 * 
	 * @param description The description to be set for this job list.
	 */
	public void setDescription(String description) {
		this.description = description;
	}

	/**
	 * Obtain the name set for this job set. 
	 * 
	 * The name is typically a short string that provides the jest
	 * about this job list. This information is typically used to
	 * display sublists in hierarchies.
	 * 
	 * @return the name The name associated with this sublist.
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set the name for this sublist.
	 * 
	 * The name is typically a short string that provides the jest
	 * about this job list. This information is typically used to
	 * display sublists in hierarchies.
	 * 
	 * @param name The name be set for this sublist.
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * The list of job objects that have been configured and added
	 * to this list. This list is populated either via the create 
	 * method or new entries are added via the {@link #add(Job)} method
	 * in this class. 
	 */
	private ArrayList<JobOrListWrapper> jobs;

	/**
	 * The name associated with the job list. The name is used to
	 * display sublist entries in a hierarchical manner.
	 */
	private String name;

	/**
	 * A brief user-entered description for this job list. The
	 * description is primarily meant for the user to help
	 * with appropriately classifying jobs for future reference.
	 */
	private String description;

	/**
	 * Sequence counter that is maintained on a per-work space basis to
	 * generate unique/valid IDs for each new job entry added to this
	 * list.
	 */
	private long seqCounter;
}
