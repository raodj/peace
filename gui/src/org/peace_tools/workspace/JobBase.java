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

import java.util.GregorianCalendar;
import javax.xml.datatype.DatatypeConfigurationException;
import javax.xml.datatype.DatatypeFactory;
import javax.xml.datatype.Duration;

/**
 * This is a common base class that contains fields that are common to both
 * the JobSummary class and the Job class. This class does not have a 
 * corresponding element in the XML representation.
 */
public abstract class JobBase {
	/**
	 * Different enumerations defining the the type of job
	 * identified by a given job entry.
	 */
	public enum JobType {
		/**
		 * This entry indicates that the job entry corresponds to a
		 * clustering type job. This type of job has a lot of 
		 * configuration information including: filters, heuristics,
		 * and other parameters. This job is created for generating
		 * MST entries and clustering files.
		 */
		CLUSTERING,
		/**
		 * This type of job entry corresponds to a assembly job.
		 * This entry is created when assembly is performed via PEACE.
		 * This type of job entry does not include filters or
		 * heuristics. However, it does have parameters used to
		 * customize the operation of EAST.
		 */
		EAST,
		/**
		 * This type of job entry corresponds to a baton job that is
		 * used for clustering and assembly. 
		 */
		BATON,
		/**
		 * This type of job entry corresponds to a DECAGON job that
		 * is used for empirical evaluation of genomic-assemblers.
		 */
		DECAGON
	};
	
	/**
	 * Different enumerations defining the last known runtime status of
	 * a given Job. These enumerations were introduced to reflect those
	 * used in the XML and to ensure that the code is overall more
	 * readable.
	 */
	public enum JobStatusType {
		/**
		 * This status indicates that the GUI is making attempt to start the
		 * job. This involves creating a folder and copying any additional
		 * configuration information to that folder.
		 */
		STARTING,
		/**
		 * This status indicates that the job is waiting on a previous
		 * job to finish. This can happen when a assembly  job is waiting
		 * on a clustering job to finish.
		 */
		WAITING,
		/**
		 * The job has been queued for running, but it has not yet started
		 * running. This can happen if a cluster is overloaded with jobs.
		 */
		QUEUED,
		/**
		 * The job has got CPU time and is currently running.
		 */
		RUNNING,
		/**
		 * The job has finished running and the GUI is downloading data files
		 * and status information from the server.
		 */
		FINISHING,
		/**
		 * The job has completed running successfully and all the necessary
		 * data files are already available in this workspace.
		 */
		SUCCESS,
		/**
		 * The job has completed running but error(s) c
		 */
		FAILED,
		/**
		 * The job that this job was waiting on failed and therefore this
		 * job has been flagged as having failed as well.
		 */
		WAIT_FAILED
	};
	
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
	 */
	public JobBase(JobType type, String jobID, String serverID) {
		this.type      = type;
		this.jobID     = jobID;
		this.serverID  = serverID;
		// Set non-final fields to null for now.
 		this.status    = JobStatusType.STARTING;
		this.runtime   = null;
		this.prevJobID = null;
	}

	/**
	 * A copy constructor to create/initialize a common Job base object using
	 * data from an existing job or job summary object.
	 * 
	 * @param jobData The source job base object from where the data for this
	 * object is to be copied.
	 */
	public JobBase(JobBase jobData) {
		this.type     = jobData.type;
		this.jobID    = jobData.jobID;
		this.serverID = jobData.serverID;
		this.status   = jobData.status;
		this.runtime  = jobData.runtime;
		this.prevJobID= jobData.prevJobID;
	}
	
	/**
	 * Obtain the work space wide unique identifier set for this job. The
	 * jobID value is created when new jobs are scheduled. These IDs are
	 * persisted in the work space configuration file and loaded when a 
	 * work space is opened in the GUI.
	 * 
	 * @return This method returns the unique identifier set for this job.
	 */
	public String getJobID() { return jobID; }
	
	/**
	 * Obtain the work space wide unique identifier set for the server
	 * on which this job is scheduled to run (or was run). The
	 * serverID value is set when new jobs are scheduled. This ID is
	 * persisted in the work space configuration file and loaded when a 
	 * work space is opened in the GUI.
	 * 
	 * @return This method returns the unique server identifier set 
	 * for this job.
	 */
	public String getServerID() { return serverID; }

	/**
	 * Change the status for this job.
	 * 
	 * @param status The new status value to be set for this job.
	 */
	public void setStatus(JobStatusType status) {
		this.status = status;
	}
	
	/**
	 * Obtain the current status set for this job. The job status value
	 * is the last known status value for this job.
	 * 
	 * @return The current status for this job.
	 */
	public JobStatusType getStatus() { return status; }

	/**
	 * This method can be used to explicitly set the actual wall clock time that
	 * this job used on the server.
	 * 
	 * @param milliseconds The actual wall clock time that this job used on the
	 * server.
	 */
	public void setRunTime(long milliseconds) {
		try {
			DatatypeFactory codec = DatatypeFactory.newInstance();
			runtime = codec.newDuration(milliseconds);
		} catch (DatatypeConfigurationException e) {
			// Cut log entry in programmer log.
			e.printStackTrace();
		}
	}

	/**
	 * Obtain the wall clock run time for this job. This value is meaningful
	 * only after the job has completed running. Until the job has completed
	 * running (and its status is updated) the runtime is returned to be -1.
	 * 
	 * @return The wall clock time (in milliseconds) taken by this job to run
	 * on the server.
	 */
	public long getRunTime() {
		long time = -1;
		if (runtime != null) {
			time = runtime.getTimeInMillis(GregorianCalendar.getInstance());
		}
		return time;
	}

	/**
	 * Convenience method to determine if job has completed.
	 * 
	 * This method returns true only if the status of the job is
	 * either SUCCESS, FAILED, or WAIT_FAILED (a job this job 
	 * was depending on failed). In all other states this method
	 * returns false.
	 * 
	 * @return Returns true to indicate if the job has completed
	 * (successfully or otherwise).
	 */
	public boolean isDone() {
		return (JobStatusType.SUCCESS.equals(status) ||
				JobStatusType.FAILED.equals(status)  ||
				JobStatusType.WAIT_FAILED.equals(status));
	}

	/**
	 * Convenience method to determine if job is waiting on another
	 * job to complete..
	 * 
	 * This method returns true only if the status of the job is
	 * WAITING. In all other states this method returns false.
	 * 
	 * @return Returns true to indicate if the job is waiting on
	 * another job to complete.
	 */
	public boolean isWaiting() {
		return (JobStatusType.WAITING.equals(status));
	}
	
	/** Obtain the job type associated with this job.
	 * 
	 * @return The job type set for this job when it was created.
	 */
	public JobType getType() {
		return type;
	}
	
	/**
	 * Set the ID of the previous job that this job is dependent on.
	 * 
	 * This method must be used set the ID of the job that this job
	 * is dependent on.
	 * 
	 * @param prevJobID The ID of the job that this job is dependent
	 * on.
	 */
	public void setPreviousJobID(String prevJobID) {
		this.prevJobID = prevJobID;
	}
	
	/**
	 * The ID (if any) of a previous job that this job is dependent on.
	 * 
	 * This method can be used to determine the ID of an earlier/previous
	 * job that this job is/was dependent on.
	 * 
	 * @return The ID of a previous/dependent job if applicable. 
	 * Otherwise this method returns null.
	 */
	public String getPreviousJobID() {
		return prevJobID;
	}
	
	/**
	 * The unique generated jobID value for this job. This value is typically
	 * created when a new job is scheduled. This value is persisted in the
	 * work space configuration. 
	 */
	protected final String jobID;
	
	/**
	 * A cross reference to the server on which this Job has been scheduled.
	 * This value is used to look up job status information and obtain the actual
	 * resulting data files for this job.
	 */
	protected final String serverID;

	/**
	 * The job type value associated with this job entry. This value
	 * is set when this class is created and is never changed after
	 * that.
	 */
	protected final JobType type;

	/**
	 * The unique generated jobID value for another/previous job that 
	 * this job is dependent on. This value is typically set when a 
	 * clustering + assembly type jobs are scheduled. The assembly job
	 * is flagged as being dependent on the clustering job.
	 */
	protected String prevJobID;
	
	/**
	 * The current runtime status of this job. This value is periodically
	 * updated until the job completes running (successfully or with erros).
	 */
	protected JobStatusType status;
	
	/**
	 * The actual CPU time that the job took to run on the server. This
	 * value is set only after the job completes. Until the job completes
	 * this value is null.
	 */
	protected Duration runtime;
}
