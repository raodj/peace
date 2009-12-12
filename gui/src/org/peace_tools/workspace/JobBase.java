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
		FAILED
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
	public JobBase(String jobID, String serverID) {
		this.jobID               = jobID;
		this.serverID            = serverID;
		// Set non-final fields to null for now.
 		this.status              = JobStatusType.STARTING;
		this.runtime             = null;
	}

	/**
	 * A copy constructor to create/initialize a common Job base object using
	 * data from an existing job or job summary object.
	 * 
	 * @param jobData The source job base object from where the data for this
	 * object is to be copied.
	 */
	public JobBase(JobBase jobData) {
		this.jobID    = jobData.jobID;
		this.serverID = jobData.serverID;
		this.status   = jobData.status;
		this.runtime  = jobData.runtime;
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
	 * either SUCCESS or FAILED. In all other states this method
	 * returns true.
	 * 
	 * @return Returns true to indicate if the job has completed
	 * (successfully or otherwise).
	 */
	public boolean isDone() {
		return (JobStatusType.SUCCESS.equals(status) ||
				JobStatusType.FAILED.equals(status));
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
