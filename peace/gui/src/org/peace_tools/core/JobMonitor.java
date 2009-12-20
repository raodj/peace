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
// Authors:   Dhananjai M. Rao              raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.core;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.SwingUtilities;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.UserLog;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * A background thread class to monitor job progress.
 *
 * This is a background thread class that is created to monitor
 * the progress of a Job. A job monitor is created either when a 
 * work space is loaded (and it has unfinished jobs pending) or
 * when a new job is added to the work space. This class is not
 * directly instantiable. Instead the create() method must be used
 * to create a monitor thread.
 */
public class JobMonitor implements Runnable {
	/**
	 * The primary API method to create a job monitor thread.
	 * 
	 * This method serves as the primary API for creating a Job
	 * monitor thread. This method creates an instance of this 
	 * class as a low priority daemon thread. The thread runs
	 * asynchronously periodically monitoring the job progress.
	 *  
	 * @param job Information regarding the job to be monitored.
	 * 
	 * @param listener The action listener to be notified when 
	 * the job has completed and this thread is about to be exited.
	 * The action event posted has the command set to "JobMonitor"
	 * and the object is set to the Job that has completed. 
	 * 
	 * @return This method returns true if the thread was started.
	 * 
	 */
	public static boolean create(Job job, ActionListener listener) {
		// Obtain the server reference for this job.
		Server server = Workspace.get().getServerList().getServer(job.getServerID());
		if (server == null) {
			// We don't have a valid server for this job. Can't
			// start up a job monitor.
			return false;
		}
		if (job.getMonitor() != null) {
			// A monitor thread already exists. Don't start another one.
			return false;
		}
		// Have a valid server and a job. Start up the job monitor.
		JobMonitor jm = new JobMonitor(job, server, listener);
		// Create a thread and start the monitor as a daemon.
		Thread jmThread = new Thread(MainFrame.getWorkerThreads(), jm, job.getJobID());
		jmThread.setDaemon(true);
		jmThread.setPriority(Thread.MIN_PRIORITY);
		// Set cross reference in job.
		job.setMonitor(jmThread);
		// Start the monitor thread.
		jmThread.start();
		// Everything went well
		UserLog.log(UserLog.LogLevel.NOTICE, "JobMonitor", 
				"Started job monitoring thread for job " + job.getJobID());
		return true;
	}
	
	/**
	 * Convenience method to interrupt a job monitor thread.
	 * 
	 * This method is a convenience method to interrupt the operations
	 * of a job monitoring thread.
	 *  
	 * @param job The job whose monitoring thread is to be interrupted.
	 * 
	 */
	public static void interrupt(Job job) {
		Thread monitor = job.getMonitor();
		if (monitor != null) {
			monitor.interrupt();
		}
	}
	
	/**
	 * Method to obtain an existing monitor thread for this job.
	 * 
	 * This method is a convenience method that can be used to obtain
	 * the job monitor associated with a given job.
	 * 
	 * @param job The job whose job monitor thread is to be returned.
	 * 
	 * @return The job thread (if any) associated with the job. If the
	 * job does not have a monitor thread associated with this then this
	 * method returns null.
	 */
	public static synchronized Thread getMonitor(Job job) {
		final String jobID   = job.getJobID();
		ThreadGroup group    = MainFrame.getWorkerThreads();
		Thread[] monitorList = new Thread[group.activeCount() + 10];
		int threadCount      = group.enumerate(monitorList);
		// Search in threads
		for(int i = 0; (i < threadCount); i++) {
			if (jobID.equals(monitorList[i].getName())) {
				return monitorList[i];	
			}
		}
		// No matching job thread found.
		return null;
	}
	
	@Override
	public void run() {
		// Keep looping away until the job status changes to a 
		// completion value.
		while ((!job.isDone()) && 
			   (!JobBase.JobStatusType.FINISHING.equals(job.getStatus()))) {
			// Create connection to the server to check job status
			if (!createSession()) {
				// The session creation was interrupted or canceled
				break;
			}
		
			// Check and update status of job.
			if (session != null) {
				if (!updateJobStatus()) {
					// No status change.
					try {
						Thread.sleep(recheckDelay);
					} catch (InterruptedException ie) {
						break;
					}
				}
			}
		}
		// When control drops here that means the job has finished
		// running or the thread ws interrupted. If job completed then
		// Either its status is FINISHING which means it
		// finished correctly and the output files are ready to be
		// copied to local machine or it failed. Either way report
		// the thread completion information through a suitable 
		// action event
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				ActionEvent ae = new ActionEvent(job, ActionEvent.ACTION_PERFORMED, "JobMonitor");
				listener.actionPerformed(ae);
			}
		});
		// Clear out reference in job
		job.setMonitor(null);
	}

	/**
	 * Helper method to run the jobRunner script on the remote machine
	 * and obtain the operational status of the job. This method is
	 * periodically invoked from the main run() method.
	 * 
	 * @return This method returns true if the status changed.
	 */
	private boolean updateJobStatus() {
		String extension = ".sh";
		try {
			if (ServerSession.OSType.WINDOWS.equals(session.getOSType())) {
				extension = ".bat";
			}
		} catch (Exception e) {
			// Exception occurred when detected OS type. This is some
			// concerning situation.
			UserLog.log(UserLog.LogLevel.WARNING,
					"JobMonitor", "Unable to determine OS type: " + e.getMessage());
			// Reset our session.
			session.disconnect();
			session = null;
			return false; // pretend status did not change
		}
		String remotePath = job.getPath() + "/jobRunner" + extension + " status";
		// Run job and get stdin and stdout
		String outputs[] = {"", ""};
		int exitCode = 0;
		try {
			exitCode = session.exec(remotePath, outputs);
		} catch (Exception e) {
			// Exception occurred when running the job. This is some
			// concerning situation.
			UserLog.log(UserLog.LogLevel.WARNING,
					"JobMonitor", "Server session error: " + e.getMessage());
			// Reset our session.
			session.disconnect();
			session = null;
			return false; // pretend status did not change
		}
		// Check if we have valid progress information.
		if ((exitCode != 0) || (outputs[1].length() > 0) || 
			(outputs[0].length() < 2)) {
			UserLog.log(UserLog.LogLevel.WARNING,
					"JobMonitor", "Job runner script error: " + outputs[1]);
			return false; // pretend status did not change
		}
		// Check if we have some progress information to report.
		// The progress information must have the following 2 lines:
		// estsAnalyzed, totalEstCount
		// running | done
		String progInfo[] = outputs[0].trim().split("\n");
		if (progInfo.length < 3) {
			// Invalid progress information.
			UserLog.log(UserLog.LogLevel.WARNING,
					"JobMonitor", "Progress information did not have" +
					" 3 lines (Raw info: "+ outputs[0] + ")");
			return false;
		}
		// Parse and process the first line.
		progInfo[0] = progInfo[0].trim();
		int comma = progInfo[0].indexOf(',');
		int estsAnalyzed = Integer.parseInt(progInfo[0].substring(0, comma));
		int estCount     = Integer.parseInt(progInfo[0].substring(comma + 1));
		// Update job progress information if progress has changed.
		int[] prevProgress = job.getProgressInfo();
		if (estsAnalyzed > prevProgress[0]) {
			// There has been advancement.
			job.setProgress(estsAnalyzed, estCount);
			// Check back again sooner so as not to miss progress
			recheckDelay = Math.max(250, recheckDelay / 2);
		} else {
			// No update. Don't bother checking so frequently then.
			recheckDelay = Math.min(60000, recheckDelay * 2);
		}
		
		// Next process and obtain overall runtime status.
		JobBase.JobStatusType status = (estCount > -1) ? JobBase.JobStatusType.RUNNING :
			JobBase.JobStatusType.QUEUED;
		progInfo[1] = progInfo[1].trim();
		if ("done".equals(progInfo[1])) {
			// The job is done running. Check and update status
			// further based on exit status.
			int exitStatus = Integer.parseInt(progInfo[2].trim());
			status = (exitStatus == 0) ? JobBase.JobStatusType.FINISHING :
				JobBase.JobStatusType.FAILED;
		}
		// Update the job status information if it has changed
		if (!status.equals(job.getStatus())) {
			job.setStatus(status);
			return true;
		}
		// no status change
		return false;
	}
	
	/**
	 * This is a helper method that is periodically invoked to 
	 * check and create a session to the remote server. This 
	 * method sets up the instance variable session once a valid
	 * connection to the local/remote server is established. 
	 * 
	 * @note If a valid session could not be established this method
	 * will simply sleep for 5 minutes before returning control.
	 * 
	 * @return This method returns false if the calling outer loop
	 * must stop because the user interrupted or canceled the session.
	 */
	private boolean createSession() {
		if (session != null) {
			return true;
		}
		
		try {
			// Synchronize on job monitor so that request for 
			// passwords don't all popup up at the same time when job
			// monitor threads are started up when a work space is loaded.
			synchronized (JobMonitor.class) {
				// First try to create a session.
				ServerSession tempSession = SessionFactory.createSession(null, server);
				// Set the purpose for this session.
				String purpose = String.format(PURPOSE_MSG, job.getJobID(), server.getName());
				tempSession.setPurpose(purpose);
				tempSession.connect();
				// Session is good. Update instance variable.
				session = tempSession;
			}
		} catch (IOException ioe) {
			UserLog.log(UserLog.LogLevel.WARNING, "JobMointor", 
					"Unable to connect to server: " + server.getName()
					+ "[Exception: " + ioe.getMessage() + "]");
			UserLog.log(UserLog.LogLevel.WARNING, "JobMointor", 
					"Unable to update status of job: " + job.getJobID());
		}
		
		if (session == null) {
			// Instead of trying again and again, this method simply
			// sleeps for a 5 minutes and try's again.
			try {
				Thread.sleep(5 * 60 * 1000L);
			} catch (InterruptedException e) {
				// We were interrupted. This is not a good sign
				return false;
			}
		}
		// Try again later
		return true;
	}
	
	/**
	 * The constructor merely initializes the instance variables to
	 * default values. The constructor has been made private to ensure
	 * this class is instantiated only via the create method.
	 * 
	 * @param job The job to be monitored by this job monitor.
	 * @param listener The action listener to be notified when this thread
	 * exits.
	 * @param server Information regarding the server on which this
	 * job is allegedly running.
 	 */
	private JobMonitor(Job job, Server server, ActionListener listener) {
		// Setup the job reference.
		this.job          = job;
		this.server       = server;
		this.listener     = listener;
		this.recheckDelay = 250; // 1/4 of a second
	}
	
	/**
	 * The job that this monitor is monitoring. This value is set
	 * in the constructor and is never changed during the life time
	 * of this object/thread.
	 */
	private final Job job;
	
	/**
	 * The action listener to be notified when this thread exists
	 * with a status of FINISHING or ERROR.
	 */
	private final ActionListener listener;
	
	/**
	 * Information regarding the server on which the job (that this 
	 * monitor is monitoring) is running. This value is set
	 * in the constructor and is never changed during the life time
	 * of this object/thread.
	 */
	private transient final Server server;
	
	/**
	 * The session (local/remote) to the server on which the job
	 * is running. The session is created when the thread starts 
	 * up.
	 */
	private transient ServerSession session;

	/**
	 * The delay in milliseconds for which the main thread must sleep
	 * before checking back for updates. This value is halved each time
	 * progress is discovered until this value falls to 250 milliseconds.
	 * Otherwise the value is doubled until this value reaches 120,000
	 * milliseconds.
	 */
	private transient int recheckDelay;
	
	/**
	 * A static formatable purpose message that is filled in and 
	 * used as the purpose for the server session initiated by this
	 * job monitor.
	 */
	private static final String PURPOSE_MSG =
		"<html>This session is being used to monitor a Job<br>" +
				"(<i>Job ID: %s</i>) running on the server named<br>" +
				"<i>%s</i>.<br>" +
				"The session will be closed once the job is complete.<br>" +
				"You may stop the monitor via the Job table.</html>";
}
