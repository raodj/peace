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

package org.peace_tools.decagon.djc;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.ProgressMonitorInputStream;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.FileInfo;
import org.peace_tools.core.JobMonitor;
import org.peace_tools.core.job.ServerPanel;
import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.decagon.DecagonVariables;
import org.peace_tools.decagon.helpers.DADXHelper;
import org.peace_tools.decagon.jaxb.AssemblerDetails;
import org.peace_tools.decagon.jaxb.Job;
import org.peace_tools.decagon.jaxb.OutputFile;
import org.peace_tools.decagon.jaxb.Process;
import org.peace_tools.generic.BackgroundTask;
import org.peace_tools.generic.BackgroundTask.UserTask;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProcessOutputDisplay.StyleKind;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DECAGONJob;
import org.peace_tools.workspace.DECAGONProcess;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.FileEntry.FileEntryType;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.JobSummary;
import org.peace_tools.workspace.Param;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * This wizard page generates the DECAGON jobs to be run on a
 * server and starts the first job.
 * 
 * This page is the last page that is run for the DECAGON Job
 * Creator (DJC) wizard. This page commences its operations
 * once it is displayed in a background Swing worker thread
 * to ensure the GUI remains operational. The background 
 * worker thread performs the following tasks:
 * 
 * Assume <i>J</i> jobs (where <i>J</i> is some fixed positive integer)
 * have to be created for the current DECAGON pipeline. This page performs 
 * the following steps (bear in mind that any of the following steps 
 * and sub-steps can generate an exception):
 * 
 * <ol>
 * 
 * <li>Establish session with the server. All jobs are going to be run on
 * the same server.</li>
 * 
 * <li>Check and copy the cDNA file and the source genes files to the
 * server.</li>
 * 
 * <li>Reserve <i>J</i> job IDs from the workspace (of course, if the jobs
 * are not successfully created then a part of the <i>J</i> IDs will be wasted.
 * But that is acceptable as the IDs are long) as jobs are created in reverse
 * order in the step below.</li>
 * 
 * <li>Create jobs starting with the <i>J</i><sup>th</sup> job and working 
 * down to the first job in the following manner (reverse order is used as 
 * the previous job needs to know the precise directory to copy files 
 * and the job to kick-off) 
 * 
 *     <ol>
 * 
 *     <li>Create directory for the current job on the server.</li>
 *     <li>Update all the DECAGON variables for the job.</li>
 *     <li>Create the job runner shell script using a template. The template
 *         has references to various DECAGON variables that may need to be 
 *         replaced.</li>
 *      <li>Create a job entry in PEACE workspace for the current job.</li>
 *      </ol>
 *      
 * </li>
 * 
 * <li>Once all the <i>J</i> jobs have been successfully created, start the
 * first job via the job runner script.</li>
 * 
 * </ol>
 */
public class JobSubmissionPage extends GenericWizardPage implements UserTask {
	/**
	 * The constructor. 
	 * 
	 * The constructor only sets up the minimal set of information
	 * associated with this wizard page (the constructor does not setup
	 * the GUI components). The core GUI components are created
	 * just before this wizard page is to be displayed (once 
	 * all the necessary information has been obtained from the
	 * user).
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public JobSubmissionPage(DecagonJobCreator dcm) {
		this.djc         = dcm;
		assert(this.djc != null);
		// Setup the title(s) for this page and border
		setTitle("Job(s) Submission", "Generate and submit jobs");
		setBorder(new EmptyBorder(5, 5, 5, 5));
	}
		
	/**
	 * Helper method to create GUI components on this wizard page.
	 * 
	 * This is a helper method that is used to create the 
	 * GUI components on this wizard page. The GUI components are
	 * created just before the wizard page is displayed. Specifically
	 * the {@link #taskRunner} object needs to be created each time
	 * the background operations are to be performed because of the way
	 * the Swing worker threads operate.
	 */
	private void createTaskRunner() {
		// Create the informational labels.
		JLabel infoLabel = new JLabel(INFO_MSG,
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		// Shortcuts to various informational components used to create
		// the wizard steps.
		final DADXHelper helper = djc.getDADXHelper();
		final List<Job> jobList = helper.getAssemblerDetails().getJob();
		final Server    srvr    = djc.getSelectedServer();
		final DataSet   ds      = djc.getDataSet();
		final DataSet   srcDS   = Workspace.get().findDataSet(ds.getSourceID());
		// Create information about each major operation that is going to
		// be performed to provide the user with some feedback.
		String[] subSteps = new String[jobList.size() + 4];
		subSteps[0]       = String.format(STEP_INFO_TMPL, "Connecting to Server", 
				Utilities.trim(srvr.getName(), 40));
		subSteps[1]       = String.format(STEP_INFO_TMPL, "Copy cDNA file (if needed)",
				Utilities.trim(ds.getPath(), 40));
		subSteps[2]       = String.format(STEP_INFO_TMPL, "Copy source Gene file (if needed)",
				Utilities.trim((srcDS != null ? srcDS.getPath() : "n/a"), 40));
		// Create step labels for jobs
		for(int i = 0; (i < jobList.size()); i++) {
			final Job jobEntry = jobList.get(i);
			subSteps[i + 3] = String.format(STEP_INFO_TMPL, jobEntry.getName(), 
					Utilities.trim(jobEntry.getDescription(), 40));			
		}
		// Create a step for starting the first job in the sequence of jobs
		subSteps[subSteps.length - 1] = String.format(STEP_INFO_TMPL, 
				"Start first job", "Other jobs are automatically run in sequence");
		// Clear out existing task runner object on this wizard page (if any)
		if (taskRunner != null) {
			this.remove(taskRunner.getTopPanel());
		}
		// Create our background helper and associated GUI components.
		taskRunner = new BackgroundTask(this, infoLabel, true, 
				subSteps, false, true, false, 1, true, false, false, null);
		// Make the log area a bit larger
		taskRunner.getLog().getTextPane().setPreferredSize(new Dimension(600, 150));
		// Ensure that the sub-labels are scrolled if there are several. 
		taskRunner.getSubStepScrollPane().setPreferredSize(new Dimension(100, 100));
		// Add the GUI panel to our our center panel.
		taskRunner.getTopPanel().setBorder(null);
		add(taskRunner.getTopPanel(), BorderLayout.CENTER);
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * 
	 * When the user is moving forward in the wizard (that is index of
	 * current page is higher than previous page), then this method 
	 * uses the {@link #taskRunner} background helper and starts up background
	 * operations. The background thread causes the {@link #run(BackgroundTask)}
	 * method in this class to be invoked. Each executable (if any) in the
	 * DADX file selected by the user is verified to ensure it is 
	 * operational. 
	 * 
	 * @param dialog The wizard dialog that logically owns this page.
	 * Currently this parameter is not used.
	 * 
	 * @param currPage The current zero-based index of this page. This
	 * parameter is used to avoid re-running checks when the user is
	 * backtracking.
	 * 
	 * @param prevPage The previous page from where the user navigated
	 * to this page. If the previous page is higher than the current page,
	 * then this method does not trigger the process of verifying
	 * executables. 
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		if (currPage < prevPage) {
			// The user is backtracking. Do not start background processing.
			return;
		}
		// Disable various button until background operations are done.
		djc.setButtonStatus(0, 0, -1);
		// Setup the various GUI components on this wizard page
		createTaskRunner();
		// Kick off the verification operation in the background after
		// the page is refreshed
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				taskRunner.start(false, JobSubmissionPage.this, "JobSubmission", null);
			}
		});
	}

	@Override
	public void run(BackgroundTask bTask) throws Exception {
		// Create local folder to store generated data files.
		File storeDir = new File(djc.getStorageDir());
		storeDir.mkdir();
		// Obtain a session (or connection) to the server.
		startSession();
		// One step done. Update progress.
		bTask.updateProgress();
		// Copy the cDNA file to be assembled to the server
		final DataSet ds = djc.getDataSet();
		copyESTFile(ds);
		bTask.updateProgress();
		// Copy the source cDNA file to be assembled to the server
		final DataSet srcDS = Workspace.get().findDataSet(ds.getSourceID());
		if (srcDS != null) {
			copyESTFile(srcDS);
		}
		bTask.updateProgress();
		// Obtain the list of DECAGON jobs to be created.
		AssemblerDetails asmInfo  = djc.getDADXHelper().getAssemblerDetails();
		final List<Job> decJobList = djc.getDADXHelper().getAssemblerDetails().getJob();
		// Shortcut to the PEACE-job-list where entries for DECAGON jobs are to be added
		final JobList peaceJobList = Workspace.get().getJobList();
		// Reserve JobIDs for jobs.
		final String reservedJobIDs[] = new String[decJobList.size()];
		for(int i = 0; (i < decJobList.size()); i++) {
			reservedJobIDs[i] = peaceJobList.reserveJobID();
		}
		// Next iterate over each job and create the job entries
		// in reverse order. Save the PEACE workspace-job entries created
		// in a job sub-list to be added to the workspace.
		final String subListName = asmInfo.getName() + " [" + 
			reservedJobIDs[0] + "-" + reservedJobIDs[decJobList.size() - 1] + "]";
		JobList jobSubList       = new JobList(subListName, djc.getDescription());
		final JobConfigPage[] jobConfigs = djc.getJobConfigs();
		DECAGONJob wsJobEntry    = null;
		for(int jobIdx = decJobList.size() - 1; (jobIdx >= 0); jobIdx--) {
			// Create job scripts (on server) and workspace job entries.
			wsJobEntry = createJob(decJobList.get(jobIdx), 
					reservedJobIDs[jobIdx],  // current job ID
					(jobIdx > 0 ? reservedJobIDs[jobIdx - 1] : null), // previous job
					jobConfigs[jobIdx].getServerInfoPanel(), 
					wsJobEntry);
			jobSubList.insert(wsJobEntry);
			// Create the list of generated files for this job.
			GeneratedFileList gfl = createGeneratedFileList(decJobList.get(jobIdx),
					wsJobEntry);
			// Add the list of generated files to the source data set.
			djc.getDataSet().insert(gfl);
			// Finished one major step in the process. Update progress.
			bTask.updateProgress();
		}
		// Add the job sublist to the workspace.
		Workspace.get().getJobList().addSubList(jobSubList);
		// Start the first job in the list of jobs.
		startJob(wsJobEntry, session);
		bTask.updateProgress();
	}
	
	/**
	 * This is a helper method that is used to create a valid
	 * directory where this job's files will be copied.
	 * 
	 * @param server The server object on which the job is to be created.
	 * The server object is used to determine the installation 
	 * location for PEACE under which job directories are created.
	 * 
	 * @param srvrSession The session on the server to be used for
	 * creating the directory. The session object must be a valid,
	 * current session on the server.
	 * 
	 * @param jobID The unique ID assigned to this job. This value
	 * is used to create the final job directory under which files
	 * for this job are to be stored.
	 * 
	 * @return This method returns the full path to the job directory
	 * on the server.
	 * 
	 * <p>NOTE: This method implicitly uses the {@link #taskRunner}
	 * object to cut logs.</p>
	 */
	private String createJobDirectory(final Server srvr, 
			final ServerSession srvrSession, 
			final String jobID) throws Exception {
		taskRunner.log("Attempting to create directory for job " + jobID + " on server...\n");
		// repeatedly try and create a work folder for the job based on its id
		int tryCount = 0;
		do {
			// Compute jobDir based on the try count
			String jobDir = srvr.getInstallPath() + "/jobs/" + jobID;
			if (tryCount > 0) {
				jobDir += "_" + tryCount;
			}
			// Check if the directory exists. If not create a new one
			// and use it.
			taskRunner.log("Checking if directory " + jobDir + " exists...");
			FileInfo fi = session.fstat(jobDir);
			if (!fi.exists()) {
				// OK, the directory does not exist
				taskRunner.log("No\n");
				taskRunner.log("Creating directory " + jobDir + " ...");
				srvrSession.mkdir(jobDir);
				taskRunner.log("Done.\n");
				taskRunner.log("The full job path is " + jobDir + "\n");
				// Break out of the while loop
				return jobDir;
			} else {
				taskRunner.log("Yes.\nTrying a different directory next.\n");
				tryCount++;
			}
		} while (tryCount < 5);
		// When control drops here we decide if we are going to
		// abort or not.
		if (tryCount >= 5) {
			// Creating work directory failed. Bail out.
			throw new IOException("Unable to create a useable work directory for job");
		}
		return "";
	}
	
	/**
	 * Helper method to run the job runner on the server.
	 * 
	 * This method is invoked from the run method to actually run the
	 * job runner on the remote machine. This method execute the job
	 * runner on the remote machine with "start" as the command line
	 * parameter.
	 */
	private void startJob(DECAGONJob job, ServerSession server) throws Exception {
		taskRunner.log("Starting job " + job.getJobID() + " on server...\n");
		job.setStatus(JobBase.JobStatusType.STARTING);
		Server.OSType os = server.getOSType();
		final String extension = (Server.OSType.WINDOWS.equals(os)) ? "bat" : "sh"; 
		String remotePath = job.getPath() + "/jobRunner." + extension
			+  " start";
		// Run job and get stdin and stdout
		String outputs[] = {"", ""};
		int exitCode = server.exec(remotePath, outputs);
		taskRunner.log("Output:\n");
		taskRunner.log(outputs[0]);
		taskRunner.log("Error Messages (if any):\n");
		taskRunner.log(outputs[1]);
		if ((exitCode != 0) || (outputs[1].length() > 0)) {
			// Error starting up job. Flag job has having failed.
			job.setStatus(JobBase.JobStatusType.FAILED);
			throw new IOException("The job runner generated an error.\n" +
			"The remote job has not started up successfully.");
		}
		taskRunner.log("Job started successfully.\n");
		// Start up the job monitor for the job.
		taskRunner.log("Starting job monitor for the job...\n");
    	JobMonitor.create(job, djc.getMainFrame());
		taskRunner.log("Done starting job monitor.\n");
		taskRunner.log("The GUI will perirodically check on the job status.\n");
	}
	
	/**
	 * Convenience method to create PEACE workspace entries for various output
	 * files associated with a job.
	 * 
	 * This method is a convenience method that is invoked from the 
	 * {@link #createJob(Job, String, String, ServerPanel)} method to build
	 * the list of generated files. This method uses the list of output
	 * files associated with the job to build the list of generated 
	 * files. The actual source for each file is set to be the server on
	 * which the files will be generated once the job runs successfully.
	 * The destination for each of the generated files is set to the
	 * storage directory (on the local machine) specified by the user.
	 * 
	 * @param dadxJob The DECAGON job entry whose list of output files are
	 * to be used to generate suitable workspace entries.
	 * 
	 * @param peaceJob The PEACE workspace job entry to be used for generating
	 * the file list. This object is used to create the necessary job summary
	 * information to be associated with the generated file list.
	 * 
	 * @return A PEACE workspace entry that contains information about
	 * the set of output files to be generated by the job.
	 * 
	 */
	private GeneratedFileList createGeneratedFileList(Job dadxJob, DECAGONJob peaceJob) {
		JobSummary summary    = new JobSummary(peaceJob);
		GeneratedFileList gfl = new GeneratedFileList(summary);
		DADXHelper dadxHelper = djc.getDADXHelper();
		for(OutputFile of: dadxJob.getOutputFile()) {
			final String fileID   = Workspace.get().reserveID();
			final String srvrPath = dadxHelper.processVariables(of.getPath());
			final File   srvrFile = new File(srvrPath);
			final String localPath= djc.getStorageDir() + File.separator + srvrFile.getName();
			FileEntry fe = new FileEntry(fileID,
					(of.isIsContigFile() ? FileEntryType.ASM : FileEntryType.OUTPUT),
					DADXHelper.toDataFileType(of.getFileType()),
					localPath, of.getDescription(), srvrPath);
					dadxHelper.processVariables(of.getPath());
			gfl.add(fe);
		}
		return gfl;
	}
	
	private DECAGONJob createJob(Job dadxJob, final String peaceJobID, 
			final String prevJobID, 
			final ServerPanel jobConfig, DECAGONJob nextJob) throws Exception {
		// Create the job directory on the server.
		final Server srvr   = jobConfig.getSelectedServer();
		final String jobDir = createJobDirectory(srvr, session, peaceJobID);
		// Update the various DECAGON variables for this job.
		djc.getDADXHelper().getVariables().setVariables(djc.getDataSet(), 
				srvr, jobDir, peaceJobID, prevJobID, 
				jobConfig.getPlatformConfiguration(), 
				(nextJob != null) ? nextJob.getPath() : null);
		// Setup the output contig file name.
		djc.getDADXHelper().setOutputContigVariable();
		// Log variables for troubleshooting potential problems when reported by users.
		taskRunner.log(djc.getDADXHelper().getVariables().toString());

		// Generate the command-line for each one of the processes and
		// setup corresponding DECAGON variables.
		ArrayList<DECAGONProcess> procEntryList = setProcessVariables(dadxJob, srvr);
		// Now load and generate the job runner script using the template
		// script that we have.
		String runnerScript = Utilities.readSmallTextFile("decagon/jobRunner.sh");
		// Substitute all DECAGON variables with suitable values.
		runnerScript = djc.getDADXHelper().processVariables(runnerScript);
		// Now copy the runnerScript to a remote file.
		String remotePath = jobDir + "/jobRunner.sh";
		taskRunner.log("Copying job runner script to: " + remotePath + "\n");
		// Create a string input stream to copy data.
		byte[] bytes = runnerScript.getBytes("UTF-8");
		ByteArrayInputStream is = new ByteArrayInputStream(bytes);
		session.copy(is, jobDir, "jobRunner.sh", 0700);
		taskRunner.log("Done copying job runner script.\n");
		// Create the workspace process entry.
		final int[] config = jobConfig.getPlatformConfiguration();
		DECAGONJob decJob = new DECAGONJob(peaceJobID, 
			(djc.getDescription() + "(" + dadxJob.getDescription() + ")"), 
			srvr.getID(), jobDir, config[0], config[1], config[2], config[3], 
			djc.getDADXHelper().getSourceDADXFilePaths().get(0),
			new ArrayList<Param>(), procEntryList);
		// Populate the actual parameters and process information.
		djc.getDADXHelper().getVariables().copyAll(decJob.getDECAGONVariables());
		// set previous job ID if supplied
		decJob.setPreviousJobID(prevJobID);
		// Set the job status to waiting by default
		decJob.setStatus(JobBase.JobStatusType.WAITING);
		// Return the workspace entry for further use.
		return decJob;
	}
	
	/**
	 * Convenience method to create process workspace-entries  and variables 
	 * for a given job.
	 * 
	 * This is an internal helper method that is invoked from the the
	 * {@link #createJob(Job, String, String, ServerPanel)} method to
	 * create the necessary process entries for each job. This method
	 * sets the {@link DecagonVariables#PROCESS_CMD_LINE_LIST} and
	 * {@link DecagonVariables#PROCESS_COUNT} variables after performing
	 * the following operations for each process associated with the job:
	 * 
	 * <ol>
	 * 	<li>It creates the final, fully resolved (with various variables
	 * substituted) process command-line via appropriate call to
	 * {@link DADXHelper#getProcessCmdLine(Process, Server)}.</li>
	 * 
	 *  <li>It concatenates the various process command-lines as a 
	 *  set of new-line separated entries.</li>
	 *  
	 * </ol>
	 * 
	 * @param dadxJob The DECAGON job for which process entries are to be
	 * created.
	 * 
	 * @param srvr The server on which the job is to run.
	 * 
	 * @return This method sets up the necessary DECAGON variables and
	 * returns a list of PEACE workspace process entries to be added to
	 * the workspace.
	 */
	private ArrayList<DECAGONProcess> setProcessVariables(Job dadxJob, Server srvr) {
		final List<Process> procList          = dadxJob.getProcess();
		ArrayList<DECAGONProcess> procEntries = new ArrayList<DECAGONProcess>();
		// Gather the command-line for each process as a string on a separate line.
		StringBuilder procCmdLineList = new StringBuilder(2048);
		for(Process proc: procList) {
			final String procCmdLine = 
				djc.getDADXHelper().getProcessCmdLine(proc, srvr);
			taskRunner.log("Command line for process: " + procCmdLine + "\n");
			// Gather the process command line into a buffer
			procCmdLineList.append("'");
			procCmdLineList.append(procCmdLine);
			procCmdLineList.append("'\n");
			// Add process entry to workspace entry.
			DECAGONProcess wsEntry = new DECAGONProcess(procCmdLine);
			procEntries.add(wsEntry);
		}
		// Remove the trailing newline
		procCmdLineList.deleteCharAt(procCmdLineList.length() - 1);
		// Setup the DECAGON variable with various process command lines
		final DecagonVariables decVars = djc.getDADXHelper().getVariables();
		decVars.addVariable(DecagonVariables.PROCESS_CMD_LINE_LIST, procCmdLineList.toString());
		decVars.addVariable(DecagonVariables.PROCESS_COUNT, Integer.toString(procList.size()));
		// Set output file list.
		StringBuilder outFileList = new StringBuilder(1024);
		for(OutputFile of: dadxJob.getOutputFile()) {
			final String outFilePath = djc.getDADXHelper().processVariables(of.getPath());
			// Gather the process command line into a buffer
			outFileList.append("'");
			outFileList.append(outFilePath);
			outFileList.append("'\n");
		}
		// Remove the trailing newline
		outFileList.deleteCharAt(outFileList.length() - 1);
		// Setup the DECAGON variable with various process command lines
		decVars.addVariable(DecagonVariables.OUTPUT_FILE_LIST, outFileList.toString());
		decVars.addVariable(DecagonVariables.OUTPUT_FILE_COUNT, 
				Integer.toString(dadxJob.getOutputFile().size()));		
		// Return the workspace process entries
		return procEntries;
	}

	@Override
	public void done(BackgroundTask bTask) {
		// Compute overall status.
		final boolean success = (bTask.getResult() == null);
		if (!success) {
			// Mark jobs created (if any) as having failed.
			for(int i = 0; (i < jobEntryList.size()); i++) {
				final DECAGONJob je = jobEntryList.get(i);
				je.setStatus(je.getPreviousJobID() == null ? JobBase.JobStatusType.FAILED :
					JobBase.JobStatusType.WAIT_FAILED);
			}
			// Display an error message to draw attention to the user
			// even though the log already has the same information.
			JPanel msg = Utilities.collapsedMessage(GENERIC_ERROR_MSG, 
					Utilities.toString(bTask.getResult()));
			JOptionPane.showMessageDialog(this, msg, 
					"DECAGON Error", JOptionPane.ERROR_MESSAGE);
		}
		// Enable various button once background operations are done.
		djc.setButtonStatus(1, (success ? 1 : 0), -1);
	}

	@Override
	public void dialogClosed(BackgroundTask bTask) {
		// Currently this method (that implements an interface from
		// BackgroundTask.UserTask interface) does not have any
		// operations to perform and is intentionally left blank.
	}
	
	/**
	 * This is a helper method that is used to establish a valid
	 * session with the remote server on which the job is to be
	 * submitted. This method sets up the {@link #session}
	 * instance variable. This method must be invoked only after
	 * a valid background task has been created and the
	 * {@link #taskRunner} object has been initialized.
	 */
	private void startSession() throws IOException {
		Server srvr = djc.getSelectedServer();
		taskRunner.log("Attempting to connect to server " + srvr.getName() + "\n");
		session = SessionFactory.createSession(djc, srvr);
		// Setup purpose so the user knows why we are connecting
		session.setPurpose("Attempting to start a new job");
		// Connect to the server.
		session.connect();
		// Done. update status.
		taskRunner.log("Successfully connected.\n");
	}
	
	/**
	 * Helper method to copy a cDNA data file to shared space if needed.
	 * This method checks to see if the cDNA data file already exists
	 * on the server. If not it copies the cDNA data file over.
	 * 
	 * @param ds The data set whose cDNA file is to be copied to the
	 * server. The server selected by the user (obtained via
	 * {@link DecagonJobCreator#getSelectedServer()}) is used along
	 * with the current session in {@link #session}.
	 */
	private void copyESTFile(final DataSet ds) throws IOException {
		// Do some basic validation on local EST file.
		String localFileName = ds.getPath();
		taskRunner.log("Checking local EST FASTA file: " + localFileName + "\n");
		File localFile = new File(localFileName);
		if (!localFile.exists() || !localFile.isFile() || 
				!localFile.canRead()) {
			throw new IOException("Unable to read local cDNA file: " + localFileName);
		}
		taskRunner.log("Local EST file exists and is readable.\n");
		// Obtain information about the corresponding remote file.
		// Construct server-specific EST file location.
		final Server srvr = djc.getSelectedServer();
		final String serverESTFileName = srvr.getServerPath(localFile); 
		
		taskRunner.log("Checking information on " + serverESTFileName + "\n");
		FileInfo serverFile = session.fstat(serverESTFileName);
		if (serverFile.exists()) {
			// The file seems to exist on the server. Do some 
			// additional validation.
			taskRunner.log("The file exists on the server. Validating...\n");
			if (!serverFile.isFile() || !serverFile.canRead()) {
				throw new IOException("The EST file seems to exist on the " +
				"server but with is not a regular file that is readable.");
			}
			IOException ioe = null;
			if (serverFile.getSize() != localFile.length()) {
				ioe = new IOException("The EST file exist on the " +
				"server but does not match the size of the local EST file.");
			}
			if (serverFile.getLastModified() < localFile.lastModified()) {
				ioe = new IOException("The EST file with same size exists " +
				"on the server but has an earlier modification date.");
			}
			if (ioe == null) {
				taskRunner.log("The remote cDNA file is up to date. Skipping file copy\n");
				return;
			}
			// There is an inconsistency in the file. Let user decide if
			// the file should be over written.
			taskRunner.log(ioe.getLocalizedMessage(), StyleKind.WARNING);
			JPanel msg = Utilities.collapsedMessage(FILE_MISMATCH, 
					Utilities.toString(ioe));
			int decision = JOptionPane.showConfirmDialog(djc, msg, 
						"cDNA File Mismatch", JOptionPane.YES_NO_OPTION);
			if (decision == JOptionPane.NO_OPTION) {
				taskRunner.log("User decided not to proceed.", StyleKind.WARNING);
				throw ioe;
			}
		} else {
			taskRunner.log("cDNA file is not present on server. Must copy it.\n", StyleKind.INFO);
		}
		// When control drops here that means we need to copy the
		// local file to the server.
		String msg = "Copying cDNA file to Server.\n";
		taskRunner.log(msg);
		ProgressMonitorInputStream pmis =
			new ProgressMonitorInputStream(this, msg, 
					new FileInputStream(localFile));
		pmis.getProgressMonitor().setNote("File: " + localFile.getAbsolutePath());
		// Let the server do the actual copying.
		taskRunner.log("Copying file: " + localFile.getAbsolutePath() +
				" to " + serverESTFileName + "\n", StyleKind.BOLD_INFO);
		session.copy(pmis, "", serverESTFileName, 0600);
		// Done copying.
		taskRunner.log("Done copying.\n");
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to obtain summary information
	 * from various pages constituting this wizard.
	 */
	private final DecagonJobCreator djc;

	/**
	 * This is a generic object that is used to run checks for
	 * executables in the background. This class provides features
	 * for displaying status information about the ongoing checks
	 * in the background thread.
	 */
	private BackgroundTask taskRunner;

	/**
	 * The server session that is used to interact with a server to
	 * copy files, create jobs, and start the first job. This object
	 * is created by the 
	 */
	private ServerSession session;
	
	/**
	 * The list of workspace job entries created by this class.
	 * Entries are added to this list by the {@link #run(BackgroundTask)}
	 * method. These entries are used to flag jobs as having failed
	 * when exceptions occur in the {@link #done(BackgroundTask)}
	 * method.
	 */
	private ArrayList<DECAGONJob> jobEntryList = 
		new ArrayList<DECAGONJob>();
	
	/**
	 * A simple static informational message that is displayed at the 
	 * top of this page.
	 */
	private static final String INFO_MSG = 
		"<html>Please wait while various jobs associated with<br/>" +
		"the genome assembler pipeline are created and submitted.</html>";

	/**
	 * A simple template string that is formatted to generate 
	 * a 2-line status messages for each step in the wizard.
	 */
	private static final String STEP_INFO_TMPL =
		"<html>%s<br/><font size=\"-2\">%s</font></html>";

	/**
	 * A generic informational message that is displayed to the user
	 * when EST file mismatch is detected
	 */
	private static final String FILE_MISMATCH = 
		"<html>Error the remote cDNA file is not consistent with the local<br>" +
		"copy of the file. Do you want to overwrite the copy on the server?</html>";

	/**
	 * A generic informational message that is displayed to the user
	 * in the {@link #done(BackgroundTask)} method if an exception
	 * was generated by various operations performed by the background
	 * thread.
	 */
	private static final String GENERIC_ERROR_MSG = 
		"<html>An unrecoverable error occurred during the process<br/>" +
		"of creating, submitting, and starting the jobs associated<br/>" +
		"with the DECAGON pipeline. Please refer to the log messages<br/>" +
		"for additional details.</html>";

	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
