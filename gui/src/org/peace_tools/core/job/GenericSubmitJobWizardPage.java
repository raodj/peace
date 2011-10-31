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

package org.peace_tools.core.job;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ProgressMonitorInputStream;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.FileInfo;
import org.peace_tools.core.JobMonitor;
import org.peace_tools.core.MainFrame;
import org.peace_tools.core.Version;
import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.JobBase.JobType;
import org.peace_tools.workspace.Server.EXEKind;
import org.peace_tools.workspace.Server.OSType;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as a common base class to create the final page 
 * in a job creation wizard. This page provides a lot of common
 * functionality that is needed for creating jobs. However, some
 * custom operations that are specific to a given type of job are
 * delegated to the derived class. When this page is displayed in
 * a wizard, it temporarily disables further navigation and performs
 * the following operations from a separate thread:
 * 
 *  <ol>
 *  
 *   <li>First a directory is created on the remote server for this
 *   job. The directory is created under the peace installation path
 *   on the remote server. This task is performed by the 
 *   {@link #createJobDirectory()} method and it uses the 
 *   {@link #getJob()} method </li>
 *   
 *   <li>New job entries and generated file entries are added to
 *   the workspace by the {@link #createEntries()} method. This
 *   method uses the {@link #getJob()} and {@link #getGeneratedFiles()}
 *   methods for its operation.</li>
 * 
 * 	 <li>Next, the existence of the source EST data file on the remote
 *   server is checked. If the file is not present then it is copied
 *   to the remote server. Otherwise the existing file is verified 
 *   and if it matches it is used. This task is performed by the
 *   {@link #copyESTFile()} method and uses the {@link #getDataSet()} call.</li>
 * 
 *   <li>After that, any additional data files used by the job are 
 *   copied from the local workspace to the remote server by the
 *   {@link #copyOtherFiles()} methods. This method uses the
 *   {@link #getFilesToCopy()} method.</li>
 *   
 * 	 <li>The job runner scripts are variable-substituted and
 *   copied to the install folder. The job runner script is executed
 *   to submit the job to run on the server. These operations are
 *   performed by the {@link #deployRunnerScript()} method. </li>
 *   
 *   <li>Finally, if all of the aforementioned steps complete
 *   successfully, the {@link #createJobUpdateThread()} method is
 *   called. It creates a job monitoring thread to periodically
 *   poll and update the status of the job.</li>
 *   
 *  </ol>
 * 
 * <p>Note that this class is not meant to be directly used. Instead
 * one of the derived classes must be instantiated and used.</p>
 * 
 * 
 * <p>Note that the process of submitting the job can take several
 * seconds (or even a minute if EST file to copied is large). 
 * Consequently the operations are run on a separate thread to ensure
 * that the GUI continues to provide feedback without appearing as if
 * it is hanging (it would appear to be hanging if the operations are
 * performed on the same thread). </p>
 */
public abstract class GenericSubmitJobWizardPage extends GenericWizardPage 
implements Runnable {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: A series of labels
	 * to show the various steps being performed and a JTextArea to
	 * display additional log information.
	 * 
	 * @param wizard The wizard dialog that logically owns this page.
	 * This object is merely used to create GUI windows (if needed).
	 * 
	 * @param title The title to be set for this wizard page.
	 * 
	 * @param subTitle The subtitle to be set for this wizard page.
	 * 
	 * @param startJob Flag to indicate if the job must be started after
	 * it is created. If the flag is true, then the job is started.
	 * Otherwise the job is not started (as it may be dependent on another
	 * job to finish). 
	 */
	public GenericSubmitJobWizardPage(WizardDialog wizard, String title,
			String subTitle, boolean startJob) {
		this.dialog   = wizard;
		this.startJob = startJob;
		assert(this.dialog != null);
		// Setup the title(s) for this page and border
		setTitle(title, subTitle);
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the job summary information text area.
		log = new JTextArea(3, 10);
		JScrollPane jsp = new JScrollPane(log);
		jsp.setMinimumSize(log.getPreferredSize());
		JComponent logBox = 
			Utilities.createLabeledComponents("Additional information:",
					null, 0, false, jsp);
		// Create the informational labels.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/16x16/Information.png"), 
				JLabel.LEFT);
		// Create the various steps this wizard is going to perform
		String StepTitle[] = {"Adding new workspace entries",
				"Create directory for job", 
				"Copy EST file (if needed)", 
				"Copy other data files (if needed)",
				"Deploy job runner script",
		"Submit/start job on server"};
		JPanel stepContainer = new JPanel(new GridLayout(StepTitle.length, 1));
		// Indent the steps a bit by 12 pixels
		stepContainer.setBorder(new EmptyBorder(0, 12, 0, 0));
		stepInfo = new JLabel[StepTitle.length];
		// If not starting the job then we skip last step.
		final int StepCount = StepTitle.length - (startJob ? 0 : 1);
		for(int i = 0; (i < StepCount); i++) {
			// Crate the label
			stepInfo[i] = new JLabel(StepTitle[i],
					Utilities.getIcon("images/16x16/Box.png"),
					JLabel.LEFT);
			stepContainer.add(stepInfo[i]);
		}
		// Back the title and step info labels into a container
		JPanel container = new JPanel(new BorderLayout(5, 5));
		container.add(info, BorderLayout.NORTH);
		container.add(stepContainer, BorderLayout.SOUTH);

		// Pack the display fields into a suitable panel
		JPanel subPanel = new JPanel(new BorderLayout(5, 5));
		subPanel.add(container, BorderLayout.NORTH);
		subPanel.add(logBox, BorderLayout.CENTER);
		// Now add the sub-panel to the main panel.
		add(subPanel, BorderLayout.CENTER);
	}

	/**
	 * This method is called just before this page is to be displayed.
	 * This page disables all the buttons in the wizard and jump
	 * starts the background thread that performs the various tasks. 
	 * 
	 * @param dialog The wizard dialog that owns this wizard page.
	 * @param currPage The current page to be displayed.
	 * @param prevPage The previous page that was displayed.
	 */
	@Override
	public final void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// First disable all the buttons in the bottom.
		dialog.setButtonStatus(0, 0, 0);
		// Now jump start the background thread.
		final Thread submitThread = new Thread(this);
		submitThread.setDaemon(true);
		// Start up the thread later on after this page is displayed
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				submitThread.start();
			}
		});
	}

	/**
	 * Helper method to update the status of a step.
	 * 
	 * This method is used to provide visual feedback to the user
	 * as to which step is being performed by this wizard.
	 * 
	 * @param stepIndex The index of the step whose status is to be
	 * updated.
	 * 
	 * @param performing If this flag is true, then this step is 
	 * underway. If it is false, it is assumed that the step has 
	 * been completed.
	 */
	private void setStepStatus(final int stepIndex, final boolean performing) {
		// The labels to update and repaint suitably
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				// Set current step font to be bold
				Utilities.adjustFont(stepInfo[stepIndex], 0, 6, performing ? 1 : -1);
				// Set the icon to be displayed 
				if (performing) {
					stepInfo[stepIndex].setIcon(Utilities.getIcon("images/16x16/RArrow.png"));
				} else {
					stepInfo[stepIndex].setIcon(Utilities.getIcon("images/16x16/CheckedBox.png"));
				}		
				dialog.getPagePane().repaint();
			}
		});
	}

	/**
	 * Method to perform various tasks involved in successfully
	 * creating and submitting a job on a remote server. This
	 * method is called just after the wizard page has been
	 * displayed from a separate thread.
	 */
	@Override
	public void run() {
		dialog.addThread(Thread.currentThread()); // Track threads
		try {
			// First add necessary work space entries.
			createEntries();
			// Next establish connection to the server
			startSession();
			// Next create a working directory for the job and update
			// the job's path.
			createJobDirectory();
			// Now check and copy the EST file to the server. The EST
			// files are all copied to a shared area to minimize 
			// unnecessary EST copying.
			copyESTFile();
			// Copy other files from the workspace to the remote server
			copyOtherFiles();
			// Now create and copy the job submission script for the
			// server.
			deployRunnerScript();
			// Finally run the job runner script if directed
			Job job = getJob();
			
			if (!startJob) {
				// The job was not started. Set its status to waiting.
				job.setStatus(JobBase.JobStatusType.WAITING);
			} else {
				// Run the job runner script to launch the job
				startJob();			
				// OK, the job has been queued. Change status
				job.setStatus(JobBase.JobStatusType.QUEUED);
				// Create a background thread to poll and update job
				// status information.
				createJobUpdateThread();
			
				// OK, all done. Let the user know the whole thing was
				// a resounding success.
				JOptionPane.showMessageDialog(dialog, SUCCESS_MESSAGE,
						"Job successfully started", JOptionPane.INFORMATION_MESSAGE);
			}
		} catch (Exception exp) {
			// Change status of the job to failed.
			Job job = getJob();
			job.setStatus(JobBase.JobStatusType.FAILED);
			JPanel msg = Utilities.collapsedMessage(JOB_SUBMIT_ERROR, 
					Utilities.toString(exp));
			JOptionPane.showMessageDialog(dialog, msg, 
					"Job Submission Error", JOptionPane.ERROR_MESSAGE);
		} finally {
			// Wind up server session if valid.
			if (server != null) {
				server.disconnect();
			}
			// Let the wizard know this thread is done
			dialog.removeThread(Thread.currentThread());
			// Enable the next/finish button in the wizard.
			dialog.setButtonStatus(0, 1, -1);
			// Setup saving the work space.
			MainFrame mf = (MainFrame) dialog.getParent();
			mf.saveDelayedWorkspace();
			// Finally do call-back method that the thread is done.
			done();
		}
	}

	/**
	 * Helper method to copy a given set of files to the job directory
	 * on the server.
	 * 
	 * This method is invoked from the run() method to adapt the local
	 * runner script (which is a template) by filling in the template
	 * parameters with actual values and copying the script to the
	 * remote machine.
	 * 
	 * @throws Exception This method throws exceptions on errors.
	 * 
	 * @see #getFilesToCopy()
	 */
	private void copyOtherFiles() throws IOException {
		setStepStatus(3, true);
		// Get the file entries to be copied.
		FileEntry[] filesToCopy = getFilesToCopy();
		if ((filesToCopy == null) || (filesToCopy.length < 1)) {
			// No files to copy. All done with this step
			setStepStatus(3, false);
			return;
		}
		// The location on the server for job-specific files.
		final String jobPath = getJob().getPath();
		// Copy each one of the files.
		for(FileEntry fe: filesToCopy) {
			// Do some basic validation on local EST file.
			String localFileName = fe.getPath();
			log.append("Checking local file: " + localFileName + "\n");
			File localFile = new File(localFileName);
			if (!localFile.exists() || !localFile.isFile() || 
				!localFile.canRead()) {
				throw new IOException("Unable to read local file: " + localFileName);
			}
			log.append("Local file exists and is readable.\n");
			String msg = "Copying file to Server.\n";
			log.append(msg);
			ProgressMonitorInputStream pmis =
				new ProgressMonitorInputStream(this, msg, 
						new FileInputStream(localFile));
			pmis.getProgressMonitor().setNote("File: " + localFile.getAbsolutePath());
			// Let the server do the actual copying.
			String serverFileName = jobPath + File.separator + localFile.getName();
			log.append("Copying file: " + localFile.getAbsolutePath() +
					" to " + serverFileName + "\n");
			server.copy(pmis, "", serverFileName, 0600);
			// Done copying.
			log.append("Done copying file" + localFileName + "\n");
		}
		setStepStatus(3, false);
	}

	
    /**
     * Helper method to create a thread to periodically poll a job.
     * 
     * This is a helper method that is used to create a 
     * daemon thread to check status of the job. The thread is created 
     * as a daemon thread so that the system can exit out of the 
     * thread in a graceful manner.
     */
    protected void createJobUpdateThread() {
    	final Job job = getJob();
    	MainFrame mf  = (MainFrame) dialog.getParent();
    	JobMonitor.create(job, mf);
    }
    
	/**
	 * Helper method to run the job runner on the server.
	 * 
	 * This method is invoked from the run method to actually run the
	 * job runner on the remote machine. This method execute the job
	 * runner on the remote machine with "start" as the command line
	 * parameter.
	 */
	private void startJob() throws Exception {
		setStepStatus(5, true);
		log.append("Submitting/starting job on server...\n");
		Server.OSType os = server.getOSType();
		final String extension = (Server.OSType.WINDOWS.equals(os)) ? "bat" : "sh"; 
		Job job = getJob();
		String remotePath = job.getPath() + "/jobRunner." + extension
			+  " start";
		// Run job and get stdin and stdout
		String outputs[] = {"", ""};
		int exitCode = server.exec(remotePath, outputs);
		log.append("Output:\n");
		log.append(outputs[0]);
		log.append("Error Messages (if any):\n");
		log.append(outputs[1]);
		if ((exitCode != 0) || (outputs[1].length() > 0)) {
			throw new IOException("The job runner generated an error.\n" +
			"The remote job has not started up successfully.");
		}
		log.append("Done starting job on server.\n");
		log.append("The GUI will perirodically check on the job status.\n");
		setStepStatus(5, false);
	}

	/**
	 * This is a helper method that is used to deploy the job runner
	 * script to the server.
	 * 
	 * This method is invoked from the run() method to adapt the local
	 * runner script (which is a template) by filling in the template
	 * parameters with actual values and copying the script to the
	 * remote machine.
	 * 
	 * @throws Exception This method throws exceptions on errors.
	 */
	private void deployRunnerScript() throws Exception {
		log.append("Generating job runner script...\n");
		setStepStatus(4, true);
		Job job = getJob();
		Server srvr = Workspace.get().getServerList().getServer(job.getServerID());
		// Now compute the PEACE executable path and arguments.
		final String cmdLineParams = job.toCmdLine();
		final String exePath       = srvr.getExecutable(job.getType());
		final String launcherPath  = srvr.getExecutable(EXEKind.WIN_LAUNCHER);
		// Setup the path to pick up the LOCAL job runner template script
		final boolean isWindows    = server.getOSType().equals(OSType.WINDOWS);
		final String jobRunnerPath = "installFiles/" + (isWindows ? "windows/" : "linux/");
		final String jobRunnerFile = "jobRunner"    + (isWindows ? ".bat" : ".sh");
		
		// Get the path to the job on the remote server via helper method.
		final String jobPath = job.getPath();
		// Extract drive information like C:\ or D:\ if job path has it.
		String drive = "";
		if ((jobPath.length() > 3) && (jobPath.substring(1, 3).equals(":\\"))) {
			drive = jobPath.substring(0, 2);
		}
		// Get EAST job path if this job is a clustering job that is supposed
		// to kick off an EAST job after clustering is done.
		String eastWorkPath = "";
		if (job.getType().equals(JobType.CLUSTERING)) {
			eastWorkPath = getEASTWorkDir();
		}
		// Load the template runner script and substitute the pre-defined
		// variables with the information.
		String runnerScript = Utilities.readSmallTextFile(jobRunnerPath + jobRunnerFile);
		runnerScript = runnerScript.replace("%workDir%",   jobPath);
		runnerScript = runnerScript.replace("%eastWorkDir%", eastWorkPath);
		runnerScript = runnerScript.replace("%workDrive%", drive);
		runnerScript = runnerScript.replace("%exe%",       exePath);
		runnerScript = runnerScript.replace("%launcher%",  launcherPath);
		runnerScript = runnerScript.replace("%cmdLine%",   cmdLineParams);
		runnerScript = runnerScript.replace("%nodes%", "" + job.getNodes());
		runnerScript = runnerScript.replace("%cpusPerNode%", "" + job.getCPUsPerNode());
		runnerScript = runnerScript.replace("%memory%", "" + job.getMemory());
		runnerScript = runnerScript.replace("%maxRunTime%", "" + job.getMaxRunTime());
		runnerScript = runnerScript.replace("%peaceVersion%", Version.GUI_VERSION.replace('\n', ' '));
		runnerScript = runnerScript.replace("%jobID%", "PEACE_" + job.getJobID());
		
		log.append("Command line arguments: " + cmdLineParams + "\n");
		// Now copy the runnerScript to a remote file.
		String remotePath = job.getPath() + "/" + jobRunnerFile;
		log.append("Copying job runner script to: " + remotePath + "\n");
		// Create a string input stream to copy data.
		byte[] bytes = runnerScript.getBytes("UTF-8");
		ByteArrayInputStream is = new ByteArrayInputStream(bytes);
		server.copy(is, job.getPath(), jobRunnerFile, 0700);
		log.append("Done copying job runner script.\n");
		setStepStatus(4, false);
	}

	/**
	 * Helper method to copy EST data file to shared space if needed.
	 * This method checks to see if the EST data file already exists
	 * on the server. If not it copies the EST data file over.
	 */
	private void copyESTFile() throws IOException {
		setStepStatus(2, true);
		// Do some basic validation on local EST file.
		String localFileName = getDataSet().getPath();
		log.append("Checking local EST FASTA file: " + localFileName + "\n");
		File localFile = new File(localFileName);
		if (!localFile.exists() || !localFile.isFile() || 
				!localFile.canRead()) {
			throw new IOException("Unable to read local EST file: " + localFileName);
		}
		log.append("Local EST file exists and is readable.\n");
		// Obtain information about the corresponding remote file.
		// Construct server-specific EST file location.
		String serverFileName = getServerESTFile();
		log.append("Checking information on " + serverFileName + "\n");
		FileInfo serverFile = server.fstat(serverFileName);
		if (serverFile.exists()) {
			// The file seems to exist on the server. Do some 
			// additional validation.
			log.append("The entry exists on the server. Validating...\n");
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
				log.append("The remote EST FASTA file seems consistent. Skipping file copy");
				setStepStatus(2, false);
				return;
			}
			// There is an inconsistency in the file. Let user decide if
			// the file should be over written.
			log.append(ioe.getLocalizedMessage());
			JPanel msg = Utilities.collapsedMessage(FILE_MISMATCH, 
					Utilities.toString(ioe));
			int decision = 
				JOptionPane.showConfirmDialog(dialog, msg, 
						"EST File Mismatch", JOptionPane.YES_NO_OPTION);
			if (decision == JOptionPane.NO_OPTION) {
				log.append("User decided not to proceed.");
				throw ioe;
			}
		} else {
			log.append("EST FASTA file is not present on server. Must copy it.\n");
		}
		// When control drops here that means we need to copy the
		// local file to the server.
		String msg = "Copying EST Data to Server.\n";
		log.append(msg);
		ProgressMonitorInputStream pmis =
			new ProgressMonitorInputStream(this, msg, 
					new FileInputStream(localFile));
		pmis.getProgressMonitor().setNote("File: " + localFile.getAbsolutePath());
		// Let the server do the actual copying.
		log.append("Copying file: " + localFile.getAbsolutePath() +
				" to " + serverFileName + "\n");
		server.copy(pmis, "", serverFileName, 0600);
		// Done copying.
		log.append("Done copying.\n");
		setStepStatus(2, false);
	}

	/**
	 * Helper method to compute the full path to the EST data file
	 * on the server.
	 * 
	 * @return The full path (based on install path) to the EST data
	 * file on the server.
	 */
	private String getServerESTFile() {
		String localFileName = getDataSet().getPath();
		File localFile = new File(localFileName);
		// Obtain information about the corresponding remote file.
		Job job = getJob();
		Server srvr = Workspace.get().getServerList().getServer(job.getServerID());
		// Construct server-specific EST file location.
		String serverESTFileName = srvr.getInstallPath() + 
		"/estData/" + localFile.getName();
		return serverESTFileName;
	}

	/**
	 * Helper method to create entries in the work space. This
	 * method is invoked from the run() method and performs
	 * the task of creating MST, Cluster, and Job entries.
	 */
	private void createEntries() {
		// First add necessary work space entries.
		setStepStatus(0, true);
		log.append("Creating workspace entries...\n");
		// Add the MST, cluster, and job entries to the work space.
		Workspace workspace = Workspace.get();
		// First create a full job entry for this job.
		Job job = getJob();
		// Add all the entries to the work space.
		workspace.getJobList().insert(job);
		log.append("\tAdded job entry to workspace.\n");

		// Get the generated files list for this job.
		GeneratedFileList gfl = getGeneratedFiles();
		DataSet dataSet = getDataSet();
		// Add GFL before any existing GFLs that are dependent on
		// this job.
		dataSet.insert(gfl);
		log.append("\tAdded generated files to workspace.\n");
		log.append("Done.\n");
		setStepStatus(0, false);
	}

	/**
	 * This is a helper method that is used to establish a valid
	 * session with the remote server on which the job is to be
	 * submitted.
	 */
	private void startSession() throws IOException {
		Job job = getJob();
		Server srvr = Workspace.get().getServerList().getServer(job.getServerID());
		log.append("Attempting to connect to server " + srvr.getName() + "\n");
		server = SessionFactory.createSession(dialog, srvr);
		// Setup purpose so the user knows why we are connecting
		server.setPurpose("Attempting to start a new job");
		// Connect to the server.
		server.connect();
		// Done. update status.
		log.append("Successfully connected.\n");
	}

	/**
	 * This is a helper method that is used to create a valid
	 * directory where this job's files will be copied.
	 */
	private void createJobDirectory() throws Exception {
		setStepStatus(1, true);
		Job job = getJob();
		Server srvr = Workspace.get().getServerList().getServer(job.getServerID());
		log.append("Attempting to create directory for job on server...\n");
		// repeatedly try and create a work folder for the job based on its id
		int tryCount = 0;
		do {
			// Compute jobDir based on the try count
			String jobDir = srvr.getInstallPath() + "/jobs/" + job.getJobID();
			if (tryCount > 0) {
				jobDir += "_" + tryCount;
			}
			// Check if the directory exists. If not create a new one
			// and use it.
			log.append("Checking if directory " + jobDir + " exists...");
			FileInfo fi = server.fstat(jobDir);
			if (!fi.exists()) {
				// OK, the directory does not exist
				log.append("No\n");
				log.append("Creating directory " + jobDir + " ...");
				server.mkdir(jobDir);
				log.append("Done.\n");
				job.setPath(jobDir);
				log.append("Set job path to " + jobDir + "\n");
				// Break out of the while loop
				break;
			} else {
				log.append("Yes.\nTrying a different directory next.\n");
				tryCount++;
			}
		} while (tryCount < 5);
		// When control drops here we decide if we are going to
		// abort or not.
		if (tryCount >= 5) {
			// Creating work directory failed. Bail out.
			throw new IOException("Unable to create a useable work directory for job");
		}
		setStepStatus(1, false);
	}

	/**
	 * Flag to indicate if the job must be started after
	 * it is created. If the flag is true, then the job is started.
	 * Otherwise the job is not started (as it may be dependent on another
	 * job to finish). This flag is used in cases where a clustering job
	 * is run first followed by an assembly job.
	 */
	private final boolean startJob;
	
	/**
	 * Method to obtain the job entry that is being created by this wizard.
	 * 
	 * This method is expected to be appropriately overridden by the
	 * derived class to return the job object on which this submission
	 * wizard is currently operating.
	 * 
	 * @return The job object on which this wizard is currently operating.
	 * The return value cannot be null.
	 */
	protected abstract Job getJob();
	
	/**
	 * Method to obtain the list of files that are being generated 
	 * by this wizard.
	 * 
	 * This method is expected to be appropriately overridden by the
	 * derived class to return the list of generated files to be
	 * created by this job..
	 * 
	 * @return The list of generated files to be created by the job
	 * being submitted. The return value cannot be null.
	 */
	protected abstract GeneratedFileList getGeneratedFiles();
	
	/**
	 * Helper method to get various other (non-data-set) files
	 * to be copied to the server for the job being submitted.
	 * 
	 * This method is expected to be appropriately overridden by the
	 * derived class to return a list of additional file entries
	 * (from the current work space) to be copied to the job
	 * directory on the server.
	 * 
	 * @return The list of additional files to be copied to the server
	 * by this job. The return value of this method can be null.
	 */
	protected abstract FileEntry[] getFilesToCopy();
	
	/**
	 * Call-back method to indicate that the thread for creating/submitting
	 * a job is done.
	 * 
	 * This method can be overridden by derived classes to perform any 
	 * additional/special operation that may be needed once a job 
	 * submission has completed.
	 */
	protected void done() {
		// Intentionally nothing to be done here.
	}
	
	/**
	 * Method to obtain the data set entry that is being 
	 * used by the job being submitted.
	 * 
	 * This method is expected to be appropriately overridden by the
	 * derived class to return the data set object on which this submission
	 * wizard is currently operating.
	 * 
	 * @return The data set object on which this wizard is 
	 * currently operating. The return value cannot be null.
	 */
	protected abstract DataSet getDataSet();
	
	/**
	 * Method to obtain the path to a dependent EAST job when a 
	 * clustering job is run.
	 * 
	 * This method is expected to be appropriately overridden by the
	 * derived class. This method must return the working directory
	 * (on the server) of the EAST job that must be started once the
	 * clustering job completes successfully. If a dependent job is
	 * not present then this method must return an empty string.
	 * 
	 * @return The work directory of an EAST job that must be started
	 * after a clustering job ends. The return value cannot be null
	 * but can be an empty string if a dependent job is absent.
	 */
	protected abstract String getEASTWorkDir();
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final WizardDialog dialog;

	/**
	 * Field to read and edit a brief description about the job.
	 * This information can be anything the user desires and is
	 * meaningful only to the user.
	 */
	private JTextArea log;

	/**
	 * The step information labels to display information about
	 * the steps to be done and those that have been completed.
	 * These labels are created in the constructor and updated
	 * periodically.
	 */
	private JLabel stepInfo[];

	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some contextual information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"<html>Please wait while the wizard verifies and submits<br>" +
		"the job on the remote server...</html>";

	/**
	 * A generic informational message that is displayed to the user
	 * when any error occurs during the job submission process.
	 */
	private static final String JOB_SUBMIT_ERROR = 
		"<html>Error occured when attempting to submit the job on<br>" +
		"the remote server. The job was not successfully started.<br>" +
		"The necessary MST and cluster data will not be computed.</html>";

	/**
	 * A generic informational message that is displayed to the user
	 * when EST file mismatch is detected
	 */
	private static final String FILE_MISMATCH = 
		"<html>Error the remote EST file is not consistent with the local<br>" +
		"copy of the EST file. Do you want to overwrite the copy on the server?</html>";

	/**
	 * A generic informational message that is displayed to the user
	 * after a job has successfully started up.
	 */
	private static final String SUCCESS_MESSAGE = 
		"<html>The job was successfully started on the server.<br>" +
		"The GUI will periodically check on the status of the job<br>" +
		"and notify you about the progress. You do not have to keep<br>" +
		"the GUI open until the job completes. You may exit the GUI<br>" +
		"and run it again at a later time to check on the status of the job.</html>";

	/**
	 * The session to the server on which the job is to be submitted.
	 * The connection is used throughout the process of job submission
	 * by various methods. It is created in the startSession() method.
	 */
	private ServerSession server;

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8643291107936550834L;
}
