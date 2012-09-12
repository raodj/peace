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
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.MainFrame;
import org.peace_tools.core.job.DataSetSelectionPage;
import org.peace_tools.core.job.ServerComboBox;
import org.peace_tools.decagon.DecagonVariables;
import org.peace_tools.decagon.helpers.DADXHelper;
import org.peace_tools.decagon.helpers.DADXSummary;
import org.peace_tools.decagon.jaxb.AssemblerDetails;
import org.peace_tools.decagon.jaxb.Job;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as the top-level class for creating/submitting
 * a new assembly job using information from a given DECAGON Assembler
 * Description XML (DADX) file.
 * 
 * In DECAGON, each genomic assembler that can be used for 
 * assembling cDNA fragments must have a suitable DADX file that
 * provides all the necessary information. This class serves as a
 * top-level class that uses a DADX file to run an assembler on
 * a given cDNA file.
 * 
 * <ul>
 * 
 * <li>This class must be instantiated using a DADXHelper class
 * that encapsulates the information about the assembler.</li>
 * 
 * </ul>
 * 
 * This top-level class merely creates the various pages and adds them
 * to this wizard. Each page performs a specific task required to 
 * create/configure a complete job(s) to be run on a server.
 *  
 * <p><b>Note:</b> In order to run an assembler, there must be at 
 * least one data set and one server entry in the work space.</p> 
 */
public class DecagonJobCreator extends WizardDialog {
	/**
	 * Helper method to validate necessary preconditions prior to creating the
	 * job wizard.
	 * 
	 * This is a helper method that is used to ensure that:
	 * 
	 * <ol>
	 * <li>That there is at least one data set in the workspace.</li>
	 * <li>That there is at least one valid server entry that can be used for
	 * assembly.</li>
	 * </ol>
	 * 
	 * @param parent The parent window for the wizard dialog box. This parameter
	 * cannot be null.
	 * 
	 * @return If the workspace has the necessary/valid data sets and server
	 * entries, then this method instantiates a valid wizard and returns
	 * the newly created wizard for further use. 
	 */
	public static DecagonJobCreator create(DADXHelper dadxHelper, 
			MainFrame parent) {
		// Check to ensure we have a valid data set that contains a 
		// FASTA file for assembly.
		String[] msg = DecagonJobCreator.NO_DATASET_MSG; // Set invalid condition
		final ArrayList<DataSet> dataSets = Workspace.get().getDataSets();
		// Check to ensure we have at least one good data set entry
		for(DataSet ds: dataSets) {
			if (ds.isGood() && (ds.getFileType().equals(DataFileType.FASTA))) {
				msg = null; // We got a good data set
				break;      // no further checks needed.
			}
		}
		// Next check to ensure that we have a server that
		// can do clustering and/or assembly.
		if ((msg == null) && 
			(ServerComboBox.getServerList(true, false, new JComboBox()) == 0)) {
			// We have valid data set/mst but no valid server
			msg = DecagonJobCreator.NO_SERVER_MSG;
		}
		// Display error message if it is set
		if (msg != null) {
			WizardDialog.showMessage(parent, msg);
			return null;
		}
		// Return new wizard only if msg is null indicating no errors.
		return new DecagonJobCreator(dadxHelper, parent);
	}
	
	/**
	 * The constructor for the DECAGON Job Creator.
	 * 
	 * The constructor lays out the wizard and creates all the
	 * wizard pages. The wizard pages are created as templates
	 * and get populated with the necessary information just before
	 * the pages get displayed to the user.
	 * 
	 * @param dadxHelper The DADX helper object that contains the 
	 * unmarshalled DADX description for the assembler pipeline for
	 * which job(s) are to be created by this wizard.
	 * 
	 * @param parent The main frame that logically owns this wizard.
	 * 
	 * @see DecagonJobCreator#create(String, MainFrame, boolean)
	 */
	private DecagonJobCreator(DADXHelper dadxHelper, MainFrame parent) {
		super(parent);
		this.dadxHelper = dadxHelper;
		this.mainFrame  = parent;
		
		setTitle("DECAGON Job Creator");
		setResizable(false);
		// Set up the title image we want to use.
		setTitleBackground("images/peace_wizard_header.png", Color.white);
		// Set up the column image we want to use.
		setSequenceBackground("images/peace_wizard_column.png");
		// First setup the overview page.
		addPage(createOverview(dadxHelper));
		// Next add the page to select data set to be used for assembly.
		dssp = new DataSetSelectionPage(this, true);
		addPage(dssp);
		// Next create a server selection page.
		ssp = new ServerSelectionPage(this);
		addPage(ssp);
		// Once the user selects a server we need to ensure that
		// the executables for the pipeline are operational on the machine.
		exeCheck = new CheckExecutables(this, ssp);
		addPage(exeCheck);
		// Next let the user set the parameters for the assembler.
		ipp = new InputParamsPage(this);
		addPage(ipp);
		// Check and create server configuration for each job entry.
		List<Job> jobList = this.dadxHelper.getAssemblerDetails().getJob();
		jobSrvrConf  = new JobConfigPage[jobList.size()];
		for(int i = 0; (i < jobList.size()); i++) {
			jobSrvrConf[i] = new JobConfigPage(this, jobList.get(i), ssp);
			addPage(jobSrvrConf[i]);
		}
		// Next add wizard-page to obtain additional information from
		// the user.
		aip = new AdditionalInfoPage(this);
		addPage(aip);
		// Next add the verification wizard page that displays the 
		// various variables.
		VerificationPage vp = new VerificationPage(this);
		addPage(vp);
		// Finally add the page that performs the job submission process
		jsp = new JobSubmissionPage(this);
		addPage(jsp);
	}
	
	/**
	 * Convenience method to update various DECAGON variables.
	 * 
	 * This method provides a convenience interface to the 
	 * wizard pages to update the DECAGON variables that
	 * can be used in a DADX file. This method sets up the 
	 * following set of variables:
	 * 
	 * <ul>
	 * <li>It sets up various variables associated with the data
	 * set and its properties via call to 
	 * {@link DecagonVariables#setVariables(DataSet)}.</li>
	 * 
	 * <li>It sets up DECAGON variables associated with the target
	 * machine via call to {@link DecagonVariables#setVariables(Server)}</li>
	 * 
	 * </ul>
	 */
	protected void updateDecagonVariables() {
		dadxHelper.getVariables().setVariables(getDataSet());
		dadxHelper.getVariables().setVariables(ssp.getSelectedServer());
	}
	
	/**
	 * Helper method to create an overview page.
	 * 
	 * This is a helper method that is invoked only once from
	 * the constructor to create an overview page. This page
	 * is a simple wizard page that contains an overview message
	 * ({@link #ASSEMBLY_OVERVIEW_MSG}) and a brief description 
	 * about the assembler obtained from the DADX file.
	 * 
	 * @param dadxHelper The helper class that encapsulates the
	 * unmarshalled information about the assembler to be run.
	 * Information about the assembler is extracted from this
	 * object. 
	 * 
	 * @return A generic wizard page that contains an overview
	 * message along with information about the genomic
	 * assembler.
	 */
	private GenericWizardPage createOverview(DADXHelper dadxHelper) {
		// Shortcut to assembler details to streamline code...
		AssemblerDetails asmInfo = dadxHelper.getAssemblerDetails();
		// Create a JLabel that contains some of the static information
		// to be displayed to the user in the overview page.
		final String msg  = String.format(ASSEMBLY_OVERVIEW_MSG, 
				asmInfo.getName(), asmInfo.getJob().size());
		JLabel genericMsg = new JLabel(msg);
		// Create a label with the assembler name and icon to make the display
		// look nice.
		DADXSummary summary = DADXSummary.getSummary(asmInfo.getName());
		JLabel asmTitle     = new JLabel("About " + asmInfo.getName(), 
				summary.getSmallIcon(), JLabel.LEFT);
		// Create a panel to hold the generic message and title.
		JPanel topPanel = new JPanel(new BorderLayout(5, 5));
		topPanel.add(genericMsg, BorderLayout.NORTH);
		topPanel.add(asmTitle, BorderLayout.SOUTH);
		
		// Place description of assembler into a HTML pane. If that fails
		// create a simple label with description about assembler.
		Color bgColor = mainFrame.getBackground();
		bgColor = new Color(bgColor.getRed(), bgColor.getGreen(), bgColor.getBlue() + 0xf);
		JComponent helpView = mainFrame.createHTMLComponent(asmInfo.getDescription(), bgColor);
		// Finally create a generic wizard page with the various GUI
		// elements in it.
		GenericWizardPage overview = new GenericWizardPage();
		overview.setTitle("Overview", "Overview of this wizard.");
		overview.setBorder(new EmptyBorder(20, 15, 10, 5));
		overview.add(topPanel, BorderLayout.NORTH);
		overview.add(helpView, BorderLayout.CENTER);
		
		// Return the overview page to be used in the main wizard
		return overview;
	}
	
	/**
	 * Helper method invoked when user clicks cancel button.
	 * 
	 * This is a helper method that is overridden in this class.
	 * This method is invoked when the user clicks the cancel
	 * button in the wizard. This method is used to display a 
	 * confirmation dialog to ensure that the user really wants to
	 * exit out of the wizard.
	 * 
	 * @return This method returns true if the user wants to quit
	 * out of the wizard dialog.
	 */
	@Override
	protected boolean cancel() {
		// Check if user really want's to quit.
		int result = JOptionPane.showConfirmDialog(this,
				"Are you sure you want to exit from this wizard?",
				"Confirm", JOptionPane.YES_NO_OPTION);
		if (result == JOptionPane.NO_OPTION) {
			// The user does not want to quit.
			return false;
		}
		// Yes, the user wants to quit. Ensure no threads are left.
		return super.cancel();
	}
	
	/**
	 * This method overrides the final notification method in this
	 * wizard. Currently this method has no specific task to 
	 * perform but is present as a place holder for future 
	 * enhancements. 
	 */
	@Override
	public void done(boolean success) {
		if (!success) {
			return;
		}
	}
	
	/**
	 * Helper method to create a temporary or permanent clustering or
	 * assembly job.
	 * 
	 * This is an internal helper method that is used by this wizard
	 * to display summary information and create a job entry to be
	 * added to the workspace.
	 * 
	 * @param reserveJobID If this flag is true, then this method
	 * requests a new job ID to be created and reserved. Otherwise
	 * the job ID is set to the string "TBD".
	 * 
	 * @param assembly If this flag is true, then this method creates
	 * a job entry for assembly via EAST. If this flag is false
	 * this method creates a job entry for a parallel clustering job
	 * via PEACE.
	 * 
	 * @return A job entry for a job corresponding to a clustering or
	 * assembly job.
	 */
	protected Job createJobEntry(boolean assembly, boolean reserveJobID) {
		/*
		final int jobIndex = (assembly ? 1 : 0);
		if (job[jobIndex] != null) {
			return job[jobIndex];
		}

		// Get parameters for assembly job
		ArrayList<Param> paramList = (assembly ? ebpwp.getParamList() :
			new ArrayList<Param>());
		// Get the server information provided by the user
		ServerWizardPage srvrPage = (assembly ? assemblySWP : clustSWP);
		Server server = srvrPage.getServerInfoPanel().getSelectedServer(); 
		// Get the platform-configuration for this job
		final int platformInfo[] = srvrPage.getServerInfoPanel().getPlatformConfiguration();
		// Compute the total memory to be used for the job by:
		// #Nodes * #CPUs/Node * #Memory/CPU
		final int totalMemory = platformInfo[0] * platformInfo[1] * platformInfo[2];
		// Setup job description
		final String rawJobDesc = (assembly ? owp.getDescription() : DEFAULT_DESCRIPTION);
		final String jobDesc    = DOMHelper.xmlEncode(rawJobDesc);
		// Get an unique ID for this job if requested
		JobList jobList = Workspace.get().getJobList();
		final String jobID = (reserveJobID ? jobList.reserveJobID() : "TBD");
		// Create the new job entry (the following variable is set by the
		// if-else block below)
		Job tmpJob = null;
		if (!assembly) {
			// Create default analyzer to be used for comparing cDNA fragments
			// for constructing MST.
			FWAnalyzer analyzer = new FWAnalyzer(FWAnalyzerType.TWOPASSD2, true, 100, 6, "heap", 128);
			// Create entries for default set of filters used by PEACE
			ArrayList<Filter> filtList = (assembly ? null : ClusteringJob.createDefaultFilters());
			// Create no heuristics as we use an adaptive two-pass-d2
			// analyzer that uses its own custom heuristics and does not
			// need any heuristics specified on the command line.
			ArrayList<Heuristic> heurList = null;
			// Create clustering job
			ClusteringJob clsJob = 
				new ClusteringJob(
					jobID, jobDesc, server.getID(),
					null,             // Path on server for job-files.
					platformInfo[0],  // number of nodes, 
					platformInfo[1],  // CPUs per node, 
					totalMemory,      // memory, 
					platformInfo[3],  // runTime
					analyzer,
					1,                // Threshold for clustering
					heurList,         // heuristic list
					filtList,         // Filter list
					paramList);       // Empty parameters for starters
			// Ensure parameters are populated
			GeneratedFileList gfl = createGFL(false, reserveJobID, clsJob);
			clsJob.setupParameters(getDataSet(), gfl, true);
			// Setup the common job variable for further use
			tmpJob = clsJob;
		} else {
			// Assembly job
			EASTJob eastJob = new EASTJob(jobID, jobDesc, server.getID(), 
					null,             // Path on server for job-files.
					totalMemory,      // memory, 
					platformInfo[3],  // runTime
					paramList);       // Cmd-line parameters
			// Add input & output files as parameters as well. For this we
			// need to figure out the MST file entry depending on whether
			// this is a clustering+assembly or just assembly job.
			FileEntry mstFE = getMSTFile();
			if (mstFE == null) {
				// This is a clustering+assembly job.
				mstFE = new FileEntry("TBD", FileEntryType.MST, 
						DataFileType.TXT, ofwp.getFilePath(1), DEFAULT_DESCRIPTION);
			}
			GeneratedFileList gfl = createGFL(true, reserveJobID, eastJob);
			eastJob.addParameters(getDataSet(), mstFE, null, gfl);
			// Setup the common job variable for further use
			tmpJob = eastJob;
		}
		assert ( tmpJob != null );
		// Save job entry for future reference only if we have reserved
		// a job ID.
		if (reserveJobID) {
			job[jobIndex] = tmpJob;
		}
		// Update the dependent job ID & summary on the EAST job.
		if ((job[0] != null) && (job[1] != null)) {
			job[1].setPreviousJobID(job[0].getJobID());
			// Update status of any job summaries in the data sets
			GeneratedFileList eastGfl = Workspace.get().getGFL(job[1].getJobID());
			if (eastGfl != null) {
				eastGfl.getJobSummary().setPreviousJobID(job[0].getJobID());
			}
		}
		// Return the newly constructed job object back to the caller
		return tmpJob;
		*/
		return null;
	}

	/**
	 * Access method to obtain the list of files to be generated by
	 * this method. This method is typically used by the
	 * SubmitClusteringJobWizard class.
	 * 
	 * @return The list of files generated by the job created by
	 * this wizard.
	 */
	protected GeneratedFileList createGFL(boolean assemblyJob, 
			boolean reserveID, Job job) {
		return null;
		/*
		// Build a job summary.
		JobSummary summary = new JobSummary(job);
		Workspace workspace = Workspace.get();
		GeneratedFileList gfl = new GeneratedFileList(summary);

		if (!assemblyJob) {
			final String mstID = (reserveID ? workspace.reserveID() : "TBD");
			FileEntry mstFE = new FileEntry(mstID, FileEntryType.MST, 
					DataFileType.TXT, ofwp.getFilePath(1), DEFAULT_DESCRIPTION);
			// Append information about cluster file by creating a 
			// temporary dummy entry.
			final String clsID = (reserveID ? workspace.reserveID() : "TBD");
			FileEntry clsFE = new FileEntry(clsID, FileEntryType.CLS,
					DataFileType.TXT, ofwp.getFilePath(2), DEFAULT_DESCRIPTION);
			// Add jobs to the generated file list
			gfl.add(mstFE);
			gfl.add(clsFE);
		} else {
			final String description = owp.getDescription();
			final String contigID = (reserveID ? workspace.reserveID() : "TBD");
			FileEntry contigFE = new FileEntry(contigID, FileEntryType.ASM, 
					ebpwp.getContigOutputFormat(), ofwp.getFilePath(3), description);
			
			final String singID = (reserveID ? workspace.reserveID() : "TBD");
			FileEntry singFE = new FileEntry(singID, FileEntryType.SINGLETONS, 
					DataFileType.FASTA, ofwp.getFilePath(4), description);
			
			final String statsID = (reserveID ? workspace.reserveID() : "TBD");
			FileEntry statsFE = new FileEntry(statsID, FileEntryType.STATS, 
					DataFileType.TXT, ofwp.getFilePath(5), description);
			
			// Add jobs to the generated file list
			gfl.add(contigFE);
			gfl.add(singFE);
			gfl.add(statsFE);
		}
		return gfl;
		*/
	}
	

	/**
	 * Helper method to return the data set on which
	 * operations are being performed by this wizard.
	 * 
	 * This method returns the the data set entry which the user
	 * has chosen to be assembled.
	 * 
	 * @return The path to the main data set file on which the user is
	 * operating on. If a valid data set file is not yet available,
	 * then this method returns null.
	 */
	protected DataSet getDataSet() {
		DataSet dataSet      = null;
		Object currSelection = dssp.getSelection();
		if (currSelection instanceof DataSet) {
			dataSet = (DataSet) currSelection;
		} 
		return dataSet;
	}
	
	/**
	 * Convenience method to obtain the server on which jobs are to
	 * be run.
	 * 
	 * This method can be used to obtain the server on which the jobs
	 * for an assembler-pipeline are to be run. This method
	 * returns the server selected by the user via call to
	 * {@link ServerSelectionPage#getSelectedServer()} method.
	 *  
	 * @return The server entry currently selected by the user for
	 * running jobs.
	 */
	protected Server getSelectedServer() {
		return ssp.getSelectedServer();
	}
	
	/**
	 * Obtain the instance of the main frame class that owns this wizard.
	 * 
	 * This method is currently used by the submit job wizard page 
	 * to create a job monitor once a job has been successfully
	 * submitted. 
	 * 
	 * @return The main frame that owns this wizard.
	 */
	protected MainFrame getMainFrame() {
		return mainFrame;
	}
	
	
	/**
	 * Obtain the DADX helper object associated with this wizard.
	 * 
	 * This is a convenience method that is used by various pages
	 * in this wizard to obtain a reference to the DADX helper 
	 * object being currently used.
	 * 
	 * @return the dadxHelper The DADX helper class associated
	 * with this wizard. The return value is never null.
	 */
	protected DADXHelper getDADXHelper() {
		return dadxHelper;
	}

	/**
	 * Obtain the list of wizard pages that contain job configuration
	 * information.
	 * 
	 * For each Job associated with this assembler pipeline this 
	 * method creates a Job configuration page. The job configuration
	 * page provides the following information about each job:
	 * server to be used, number of nodes to be used, CPUs per node,
	 * memory to be reserved for the job, and job run time set by
	 * the user.
	 * 
	 * @return The list of wizard pages (in this wizard) that contain
	 * job configuration information.
	 */
	protected JobConfigPage[] getJobConfigs() {
		return jobSrvrConf;
	}
	
	/**
	 * Obtain the description entered by the user for this pipeline run.
	 * 
	 * This method can be used to obtain the description entered by the
	 * user for this pipeline run. This method returns correct information
	 * only after the corresponding wizard page has been displayed and
	 * the user has entered the description.
	 * 
	 * @return The description entered by the user.
	 */
	protected String getDescription() {
		return aip.getDescription();
	}
	
	/**
	 * Obtain the directory on local machine where generated files
	 * are to be copied.
	 * 
	 * This method must be used by various wizard pages to determine
	 * the directory where the various files generated by running the
	 * DADX pipeline are to be stored. This method returns correct information
	 * only after the corresponding wizard page has been displayed and
	 * the user has entered the description.
	 * 
	 * @return This directory on the local machine where the generated
	 * output files are to be copied. 
	 */
	protected String getStorageDir() {
		return aip.getStorageDir();
	}
	
	/**
	 * The DADX helper object that contains information about the
	 * assembler's software pipeline. This information is used
	 * by the various wizard pages to perform the operations
	 * necessary to create and submit a suitable DECAGON job.
	 */
	private final DADXHelper dadxHelper;

	/**
	 * The wizard page that contains the description of the data
	 * set selected for the job along with a brief job description.
	 */
	private final DataSetSelectionPage dssp;
	
	/**
	 * The wizard page that permits the user to select the server
	 * to be used for running the various jobs in the pipeline
	 * (associated with the selected assembler/DADX file).
	 */
	private final ServerSelectionPage ssp;
	
	/**
	 * The set of wizard pages that contain the server and associated
	 * job configuration information for each job.
	 */
	private final JobConfigPage[] jobSrvrConf;

	/**
	 * The wizard page that checks to ensure that the various
	 * executables referenced in the assembler configuration
	 * (in {@link #dadxHelper}) are operational on the server
	 * selected by the user. This page is run after the job
	 * configurations are obtained.
	 */
	private final CheckExecutables exeCheck;
	
	/**
	 * The wizard page that obtains input from the users for
	 * various parameters associated with the given DADX 
	 * assembler pipeline. This page is displayed to the user
	 * before the individual job configurations are obtained.
	 */
	private final InputParamsPage ipp;

	/**
	 * The wizard page that obtains some additional information
	 * from the user regarding the use of the given DADX
	 * assembler pipeline. This page is displayed as the
	 * last page to obtain user inputs.
	 */
	private final AdditionalInfoPage aip;

	/**
	 * The wizard page that performs the actual task of creating
	 * DECAGON jobs on the server. This page is the last page
	 * on the wizard. The process of submitting jobs on the server
	 * is run as a background thread with the GUI continuously
	 * displaying updated status
	 */
	private final JobSubmissionPage jsp;

	/**
	 * The main frame that logically owns this job. This value is
	 * used to create a job monitor if a job is successfully 
	 * submitted to run on a server.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * A simple static message that is displayed to the user if
	 * a valid data set is not found in the workspace for a
	 * clustering + assembly job.  
	 */
	private static final String[] NO_DATASET_MSG = {
		"The cDNA/EST fragments (called a <i>dataset</i>) to be<br/>" +
		"clustered and assembled must be present in the workspace.",
		"<b>Currently, the workspace does not have a dataset with<br/>" +
		"a FASTA file (EAST can work only with FASTA files)</b>",
		"Add the cDNA data set to be clustered and assembled<br/>" +
		"to the workspace prior to launching this wizard.<br/>" +
		"You may add a new data set using the data set wizard."
	};
	
	/**
	 * A simple static message that is displayed to the user if
	 * a valid MST file is not found in the workspace for a
	 * assembly-only job.  
	 */
	private static final String[] NO_SERVER_MSG = {
		"In order to run an assembly job a server in good operating<br/>" +
		"condition must be present in the workspace.",
		"<b>Could not find a valid server in good operating condition<br/>" +
		"that has both PEACE and EAST installed on it.</b>",
		"This issue can resolved either by:<br/>" +
		"&nbsp; &#149; Adding a new server entry using the server wizard<br/>" +
		"&nbsp; &#149; or by troubleshooting an existing server entry, if any<br/>" +
		"</ul>"
	};
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user.
	 */
	private static final String ASSEMBLY_OVERVIEW_MSG = "<html>"+
		"This wizard guides you through the process of assembling cDNA<br/>" +
		"reads from a given DataSet using %s<br/><br/>" + 
		"The assembly process will consist of %d job(s) to be run on a server<br/>" +
		"you select in this wizard. An outline of steps to be performed by<br/>" + 
		"this wizard is shown to the left. Additional details about the assembler<br/>" +
		"is available below." +
		"</html>";
	
	/**
	 * The generated serialization UID (need to keep the compiler happy) 
	 */
	private static final long serialVersionUID = -7424095401537856423L;
}
