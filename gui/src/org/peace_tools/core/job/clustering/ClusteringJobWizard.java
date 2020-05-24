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

package org.peace_tools.core.job.clustering;

import java.awt.BorderLayout;
import java.awt.Color;
import java.util.ArrayList;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.MainFrame;
import org.peace_tools.core.job.JobInfoWizardPage;
import org.peace_tools.core.job.ServerPanel;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.ClusteringJob;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FWAnalyzer;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Heuristic;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.JobSummary;
import org.peace_tools.workspace.Param;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as the top-level class for creating a new 
 * job. A job is used to compute MST data and additionally clustering
 * information for a given EST. This top-level class merely
 * creates the various pages and adds them to the wizard. Each
 * page performs a specific task required to create/configure
 * a complete job to be run on a server.
 *  
 * <p><b>Note:</b> In order to create a job, there must be at least one
 * data set and one server entry in the work space.</p> 
 */
public class ClusteringJobWizard extends WizardDialog {
	/**
	 * The constructor for the job wizard.
	 * 
	 * The constructor lays out the wizard and creates all the
	 * wizard pages. The wizard pages are created as templates
	 * and get populated with the necessary information just before
	 * the pages get displayed to the user.
	 * 
	 * @param title The title to be set for the main wizard frame.
	 * @param parent The main frame that logically owns this wizard.
	 */
	public ClusteringJobWizard(String title, MainFrame parent) {
		super(parent);
		this.mainFrame = parent;
		setTitle(title);
		setResizable(false);
		// Set up the title image we want to use.
		setTitleBackground("images/peace_wizard_header.png", Color.white);
		// Set up the column image we want to use.
		setSequenceBackground("images/peace_wizard_column.png");
		// First setup the overview page.
		createOverview();
		// Create page to permit the user to select data set
		jiwp = new JobInfoWizardPage(this);
		addPage(jiwp);
		// Create page to permit the user to choose tree type
		twp = new TreeWizardPage();
		addPage(twp);
		// Create page to permit the user to setup the analyzer params
		awp = new AnalyzerWizardPage(this);
		addPage(awp);
		// Create heuristics information page
		hwp = new HeuristicsWizardPage(this, awp);
		addPage(hwp);
		// Create prime heuristic page
		phwp = new PrimesHeuristicWizardPage();
		addPage(phwp);
		// Create filter information page.
		fwp = new FiltersWizardPage(this);
		addPage(fwp);
		// Create the MST file information file.
		mwp = new MSTWizardPage(this);
		addPage(mwp);
		// Create the cluster file information.
		cwp = new ClusterWizardPage(this);
		addPage(cwp);
		// Create the verification page to display summary
		VerifyWizardPage vwp = new VerifyWizardPage(this);
		addPage(vwp);
		// Finally add the page that does all the hard work
		SubmitClusteringJobWizardPage sjwp = new SubmitClusteringJobWizardPage(this);
		addPage(sjwp);
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
	 * Helper method to create the overview page. This method was
	 * introduced to keep the code clutter in the constructor to
	 * a bare minimum.
	 */
	private void createOverview() {
		JLabel message = new JLabel(OVERVIEW_MSG);
		Utilities.adjustFont(message, 0, 10, -1);
		GenericWizardPage overview = new GenericWizardPage();
		overview.add(message, BorderLayout.NORTH);
		overview.setTitle("Overview", "Overview of tasks in this wizard.");
		overview.setBorder(new EmptyBorder(20, 15, 10, 5));
		addPage(overview);
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
	 * Obtain the job entry created and added by this wizard.
	 * 
	 * The return value from this method is valid only after the
	 * last wizard page has been displayed and it has called
	 * the createWorkspaceEntries() method to actually create
	 * this entry.
	 * 
	 * @return The job entry created by this wizard.
	 */
	protected Job getJob() { return job; }

	/**
	 * Access method to obtain the list of files to be generated by
	 * this method. This method is typically used by the
	 * SubmitClusteringJobWizard class.
	 * 
	 * <p><b>Note:</b> It is important to create a new GFL each
	 * time this method is called (to ensure job summary is 
	 * updated consistently).
	 * PrimesHeuristicWizardPage 
	 * @return The list of files generated by the job created by
	 * this wizard.
	 */
	protected GeneratedFileList createGFL() {
		// Build a job summary.
		JobSummary summary = new JobSummary(job, getAnalyzerType());
			
		Workspace workspace = Workspace.get();
		GeneratedFileList gfl = new GeneratedFileList(summary);
		gfl.add(mwp.getMSTFileEntry(workspace.reserveID(), job.getDescription()));
		gfl.add(cwp.getClusterFileEntry(workspace.reserveID()));

		return gfl;
	}
	
	/**
	 * Helper method to create the job entry to be added to the workspace.
	 * 
	 * This method is a helper method that is used to create the actual
	 * job entry in the workspace. This method is invoked from the 
	 * createWorkspaceEntries() method. This method obtains information
	 * from various wizard pages and uses them to create the job entry.
	 * 
	 * <p><b>Note</b>: This method performs the job creation operation
	 * the first time it is called. After that the same job object
	 * is returned when this method is called.</p>
	 * 
	 * @return The workspace object that represents the new job created
	 * by this method.
	 */
	protected ClusteringJob createJobEntry() {
		if (job != null) {
			return job;
		}
		ServerPanel sp = mwp.getServerInfoPanel();
		// Determine ID of the selected server.
		Server server   = sp.getSelectedServer();
		// Get the nodes and cpus/node..
		int platformInfo[] = sp.getPlatformConfiguration();
		// Get an unique ID for this job.
		JobList   jobList   = Workspace.get().getJobList();
		// Compute the total memory to be used for the job by:
		// #Nodes * #CPUs/Node * #Memory/CPU
		final int totalMemory = platformInfo[0] * platformInfo[1] * platformInfo[2];
		// Get the frame-word analyzer information
		FWAnalyzer analyzer = awp.getAnalyzer();
		// Setup a empty list of parameters
		ArrayList<Param> parameters = new ArrayList<Param>();
		// Setup --maxUse parameter if user chose AST instead of MST
		if(twp.getTreeSelection().equals(TreeWizardPage.AST)) {
			parameters.add(new Param("--maxUse",  twp.getASTThreshold()));
		}
		// Get the heuristics, including primes heuristic
		ArrayList<Heuristic> heuristics = hwp.getHeuristics();
		Heuristic primesHeuristic = phwp.getPrimesHeuristic();
		if(primesHeuristic != null) {
			if(heuristics == null) {
				heuristics = new ArrayList<Heuristic>();
			}
			heuristics.add(primesHeuristic);
		}
		// Now we have all the info to create the new job entry. The
		// only thing critical thing that is not yet finalized is
		// the path where the job files are stored on the remote 
		// server. That will happen once the job is actually submitted.
		job = new ClusteringJob(jobList.reserveJobID(),
				jiwp.getDescription(),
				server.getID(),	null, 
				platformInfo[0],     // number of nodes, 
				platformInfo[1],     // CPUs per node, 
				totalMemory,         // memory, 
				platformInfo[3],     // runTime
				analyzer,            // FWAnalyzer
				cwp.getThreshold(),  // clustering threshold
				heuristics, 		 // heuristic list
				fwp.getFilters(),    // Filter list
				parameters);
		// Ensure parameters are populated
		job.setupParameters(jiwp.getDataSet(), createGFL(), jiwp.isMaksBasesSet());
		// Return the newly created job entry.
		return job;
	}
	
	/**
	 * This is a helper method that is invoked from the last wizard
	 * page (SubmitJobWizardPage) to obtain summary data. This method
	 * creates a job summary string and returns it back to the caller.
	 * This information is used to provide the user with a brief
	 * summary on the job that will be run by the wizard.
	 */
	protected String getSummary() {
		StringBuilder sb = new StringBuilder();
		// Obtain the necessary information into a simple string.
		sb.append("SUMMARY OF JOB\n");
		// Append information about the server and nodes.
		sb.append(mwp.getSummary("\t"));
		// Append information about the analyzer
		FWAnalyzer analyzer = awp.getAnalyzer();
		// Summarize the analyzer job.
		sb.append(analyzer.getSummary());
		// Append information about the heuristics being used.
		sb.append("\n\nHEURISTICS INFORMATION:\n");
		if (this.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2)) {
			sb.append("\tAutomatically configured u/v and t/v heuristics\n");
		} else {
			sb.append(hwp.getHeuristicsSummary());
		}
		if(phwp.getPrimesHeuristic() != null) {
			sb.append(phwp.getPrimesHeuristic().getSummary());
		}
		// Append information about filters being used.
		sb.append("\nFILTER INFORMATION:\n");
		sb.append(fwp.getSummary());
		
		// Append information about the MST and cluster data files
		// by creating temporary dummy entries.
		FileEntry fe = mwp.getMSTFileEntry("TBD1", jiwp.getDescription());
		sb.append("\nMST File Summary:\n");
		sb.append(fe.getSummary("\t"));
		// Append information about cluster file by creating a 
		// temporary dummy entry.
		fe = cwp.getClusterFileEntry("TBD2");
		sb.append("\nCluster File Summary:\n");
		sb.append(fe.getSummary("\t"));
		// Return the summary as a string.
		return sb.toString();
	}
	
	/**
	 * Helper method to obtain the MST file specified by the user.
	 * 
	 * This is a helper method that is used by the ClusterWizardPage
	 * to obtain the MST file specified by the user. The cluster
	 * wizard page uses this information to automatically provide a 
	 * default cluster file name to the user.
	 * 
	 * @return The MSTFile specified by the user.
	 */
	protected String getMSTFileName() {
		return mwp.getMSTFile();
	}
	
	/**
	 * Obtain the currently selected analyzer.
	 * 
	 * This method is a convenience method that can be used by other
	 * wizard pages to determine the type of the wizard that the user
	 * has currently selected. Currently this method is used by the
	 * HeuristicsWizardPage and ClusterWizardPage to enable or disable 
	 * configuration.
	 * 
	 * @return A predefined enumeration defining the type of analyzer that
	 * the user has currently chosen.
	 */
	protected FWAnalyzer.FWAnalyzerType getAnalyzerType() {
		return awp.getAnalyzerType();
	}
	
	/**
	 * Helper method to obtain the currently selected data set.
	 * 
	 * This is a helper method that is used by the MSTWizardPage
	 * to determine the EST file specified by the user. The MST
	 * wizard page uses this information to automatically provide a 
	 * default MST file name to the user.
	 * 
	 * @return The full absolute path to the ESTFile specified by 
	 * the user.
	 */
	protected DataSet getDataSet() {
		return jiwp.getDataSet();
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
	 * The job information wizard page that contains the description
	 * of the job.
	 */
	private final JobInfoWizardPage jiwp;
	
	/**
	 * The wizard page contains the tree selection to be used for clustering.
	 */
	private final TreeWizardPage twp;

	/**
	 * The analyzer wizard page contains the analyzer information to
	 * be used for constructing the MST for this job.
	 */
	private final AnalyzerWizardPage awp;
	
	/**
	 * The analyzer wizard page contains the heuristic information to
	 * be used for constructing the MST for this job.
	 */
	private final HeuristicsWizardPage hwp;
	
	/**
	 * The wizard page contains the prime heuristic information to
	 * be used for constructing the MST for this job.
	 */
	private final PrimesHeuristicWizardPage phwp;
	
	/**
	 * The filter wizard page that contains information about the various
	 * filters configured by the user to weed out problematic ESTs that
	 * could negatively impact the quality of clustering. 
	 */
	private final FiltersWizardPage fwp;
	
	/**
	 * The wizard page that collects information about the MST file
	 * and the server on which the job is to be run.
	 */
	private MSTWizardPage mwp;
	
	/**
	 * The wizard page that collects information about the cluster
	 * file in which clustering information is to be stored.
	 */
	private ClusterWizardPage cwp;

	/**
	 * The job entry that is finally created by this job wizard
	 * when the final wizard page is displayed. This entry is
	 * used by the final wizard page to perform its operations. 
	 */
	private ClusteringJob job;
	
	/**
	 * The main frame that logically owns this job. This value is
	 * used to create a job monitor if a job is successfully 
	 * submitted to run on a server.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user.
	 */
	private static final String OVERVIEW_MSG = "<html>" +
		"This wizard guides you through the process of creating a job.<br>" +
		"A job runs on a selected server (in parallel is so configured)<br>"+
		"and computes Minimum Spanning Tree (MST) data as well as<br>" +
		"clustering information for a given EST file. The EST file is<br>" +
		"chosen from a data set added to the workspace.<br><br>" + 
		"The time taken to compute MST data is O(<i>n</i><sup>2</sup>), where <i>n</i><br>" +
		"is the number of EST sequences in the EST file. It is best to<br>" +
		"run even medium (more than 25,0000) sets of ESTs in parallel<br>" +
		"to obtain the results quickly.<br><br>" +		
		"<b>Note</b>: Existing MST and cluster files can be added to the data<br>" +
		"set and do not have to be generated.<br></html>";

	/**
	 * The generated serialization UID (need to keep the compiler happy) 
	 */
	private static final long serialVersionUID = 804993573751301886L;
}
