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

package org.peace_tools.core;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.AbstractButton;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.ProgressMonitor;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TreeSelectionListener;

import org.peace_tools.core.job.clustering.ClusteringJobWizard;
import org.peace_tools.core.job.east.EASTJobWizard;
import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.data.JobListTableModel;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * The file menu helper for the application.
 * 
 * This class encapsulates the code related to the operations performed by
 * various menu items in the "Job" menu. This class is typically created once
 * from the MainFrame.createMenus() method. The primary motivation  for
 * introducing sub-menu handler classes is to improve code organization and
 * minimize code clutter. Note that the JobMenuHelper is essentially an event
 * handler that is set on the various menu items in this class. This helper
 * class also provides a  createJobMenu() method that actually creates
 * the file menu.
 */
public class JobMenuHelper extends AbstractMenuHelper
	implements ActionListener, ListSelectionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "File" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public JobMenuHelper(MainFrame mainFrame) {
		super(AbstractMenuHelper.HelperType.JOB_MENU, mainFrame);
	}
	
	/**
	 * Helper method to create the job menu.
	 * 
	 * This method performs the actual task of creating the "Job" menu. This
	 * method has been introduced to organize all the methods related to the
	 * "Job" menu into a single class. This method is invoked from the MainMenu
	 * class to create the "Job" menu.
	 * 
	 * @param toolbar The tool bar to which frequently used shortcuts can
	 * be added typically in the form of icons. If the tool bar is null, then
	 * shortcuts are not added.
	 * 
	 * @param vmh The view menu helper associated with the main frame. This
	 * entry is used to create menu entry for viewing the "Jobs" tab.
	 */
	public JMenu createJobMenu(JToolBar toolbar, AbstractMenuHelper vmh) {
		// Create the actual menu.
		JMenu jobMenu = new JMenu("Job  ");
		// Add various menu options to the main menu.
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.COMPUTE_MST, true));
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.COMPUTE_CLUSTERS, true));
		jobMenu.addSeparator();
		// Add menus for EAST related jobs
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.EAST_ASSEMBLY, true));
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.CLUSTER_AND_EAST_ASSEMBLY, true));
		jobMenu.addSeparator();
		
		// Context sensitive job operation menu items.
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.START_JOB_MONITOR, true));
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.STOP_JOB_MONITOR, true));
		jobMenu.addSeparator();
		
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.SHOW_JOB_DETAILS, true));
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.ABORT_JOB, true));
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.REMOVE_JOB, true));
		jobMenu.addSeparator();
		
		// Finally create the server entries in the menu.
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.SHOW_JOBS_ON_SERVER, true));
		jobMenu.add(getMenuItem(AbstractMenuHelper.ActionType.SHOW_MY_JOBS_ON_SERVER, true));
		jobMenu.addSeparator();
		jobMenu.add(vmh.getMenuItem(AbstractMenuHelper.ActionType.VIEW_JOBS, true));
		
		if (toolbar != null) {
			// Add tools to the toolbar.
			toolbar.add(getTool(AbstractMenuHelper.ActionType.COMPUTE_MST, true));
			toolbar.add(getTool(AbstractMenuHelper.ActionType.COMPUTE_CLUSTERS, true));
			// EAST assembly tools
			toolbar.add(getTool(AbstractMenuHelper.ActionType.EAST_ASSEMBLY, true));
			toolbar.add(getTool(AbstractMenuHelper.ActionType.CLUSTER_AND_EAST_ASSEMBLY, true));
		}
		// Return the job menu to the caller.
		return jobMenu;
	}

	@Override
	public void actionPerformed(ActionEvent event) {
		// Handle various action commands from popup menus or from the
		// tool bar associated with the jobs tab.
		final String cmd = event.getActionCommand();
		if ("NewMST".equals(cmd) || "NewCluster".equals(cmd)) {
			 // Display the wizard to create a new data set
			 ClusteringJobWizard jobWizard = 
				 new ClusteringJobWizard("Create New Clustering", mainFrame);
			 Utilities.centerPanel(mainFrame, jobWizard);
			 jobWizard.showWizard("http://www.peace-tools.org/downloads/manual.pdf#page=30");
			 return;
		}
		if ("east".equals(cmd) || "peace_east".equals(cmd)) {
			 EASTJobWizard ejw = EASTJobWizard.create("Assemble Clustered Data Set via EAST", mainFrame, 
					 "peace_east".equals(cmd));
			 if (ejw != null) {
				 Utilities.centerPanel(mainFrame, ejw);
				 ejw.showWizard("http://www.peace-tools.org/downloads/manual.pdf#page=30");
			 }
			 return;
		}
		// Rest of the menu options are context sensitive.
		if (job == null) {
			// Huh! no job is really selected. can't perform action
			return;
		}
		if ("startMonitor".equals(cmd)) {
			boolean result = JobMonitor.create(job, mainFrame);
			String msg = "Job monitoring thread was " + 
						 (result ? "" : "not ") + "started.";
			JOptionPane.showMessageDialog(mainFrame, msg, "Job Monitor", 
					result ? JOptionPane.INFORMATION_MESSAGE : JOptionPane.WARNING_MESSAGE);
		} else if ("stopMonitor".equals(cmd)) {
			// Try and interrupt the job monitor thread. But no guarantee
			// that the thread is in a position to actually stop.
			JobMonitor.interrupt(job);
		} else if ("deleteJob".equals(cmd)) {
			// Use the delete dialog to delete the job entry and 
			// remove associated files.
			DeleteDialog delDiag = new DeleteDialog(mainFrame, job);
			delDiag.setVisible(true);
		} else if ("viewJob".equals(cmd)) {
			Server server = Workspace.get().getServerList().getServer(job.getServerID());
			mainFrame.getViewFactory().createView(job, server);
		} else if ("abortJob".equals(cmd)) {
			// Use helper method to keep this method streamlined
			abortJob(job);
		} else if ("showJobsOnServer".equals(cmd)) {
			Server server = Workspace.get().getServerList().getServer(job.getServerID());
			mainFrame.getViewFactory().createView(server, null);
		} else if ("showMyJobsOnServer".equals(cmd)) {
			Server server = Workspace.get().getServerList().getServer(job.getServerID());
			mainFrame.getViewFactory().createView(server, server.getUserID());
		}
	}

	/**
	 * Method to abort a job running on a server.
	 * 
	 * This method is a convenience method that can be used to abort
	 * a job that is currently running on a server. This method 
	 * performs the following tasks:
	 * 
	 * 
	 * @param job The job to be aborted. This entry cannot be null.
	 */
	public void abortJob(Job job) {
		// First determine the server on which the job is running.
		Server server = Workspace.get().getServerList().getServer(job.getServerID());
		if (server == null) {
			// Can't find server entry. Check to see if user want's to
			// just tag job as being aborted.
			final String msg = String.format(CANT_FIND_SERVER, 
					job.getServerID(), job.getJobID());
			int choice = 
				JOptionPane.showConfirmDialog(mainFrame, msg,
					"Cannot abort job", JOptionPane.YES_NO_OPTION,
					JOptionPane.ERROR_MESSAGE);
			if (choice == JOptionPane.YES_OPTION) {
				// Simply tag the job as having finished unsuccessfully.
				job.setStatus(JobBase.JobStatusType.FAILED);
			}
			// We did all we could. 
			return;
		}
		// Double check to ensure that the user want's to really abort
		final String msg = String.format(SURE_ABORT_JOB, job.getJobID(),
				server.getName());
		// Check with the user
		int choice = 
			JOptionPane.showConfirmDialog(mainFrame, msg,
				"Abort job?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
		if (choice == JOptionPane.NO_OPTION) {
			// User changed his/her mind
			return;
		}
		// Get the job runner script to kill the job. 
		ProgressMonitor monitor = null;
		try {
			// Use a progress monitor to give user feedback if it takes a while.
			String waitMsg = "Please wait while job (Job ID: " + job.getJobID() +
				") is aborted...";
			monitor = new ProgressMonitor(mainFrame, waitMsg, "", 1, 5);
			monitor.setProgress(1);
			// First try and get a connection to the remote server.
			ServerSession session = SessionFactory.createSession(mainFrame, server);
			session.setPurpose(String.format(SESSION_PURPOSE, 
					job.getJobID(), server.getName()));
			session.connect();
			monitor.setProgress(2);
			// Now request the job runner script to abort the job.
			Server.OSType os       = session.getOSType();
			final String extension = (Server.OSType.WINDOWS.equals(os)) ? "bat" : "sh"; 
			String remoteCmd = job.getPath() + "/jobRunner." + extension
				+  " abort";
			String streams[] = {"", ""};
			ProgrammerLog.log("Attempting to run command '" + remoteCmd + 
					"' on server: " + server.getName() + "\n");
			monitor.setProgress(3);
			int exitCode = session.exec(remoteCmd, streams);
			monitor.setProgress(4);
			ProgrammerLog.log("Exit code: " + exitCode + "\n"); 
			if (exitCode != 0) {
				throw new IOException("The job runner script terminated" +
						"with exit code: " + exitCode + ".\n" +
						"The error message reported was:\n" + streams[1] + "\n.");
			}
			monitor.setProgress(5);
			// The job was successfully aborted.
			String successMsg = String.format(SUCCESS_MSG, 
					job.getJobID(), server.getName());
			JOptionPane.showMessageDialog(mainFrame, successMsg,
					"Job Aborted", JOptionPane.INFORMATION_MESSAGE);
		} catch (Exception e) {
			monitor.setProgress(5);
			final String errMsg = String.format(CANT_ABORT_JOB, 
					job.getJobID(), server.getName());
			JPanel fullInfo = Utilities.collapsedMessage(errMsg, Utilities.toString(e));
			JOptionPane.showMessageDialog(mainFrame, fullInfo,
					"Unable to Abort job", JOptionPane.ERROR_MESSAGE);
		}
	}
	
	@Override
	public ActionListener getActionListener() {
		return this;
	}

	@Override
	public ListSelectionListener getListSelectionListener(JTable table) {
		this.table = table;
		return this;
	}

    /**
     * The selection listener/handler for a table.
     * 
     * This method is invoked by the core Swing classes whenever the
     * user selects a specific entry in a list or table. This
     * method essentially enables/disables various tool bar buttons
     * and menu items based on the current job selection.
     * 
     * <p><b>Note:</b>  Currently we only handle JTable and not JList.</p>
     * 
     * @param event The selection event associated with this method.
     * This event is not really used for any major information. So
     * it can be null.
     */
	@Override
	public void valueChanged(ListSelectionEvent event) {
		if ((event != null) && (event.getValueIsAdjusting())) {
			return;
		}
		assert ( table != null );
		if (!(table.getModel() instanceof JobListTableModel)) {
			// The table model is not really compatible.
			return;
		}
		JobListTableModel model = (JobListTableModel) table.getModel();
		final int row           = table.getSelectedRow();
		// A couple of flags to help make checks below easier
		job = model.getJob(row);

		final boolean haveMonitor = ((job != null) && (job.getMonitor() != null));
        final boolean haveServer  = !"<n/a>".equals(model.getValueAt(row, 3));
		// Enable / disable tools and menus based on the job status.
		setEnabled("startMonitor", (job != null) && (!haveMonitor) && !job.isDone() && !job.isWaiting());
		setEnabled("stopMonitor",  (job != null) && (haveMonitor));
		setEnabled("viewJob",      (job != null));
		setEnabled("abortJob",     (job != null) && !job.isDone() && !job.isWaiting());
		setEnabled("deleteJob",    (job != null) && job.isDone());
		
		// If we have a server enabled then enable job view options
		setEnabled("showJobsOnServer",   haveServer);
		setEnabled("showMyJobsOnServer", haveServer);
	}

	@Override
	public JMenuItem getMenuItem(ActionType actionType, boolean mainMenu) {
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainMenu ? "24x24/" : "16x16/");
		// Create the actual menu item if we support this action. Note
		// that we immediately return menu options that are not context
		// sensitive.
		JMenuItem item = null;
		if (ActionType.COMPUTE_MST.equals(actionType)) {
			return Utilities.createMenuItem(Utilities.MENU_ITEM, "Create Job to Compute MST",
					(mainMenu ? MenuSubTitles[0] : null), 
					"NewMST", this, IconPath + "NewMST.png", 
					null, false, false);
		} else if (ActionType.COMPUTE_CLUSTERS.equals(actionType)) {
			return Utilities.createMenuItem(Utilities.MENU_ITEM, "Job to Compute Clusters",
					(mainMenu ? MenuSubTitles[1] : null),
					"NewCluster", this, IconPath + "NewCluster.png", 
					null, true, false);
		} else if (ActionType.EAST_ASSEMBLY.equals(actionType)) {
			return Utilities.createMenuItem(Utilities.MENU_ITEM, "Job to Assemble via EAST",
					(mainMenu ? MenuSubTitles[2] : null),
					"east", this, IconPath + "EAST.png", 
					null, true, false);
		} else if (ActionType.CLUSTER_AND_EAST_ASSEMBLY.equals(actionType)) {
			return Utilities.createMenuItem(Utilities.MENU_ITEM, "Job to Cluster & Assemble via EAST",
					(mainMenu ? MenuSubTitles[3] : null),
					"peace_east", this, IconPath + "PEACE_EAST.png", 
					null, true, false);
		} else if (ActionType.START_JOB_MONITOR.equals(actionType)) {
			item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Start Job Monitor",
					(mainMenu ? MenuSubTitles[4] : null),
					"startMonitor", this, IconPath + "StartJobMonitor.png", 
					null, false, false);
		} else if (ActionType.STOP_JOB_MONITOR.equals(actionType)) {
			item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Stop Job Monitor",
					(mainMenu ? MenuSubTitles[5] : null),
					"stopMonitor", this, IconPath + "StopJobMonitor.png", 
					null, false, false);
		} else if(ActionType.SHOW_JOB_DETAILS.equals(actionType)) {
			item = Utilities.createMenuItem(Utilities.MENU_ITEM, "View Job Details",
					(mainMenu ? MenuSubTitles[6] : null),
					"viewJob", this, IconPath + "ViewJob.png", 
					null, false, false);
		} else if (ActionType.ABORT_JOB.equals(actionType)) {
			item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Abort Job",
					(mainMenu ? MenuSubTitles[7] : null),
					"abortJob", this, IconPath + "AbortJob.png", 
					null, false, false);
		} else if (ActionType.REMOVE_JOB.equals(actionType)) {
			item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Remove Job",
					(mainMenu ? MenuSubTitles[8] : null),
					"deleteJob", this, IconPath + "Delete.png", 
					null, false, false);
		} else if (ActionType.SHOW_JOBS_ON_SERVER.equals(actionType)) {
			item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Show all Jobs on Server",
					(mainMenu ? MenuSubTitles[9] : null),
					"showJobsOnServer", this, IconPath + "ServerInfo.png", 
					null, false, false);
		} else if (ActionType.SHOW_MY_JOBS_ON_SERVER.equals(actionType)) {
			item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Show my jobs on Server",
					(mainMenu ? MenuSubTitles[10] : null),
					"showMyJobsOnServer", this, IconPath + "ServerMyJobs.png", 
					null, false, false);
		}
		// When control drops here that means the menu item is
		// context sensitive. So track these menu items.
		if (item != null) {
			contextItemList.add(item);
			Utilities.setEnabled(item, false);
		}
		return item;
	}

	@Override
	public AbstractButton getTool(ActionType actionType, boolean mainToolBar) {
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainToolBar ? "24x24/" : "16x16/");
		// Create the actual menu item if we support this action. Note
		// that we immediately return menu options that are not context
		// sensitive.
		AbstractButton tool = null;
		if (ActionType.COMPUTE_MST.equals(actionType)) {
			return Utilities.createToolButton(IconPath + "NewMST.png", null,
 					"NewMST", this, MenuSubTitles[0], false);
		} else if (ActionType.COMPUTE_CLUSTERS.equals(actionType)) {
			return Utilities.createToolButton(IconPath + "NewCluster.png", null,
					"NewCluster", this, MenuSubTitles[1], true);
		} else if (ActionType.EAST_ASSEMBLY.equals(actionType)) {
			return Utilities.createToolButton(IconPath + "EAST.png", null,
					"east", this, MenuSubTitles[2], true);
		} else if (ActionType.CLUSTER_AND_EAST_ASSEMBLY.equals(actionType)) {
			return Utilities.createToolButton(IconPath + "PEACE_EAST.png", null,
					"peace_east", this, MenuSubTitles[3], true);
		} else if (ActionType.START_JOB_MONITOR.equals(actionType)) {
			tool = Utilities.createToolButton(IconPath + "StartJobMonitor.png", null,					
					"startMonitor", this, MenuSubTitles[4], false);
		} else if (ActionType.STOP_JOB_MONITOR.equals(actionType)) {
			tool = Utilities.createToolButton(IconPath + "StopJobMonitor.png", null, 
					"stopMonitor", this, MenuSubTitles[5], false);
		} else if(ActionType.SHOW_JOB_DETAILS.equals(actionType)) {
			tool = Utilities.createToolButton(IconPath + "ViewJob.png", null, 
					"viewJob", this, MenuSubTitles[6], false);
		} else if (ActionType.ABORT_JOB.equals(actionType)) {
			tool = Utilities.createToolButton(IconPath + "AbortJob.png", null,
					"abortJob", this,  MenuSubTitles[7], false);
		} else if (ActionType.REMOVE_JOB.equals(actionType)) {
			tool = Utilities.createToolButton(IconPath + "Delete.png", null,
					"deleteJob", this, MenuSubTitles[8], false);
		} else if (ActionType.SHOW_JOBS_ON_SERVER.equals(actionType)) {
			tool = Utilities.createToolButton(IconPath + "ServerInfo.png", null,
					"showJobsOnServer", this, MenuSubTitles[9], false);
		} else if (ActionType.SHOW_MY_JOBS_ON_SERVER.equals(actionType)) {
			tool = Utilities.createToolButton(IconPath + "ShowMyJobs.png", null,
					"showMyJobsOnServer", this, MenuSubTitles[10], false);
		}
		// When control drops here that means the menu item is
		// context sensitive. So track these menu items.
		if (tool != null) {
			contextItemList.add(tool);
			tool.setEnabled(false);
		}
		
		return tool;
	}

	@Override
	public TreeSelectionListener getTreeSelectionListener(JTree tree) {
		return null;
	}
	
	/**
	 * The currently selected job entry (if any). This value is
	 * updated by the list selection listener whenever a valid
	 * job entry is selected. If a valid job entry is not selected
	 * then this entry is set to null.
	 */
	private Job job;
	
	/**
	 * The various sub menu titles that are used in the main menu. The
	 * sub menu titles are used to provide the user with a bit more 
	 * verbose description on the action that will be performed by a given
	 * menu item.
	 */
	private static final String MenuSubTitles[] = {
		"Starts wizard to schedule job to compute MST and clustering for a ESTs in a data set",
		"Computes new clustering using existing MST data (quick operation)",
		"Assemble via EAST using existing clustering solution in workspace",
		"Cluster & Assemble via EAST for cDNA fragments in a data set",
		"Starts job monitor (daemon thread) for selected job",
		"Stop job monitor (daemon thread) for selected job",
		"Display ouputs, errors, and scripts for selected job",
		"Try and abort the selected job",
		"Remove currently selected job from work space",
		"Show all jobs running on the same server",
		"Show only my jobs running on the same server"
	};
	
	/**
	 * An error message that is formatted (to fill-in missing
	 * information) and displayed to the user when a server
	 * entry could not be successfully located to abort a job.
	 */
	private static final String CANT_FIND_SERVER = "<html>" +
			"Unable to find entry for server ID %s<br>" +
			"in workspace. Without server entry the job ID: %s<br>" +
			"cannot be really aborted. If a server entry was not<br>" +
			"found, there is not a whole lot you can do with this job." +
			"<br><br>Do you want to simply flag it as aborted?<br>" +
			"</html>";
	
	/**
	 * An verification message that is formatted (to fill-in missing
	 * information) and displayed to the user just before a
	 * job is aborted on a server.
	 */
	private static final String SURE_ABORT_JOB = "<html>" +
			"Are you sure you want to abort the job with Job ID %s<br>" +
			"that is running on server: %s?<br><br>" +
			"<i>Note that if the job is aborted work done by the job<br>" +
			"will be lost and no data files will be generated." +
			"</html>";
	
	/**
	 * An error message that is formatted (to fill-in missing
	 * information) and displayed to the user if aborting a job
	 * fails.
	 */
	private static final String CANT_ABORT_JOB = "<html>" +
			"Unable to abort the job with Job ID %s<br>" +
			"that is running on server: %s.<br>" +
			"because an error occured. Refer to the details below<br>" +
			"for additional information about the error." +
			"</html>";
	
	/**
	 * An error message that is formatted (to fill-in missing
	 * information) and displayed to the user as the reason
	 * for making a remote server connection.
	 */
	private static final String SESSION_PURPOSE = "<html>" +
			"Attempting to abort the Job with ID %s<br>" +
			"that is running on server: %s." +
			"</html>";

	/**
	 * A success message that is formatted (to fill-in missing
	 * information) and displayed to the user indicating that
	 * the job was aborted.
	 */
	private static final String SUCCESS_MSG = "<html>" +
		"The job with job ID %s that was running on<br>" +
		"server %s has been aborted." +
		"</html>";
}
