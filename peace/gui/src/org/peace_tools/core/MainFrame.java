//--------------------------------------------------------------------

//This file is part of PEACE.

//PEACE is free software: you can redistribute it and/or modify it
//under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//PEACE is distributed in the hope that it will be useful, but
//WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with PEACE.  If not, see <http://www.gnu.org/licenses/>.

//Miami University makes no representations or warranties about the
//suitability of the software, either express or implied, including
//but not limited to the implied warranties of merchantability,
//fitness for a particular purpose, or non-infringement.  Miami
//University shall not be liable for any damages suffered by licensee
//as a result of using, result of using, modifying or distributing
//this software or its derivatives.

//By using or copying this Software, Licensee agrees to abide by the
//intellectual property laws, and all other applicable laws of the
//U.S., and the terms of GNU General Public License (version 3).

//Authors:   Dhananjai M. Rao          raodm@muohio.edu

//---------------------------------------------------------------------

package org.peace_tools.core;

import java.awt.BorderLayout;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URI;
import java.util.HashMap;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JToolBar;
import javax.swing.SwingUtilities;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.AbstractMenuHelper.HelperType;
import org.peace_tools.generic.CustomBorder;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.dndTabs.DnDTabbedPane;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

public class MainFrame extends JFrame implements ActionListener {
	/**
	 * The permanent tab that is always present in in the main frame.
	 * This pane is created in the constructor and is never changed. 
	 */
	private final DnDTabbedPane centerPane;

	/**
	 * This is the parent component of all the tabs and tabbed panes
	 * that constitute the visual layout of the main frame. This
	 * component is always present and is never removed from the
	 * main frame.
	 */
	private final JPanel desktop;

	/**
	 * Sets up the default layout for the main frame.
	 * The constructor uses the data in the current workspace to setup
	 * the various panels and components. 
	 * 
	 *  @note Prior to creating the main frame ensure that a valid
	 *  operational work space has been created.
	 *  
	 *  @param firstLaunch If this flag is true, that indicates that
	 *  PEACE GUI is being launched for the first time. In this case,
	 *  display the "Welcome" view to the user.
	 */
	public MainFrame() {
		super();
		assert ( Workspace.get() != null );
		setTitle("PEACE: " + Workspace.get().getDirectory());
		setIconImage(Utilities.getIcon("images/16x16/PEACE.png").getImage());
		// Setup the layout.
		setLayout(new BorderLayout(0, 0));
		setPreferredSize(new Dimension(1024, 768));
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		// Create menus and tool bars.
		createMenus();
		// First create our main center pane.
		centerPane = new DnDTabbedPane();
		centerPane.setPermanent(true);
		centerPane.setPreferredSize(getPreferredSize());
		// Create our desktop and add it to the main frame.
		desktop = new JPanel(new BorderLayout(0, 0));
		desktop.setBorder(new EmptyBorder(2, 2, 2, 2));
		desktop.add(centerPane);
		// Add the permanent tab to our frame. This must be done
		// before splitting the centerPane below.
		add(desktop, BorderLayout.CENTER);
		// Create the standard views via the view factory.
		defaultViewFactory = new DefaultViewFactory(this);
		defaultViewFactory.createStaticViews();
	}

	/**
	 * This is a helper method that invoked from the constructor to
	 * create the main menu and associated toolbars. This method was
	 * introduced to keep the code clutter down in the constructor.
	 */
	private void createMenus() {
		// First create the various menu action listeners so that
		// they can be cross registered as needed.
		FileMenuHelper   fm  = new FileMenuHelper(this);
		JobMenuHelper    jmh = new JobMenuHelper(this);
		ServerMenuHelper smh = new ServerMenuHelper(this);
		ViewMenuHelper   vm  = new ViewMenuHelper(this);
		HelpMenuHelper   hm  = new HelpMenuHelper(this);
		// Save references for quick look up.
		menuHelpers.put(HelperType.FILE_MENU, fm);
		menuHelpers.put(HelperType.JOB_MENU, jmh);
		menuHelpers.put(HelperType.SERVER_MENU, smh);
		menuHelpers.put(HelperType.VIEW_MENU, vm);
		menuHelpers.put(HelperType.HELP_MENU, hm);
		
		// Next create tool bar and the main menu using various
		// helpers.
		toolbar  = new JToolBar();
		toolbar.setFloatable(false);
		toolbar.setBorder(new CompoundBorder(new CustomBorder("ssds"),
				new EmptyBorder(2, 6, 2, 6)));
		// Create all the menu items.
		JMenuBar mainmenu = new JMenuBar();
		mainmenu.add(fm.createFileMenu(toolbar, jmh));
		// Create the Server menu.
		mainmenu.add(smh.createServerMenu(toolbar));
		// The job menu goes next
		mainmenu.add(jmh.createJobMenu(toolbar, vm));
		// Create view menu next
		mainmenu.add(vm.createViewMenu(toolbar, hm));
		// Create the help menu and set it up
		mainmenu.add(hm.createHelpMenu(toolbar));
		// Setup the main menu and tool bar appropriately
		setJMenuBar(mainmenu);
		add(toolbar, BorderLayout.NORTH);
	}

	/**
	 * Obtain the permanent desktop panel for this frame.
	 * 
	 * @return The permanent desktop panel for this frame.
	 */
	JPanel getDesktop() { return this.desktop; }

	/**
	 * Obtain the permanent, center pane for this frame.
	 * 
	 * @return The permanent center pane for this frame.
	 */
	public DnDTabbedPane getCenterPane() { return this.centerPane; }

	/**
	 * Obtain the toolbar associated with this frame.
	 * 
	 * @return The toolbar associated with this frame. The return value
	 * is never null. However, the toolbar may not be visible at all
	 * times. 
	 */
	public JToolBar getToolBar() { return toolbar; }

	/**
	 * The action performed method that handles various application-level
	 * events. One event that this method processes is notifications from
	 * JobMonitor threads when they exit. 
	 * 
	 * @see JobMonitor
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		if ("JobMonitor".equals(e.getActionCommand())) {
			Job job = (Job) e.getSource();
			// Get server information about the job for messages below.
			Server srvr = Workspace.get().getServerList().getServer(job.getServerID());
			String srvrName = (srvr != null) ? srvr.getName() : job.getServerID();
			// Perform various tasks based on job exit status.
			if (JobBase.JobStatusType.FINISHING.equals(job.getStatus())) {
				// The job has finished successfully. Pull down the data
				// files from the remote to the local machine for viewing.
				String msg = String.format(JOB_DONE_MSG, job.getJobID(), srvrName);
				int decision = JOptionPane.showConfirmDialog(this, msg,
						"Job completed", JOptionPane.YES_NO_OPTION);
				if (decision == JOptionPane.YES_OPTION) {
					// Copy the files generated for this job now.
					FileCopyDialog fcd = new FileCopyDialog(this, job);
					fcd.showDialog();
				}
			} else if (JobBase.JobStatusType.FAILED.equals(job.getStatus())) {
				// The job did not complete successfully. Let the user know
				// about it either ways.
				String msg = String.format(JOB_ERR_MSG, job.getJobID(), srvrName);
				JOptionPane.showMessageDialog(this, msg, "Job error report", 
						JOptionPane.ERROR_MESSAGE);
			} else {
				// Hmmm...This is a weird situation. log it
				UserLog.log(UserLog.LogLevel.WARNING, "JobMonitor",
						"The job monitor for job " + job.getJobID() + 
						" running on server " + srvrName + 
						" has unexpectedly terminated. You may want to " +
				"restart the job monitor for this job.");
			}
		}
	}

	/**
	 * Save the workspace data directly.
	 * 
	 * This method saves the work space data via the work space
	 * file.
	 * 
	 * @param dialog If the dialog is not null, then this method
	 * hides and disposes the dialog. This parameter can be null.
	 */
	public void saveWorkspace(final Dialog dialog) {
		try {
			Thread.sleep(100); // Maybe this sleep is not needed. but helps visually
			Workspace.get().saveWorkspace();
		} catch (Exception e) {
			ProgrammerLog.log(e);
			UserLog.log(UserLog.LogLevel.ERROR, 
					"PEACE", e.getMessage());
			JPanel msg = Utilities.collapsedMessage(SAVING_ERROR,
					Utilities.toString(e));
			JOptionPane.showMessageDialog(MainFrame.this, msg,
					"Error saving workspace", JOptionPane.ERROR_MESSAGE);
		} finally {
			// If dialog is not null, the arrange to have the dialog closed
			// from the swing thread. Otherwise closing the dialog may hang
			// in some cases.
			if (dialog != null) {
				SwingUtilities.invokeLater(new Runnable() {
					@Override
					public void run() {
						dialog.setVisible(false);
						dialog.dispose();
					}
				});
			}
			setCursor(Cursor.getDefaultCursor());
		}
	}

	/**
	 * This method is a helper method to save the workspace from a 
	 * separate thread. This method creates a dialog indicating that
	 * the work space is being saved and then saves the work space
	 * from another thread so as not to hold up the GUI thread.
	 */
	public synchronized void saveDelayedWorkspace() {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				// First create a dialog to display the message.
				final JDialog dialog = new JDialog((Frame) null, 
						"Saving workspace", true);
				JLabel info = new JLabel("Saving workspace. Please wait...",
						Utilities.getIcon("images/16x16/HourGlass.png"), JLabel.LEFT);
				info.setBorder(new EmptyBorder(20, 20, 20, 20));
				dialog.setLayout(new BorderLayout(0, 0));
				dialog.add(info, BorderLayout.CENTER);
				dialog.setIconImage(Utilities.getIcon("images/16x16/PEACE.png").getImage());
				dialog.pack();
				Utilities.centerPanel(MainFrame.this, dialog);
				dialog.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
				dialog.setAlwaysOnTop(true);
				dialog.setResizable(false);
				dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
				// Create a thread to do the actual saving so as not
				// to block the GUI thread.
				Thread saver = new Thread(new Runnable() {
					@Override
					public void run() {
						saveWorkspace(dialog);
					}
				});
				// Setup a wait cursor
				setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
				// Start thread in the background.
				saver.start();
				// Finally make the dialog visible.
				dialog.setVisible(true);
			}
		});
	}

	/**
	 * Obtain the view factory associated with this main frame.
	 * 
	 * @return The view factory associated with this main frame.
	 */
	public ViewFactory getViewFactory() {
		return defaultViewFactory;
	}

	/** Displays a given URL.
	 * 
	 * This method is used when the user chooses to view one of
	 * the selected help topics.  This method is actually invoked 
	 * from various classes to display help. The help is actually
	 * redirected to the PEACE web site.  Having users visit a web
	 * site permits help content to be developed independently
	 * and updated consistently.
	 * 
	 * @param url The complete web site URL that must be displayed in 
	 * the default system browser.
	 */
	public void showHelp(String url) {
		if (Desktop.isDesktopSupported() && 
				(Desktop.getDesktop().isSupported(Desktop.Action.BROWSE))) {
			// Desktop is supported with browse to launch browser
			try {
				Desktop.getDesktop().browse(new URI(url));
			} catch (Exception e) {
				// Log the exception
				ProgrammerLog.log(e);
				// Display error to the user.
				UserLog.log(UserLog.LogLevel.ERROR, 
						"PEACE", e.getMessage());
				JPanel msg = Utilities.collapsedMessage(HELP_ERROR,
						Utilities.toString(e));
				JOptionPane.showMessageDialog(this, msg,
						"Error Displaying Help", JOptionPane.ERROR_MESSAGE);
			}
		} else {
			JOptionPane.showMessageDialog(this, 
					"Your Java version does not provide the necessary features\n" +
					"to launch the online help directly from PEACE GUI.\n" +
					"You may view the requested help information via the URL:\n" +
					url, "Java Feature Unavailable", JOptionPane.WARNING_MESSAGE);
		}
	}

	/**
	 * Method to create threads for all jobs whose status is to be monitored.
	 * 
	 * This method is typically invoked from PEACE.launchMainFrame() method
	 * after the main frame is created. This method performs the task of 
	 * checking job status and creating a job thread if a job needs to be
	 * monitored.
	 * 
	 * @note Invoking this method twice will cause unnecessary monitoring
	 * threads to start up. So avoid duplicate calls.
	 */
	public void createJobThreads() {
		JobList jobList = Workspace.get().getJobList();
		for(Job job: jobList.getJobs()) {
			if (!JobBase.JobStatusType.FAILED.equals(job.getStatus()) &&
				!JobBase.JobStatusType.SUCCESS.equals(job.getStatus())) {
				// This job requires some monitoring. Start thread for this job.
				JobMonitor.create(job, this);
			}
		}
	}

	/**
	 * The default view factory that is used to create all the
	 * views associated with this main frame.
	 */
	private DefaultViewFactory defaultViewFactory = null;

	/**
	 * A default implementation for a view factory. This approach
	 * is a feeble attempt at trying to control the instantiation
	 * of a view factory.
	 */
	private class DefaultViewFactory extends ViewFactory {
		/**
		 * The constructor. 
		 * 
		 * The constructor merely passes the reference to the main frame
		 * to the base class.
		 * 
		 * @param mf The main frame that logically owns this view factory.
		 */	
		public DefaultViewFactory(MainFrame mf) {
			super(mf);
		}
	}

	/**
	 * Obtain the thread group that contains the job status monitoring
	 * threads.
	 * 
	 * @return The thread group that logically contains all the job
	 * status monitoring threads that are currently active.
	 */
	protected static ThreadGroup getWorkerThreads() { return workerThreads; }


	/**
	 * Obtain a menu helper for a given helper type.
	 * 
	 * This method can be used to obtain a menu helper for a given
	 * type of action to be performed.
	 * 
	 * @param helperType The class (or category) of operation for 
	 * which the menu helper object is desired.
	 * 
	 * @return The menu helper object for the specified category. If
	 * the helper object is not found, then this method returns null.
	 */
	public AbstractMenuHelper getMenuHelper(HelperType helperType) {
		return menuHelpers.get(helperType);
	}

	/**
	 * The list of menu helpers associated with this main frame.
	 * This hash map is used to maintain and rapidly access various
	 * menu helpers associated with the main frame. The list is 
	 * populated in the createMenus() method.
	 */
	private final HashMap<AbstractMenuHelper.HelperType, AbstractMenuHelper>
		menuHelpers = new HashMap<AbstractMenuHelper.HelperType, AbstractMenuHelper>();
	
	/**
	 * A reference to the toolbar associated with this main frame. The
	 * tool bar may be visible or invisible at any given time.
	 */
	private JToolBar toolbar;

	/**
	 * A thread group that contains various worker threads. This
	 * thread group is created when a main frame is instantiated.
	 * Threads are created whenever the status of a job needs to
	 * be asynchronously monitored and updated by the JobWizard.
	 */
	private static transient final ThreadGroup workerThreads 
		= new ThreadGroup("JobStatusUpdateThreads");

	/**
	 * This message is displayed by the main frame when a job being
	 * monitored does not complete successfully. This message is
	 * formatted with suitable values before it is displayed to
	 * the user.
	 */
	private static final String JOB_ERR_MSG = "<html>" +
	"The job (ID: %1$s) that was scheduled to run on<br>" +
	"server: %2$s<br>" +
	"did not complete successfully. You can review the outputs<br>" +
	"from the job to determine the root cause of failure.</html>";

	/**
	 * This message is displayed by the main frame when a job being
	 * monitored does complete successfully and the files need to be
	 * copied over. This message is formatted with suitable values 
	 * before it is displayed to the user.
	 */
	private static final String JOB_DONE_MSG = "<html>" +
	"The job (ID: %1$s) that was scheduled to run on<br>" +
	"server: %2$s<br>" +			
	"has completed successfully. The generated data files<br>" +
	"need to copied to be viewed in the GUI. <br><br>" +
	"<b>Would you like to copy the generated file(s) now?</b></html>";

	/**
	 * This message is displayed when the work space file could not
	 * be saved successfully.
	 */
	private static final String SAVING_ERROR = "<html>" +
	"The workspace file could not be saved. This is a serious issue.<br>" +
	"Plese refer to the details to identify the root cause. You may try<br>" +
	"to save the workspace via the menus." +
	"</html>";

	/**
	 * This message is displayed when the default system browser could
	 * not be launched to display help.
	 */
	private static final String HELP_ERROR = "<html>" +
	"The default system browser could not be launched to display help.<br>" +
	"You may view the help web site directly from your system browser by<br>" +
	"navigating to the appropriate help section from http://www.peace-tools.org/" +
	"</html>";

	/**
	 * A generated serial version ID.
	 */
	private static final long serialVersionUID = 1414403952944170343L;
}
