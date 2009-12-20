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

import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.UIManager;

import org.peace_tools.generic.Log;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Workspace;

/**
 * The top-level class for the PEACE GUI application.
 * 
 * This class merely contains the main method that jump starts the various
 * operations including creation of the top-level frame for normal use.
 */
public class PEACE {
	/**
	 * The main method launches the core dialogs in the system to start
	 * up on the swing thread. Starting the GUI on the swing thread is
	 * important to ensure that any custom UI installed with PEACE
	 * starts up correctly.

	 * @param args The command line arguments are currently unused.
	 */
	public static void main(String[] args) {
		 // Schedule a job for the event-dispatching thread:
        // creating and showing this application's GUI.
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                new PEACE();
            }
        });
	}
	
	/**
	 * Call back method from the work space chooser. This method is
	 * called from the work space chooser after the user has selected
	 * a valid work space to work. This method launches the main frame
	 * that takes the place of the work space chooser.
	 * 
	 * @param workspacePath The full absolute path to the work space.
	 * 
	 * @param firstLaunch If this flag is true, then it is assumed that
	 * PEACE GUI is being run for the first time and a welcome screen 
	 * with some basic instructions is displayed to the user.
	 */
	protected void launchMainFrame(String workspacePath, final boolean firstLaunch) {
		if ((workspacePath != null) && (Workspace.get() != null)) {
			// Launch the main frame and we are done.
			final MainFrame mf = new MainFrame();
			mf.pack();
			mf.setVisible(true);
			// Cut logs in the user log to provide core information.
			UserLog.log(Log.LogLevel.NOTICE, "PEACE", Version.GUI_VERSION);
			UserLog.log(Log.LogLevel.NOTICE, "PEACE", Version.COPYRIGHT);
			UserLog.log(Log.LogLevel.NOTICE, "PEACE", Version.DISCLAIMER);
			// Cut information in the programmer log about environment
			ProgrammerLog.log(Version.GUI_VERSION + "\n");
			ProgrammerLog.log(Version.COPYRIGHT + "\n");
			// Create welcome tab and the job threads for monitoring
			// jobs after some delay
			// to let the GUI settle in first. So we do it from a 
			// separate thread.
			Thread bgThread = new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						// Give a couple of seconds for the GUI to settle
						Thread.sleep(2000);
						// Create threads for monitoring job status from
						// the main Swing thread.
						SwingUtilities.invokeLater(new Runnable() {
							@Override
							public void run() {
								// First create the welcome tab if needed
								// Show the welcome screen if it is not already there.
								if (firstLaunch) {
									ViewFactory vf = mf.getViewFactory();
									vf.createView("installFiles/welcome.html", null,
											ViewFactory.ViewType.HTML_VIEW, false, false);
								}
								// Create job threads for all jobs that need them.
								mf.createJobThreads();
							}
						});
					} catch (InterruptedException e) {
						ProgrammerLog.log(e);
					}
				}
			});
			bgThread.start();
		}
	}
	
	/**
	 * The constructor. This is the main PEACE application method 
	 * that launches the core dialogs in the system. Specifically it
	 * performs the following tasks:
	 * 
	 * <ol>
	 * 
	 * <li>It launches a the WorkspaceChooser dialog to permit the user to
	 * choose a work space directory. This dialog also handles the task of
	 * asking the user to accept the license during first launch of PEACE.</li>
	 * 
	 * <li>If the user chooses a valid work space this method launches the main
	 * frame that does rest of the normal operational tasks.</li>
	 * 
	 * </ol>
	 * 
	 * @param args
	 */
	private PEACE() {
		// Turn off metal's use of bold fonts
		UIManager.put("swing.boldMetal", Boolean.FALSE);
		// try {
		//	UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		// } catch (Exception e) {}
		
		// Make tool tips show up faster than normal
		ToolTipManager ttm = ToolTipManager.sharedInstance();
		ttm.setInitialDelay(3 * ttm.getInitialDelay() / 4);
		
		// Create and show the top-level work space chooser class.
		WorkspaceChooser wsChooser = new WorkspaceChooser(null, this);
		if (!wsChooser.createTabs()) {
			// Error occurred during initialization. Bail out.
			wsChooser.dispose();
			return;
		}
		// The chooser was created successfully.
		wsChooser.pack();
		// Center the dialog box on the screen.
		Utilities.centerPanel(null, wsChooser);
		wsChooser.setVisible(true);
		// Dispose the dialog box.
		wsChooser.dispose();
	}
}
