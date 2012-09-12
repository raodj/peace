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
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JLabel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.decagon.helpers.DADXHelper;
import org.peace_tools.decagon.jaxb.Executable;
import org.peace_tools.decagon.jaxb.OutputCheckData;
import org.peace_tools.generic.BackgroundTask;
import org.peace_tools.generic.BackgroundTask.UserTask;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProcessOutputDisplay.StyleKind;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;

/**
 * This class is used to verify that the executables used in a job
 * description are operational. This page uses the list of
 * executables specified in a DADX description and verifies that
 * each one of them are operational. 
 */
public class CheckExecutables extends GenericWizardPage implements UserTask {
	/**
	 * The constructor. 
	 * 
	 * The constructor sets up the various components
	 * on this wizard page. The constructor does not setup
	 * the GUI components. The GUI components are created
	 * just before this wizard page is to be displayed.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param ssp The server selection page from where the currently
	 * selected server entry is to be obtained to display as the
	 * chosen server to be displayed. 
	 */
	public CheckExecutables(DecagonJobCreator dcm, ServerSelectionPage ssp) {
		this.djc         = dcm;
		this.srvrSelPage = ssp;
		assert(this.djc != null);
		// Setup the title(s) for this page and border
		setTitle("Check Executables", 
				"Verify executables for job(s) are operational");
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
		// Create information about each executable as sub-steps.
		DADXHelper helper = djc.getDADXHelper();
		final List<Executable> exeList = helper.getAssemblerDetails().getExecutable();
		String[] subSteps = new String[(exeList == null) ? 0 : exeList.size()];
		for(int i = 0; (i < subSteps.length); i++) {
			final Executable exeEntry = exeList.get(i);
			final String     exeName  = (exeEntry.getFileName() != null ? exeEntry.getFileName() : 
				exeEntry.getInternal().toString());
			subSteps[i] = "<html>" + exeName +
				(exeEntry.getFileName() == null ? " (<i>internal</i>)" : "") +
				"<br/><font size=\"-2\">" + 
				Utilities.wrapStringToHTML(exeEntry.getDescription(), 60) +  
				"</font></html>";
		}
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
				taskRunner.start(false, CheckExecutables.this, "ExeCheck", null);
			}
		});
	}


	@Override
	public void run(BackgroundTask bTask) throws Exception {
		// Obtain a session (or connection) to the server. 
		final Server srvr = srvrSelPage.getSelectedServer();
		bTask.log("Creating session for " + srvr + "\n");
		final ServerSession session = 
			SessionFactory.createSession(this, srvrSelPage.getSelectedServer());
		// Establish a connection with the server. This call 
		// prompts the user for password if needed. If connection fails,
		// the following method throws an exception.
		bTask.log("Attempting to connect to " + srvr + "\n", StyleKind.INFO);
		session.connect();
		bTask.log("Successfully connected to " + srvr + "\n", StyleKind.INFO);
		// Obtain the list of executables to be verified.
		DADXHelper helper = djc.getDADXHelper();
		final List<Executable> exeList = helper.getAssemblerDetails().getExecutable();
		// Check each of the executables are operational.
		for(Executable exe: exeList) {
			// Obtain the command to be executed to verify operational
			// status of a given executable.
			bTask.log("Verifying executable with name: " + exe.getName() + "\n", StyleKind.BOLD_INFO);
			final String command = helper.getExeName(exe, srvr) + " " + 
				helper.getCmdLineArgs(exe.getDetectArgsList().getArgument());
			bTask.log("\tRunning command: " + command + "\n", StyleKind.INFO);
			// Attempt to run the executable on the server.
			final String outputs[] = new String[2];
			final int exitCode     = session.exec(command, outputs);
			// Log the information we got.
			bTask.log("\tExit code      : " + exitCode   + "\n", StyleKind.INFO);
			bTask.log("\tStandard output: " + outputs[0] + "\n", StyleKind.INFO);
			bTask.log("\tStandard error : " + outputs[1] + "\n", StyleKind.INFO);
			// Check if the exit code and outputs match the values specified
			// by the user.
			boolean success = (exe.getDetectExitCode().intValue() == exitCode);
			if (success && (exe.getDetectOutputCheck() != null)) {
				// In addition to exit code we need to check outputs
				final OutputCheckData ocd = exe.getDetectOutputCheck();
				Pattern regEx   = Pattern.compile(ocd.getValue().trim());
				Matcher checker = regEx.matcher(outputs[ocd.isStdout()? 0 : 1]);
				success         = checker.find();
			}
			// Check if the output matched.
			if (success) {
				bTask.log("Executable successfully checked.\n\n", StyleKind.BOLD_INFO);
			} else {
				bTask.log("Executable did not validate successfully\n", StyleKind.WARNING);
				bTask.log("Expected exit code: " + exe.getDetectExitCode() + "\n", StyleKind.WARNING);
				// Log any expected output values
				final OutputCheckData ocd = exe.getDetectOutputCheck();
				if (ocd != null) {
					bTask.log("Expected output on "  + (ocd.isStderr() ? "stderr" : "stdout") +
							" is: " + ocd.getValue(), StyleKind.WARNING);
				}
				throw new Exception("Executable validation failed. Please see logs above");
			}
			// One step done. Update progress.
			bTask.updateProgress();
		}
	}
	
	@Override
	public void done(BackgroundTask bTask) {
		// Enable various button once background operations are done.
		final boolean success = (bTask.getResult() == null);
		djc.setButtonStatus(1, (success ? 1 : 0), -1);
	}

	@Override
	public void dialogClosed(BackgroundTask bTask) {
		// Currently this method (that implements an interface from
		// BackgroundTask.UserTask interface) does not have any
		// operations to perform and is intentionally left blank.
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to obtain summary information
	 * from various pages constituting this wizard.
	 */
	private final DecagonJobCreator djc;

	/**
	 * This object is initialized in the constructor to point to the 
	 * page in the wizard in which the server to be used for running
	 * the pipeline was setup by the user. This object is never
	 * null.
	 */
	private final ServerSelectionPage srvrSelPage;
	
	/**
	 * This is a generic object that is used to run checks for
	 * executables in the background. This class provides features
	 * for displaying status information about the ongoing checks
	 * in the background thread.
	 */
	private BackgroundTask taskRunner;

	/**
	 * A simple static informational message that is displayed at the 
	 * top of this page.
	 */
	private static final String INFO_MSG = 
		"<html>Please wait while operational status of various<br/>" +
		"executables used in the job(s) are verified.</html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
