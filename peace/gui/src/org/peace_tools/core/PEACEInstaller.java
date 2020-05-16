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

import java.awt.BorderLayout;
import java.awt.Color;
import java.io.InputStream;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import javax.swing.border.EtchedBorder;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.generic.ProcessOutputDisplay;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

public class PEACEInstaller extends SwingWorker<Void, Void> {
	/**
	 * The constructor lays out the main components in this tab.
	 * 
	 * @param server The server to be installed.
	 */
	public PEACEInstaller(Server server, ServerSession serverSession) {
		// Save the server reference.
		this.server = server;
		this.serverSession = serverSession;

		// Create the styled document used to display output
		// from the installer process.
		pod = new ProcessOutputDisplay();

		// Create the top-level container.
		container = new JPanel(new BorderLayout(0, 0));
		container.setOpaque(false);

		// Create and add the top panel.
		createTopPanel();

		// Setup the text pane to display output logs.
		createOutputPane();
	}

	/**
	 * Method to obtain the top-level GUI container for this installer.
	 * 
	 * @return The top-level GUI container that contains all the GUI elements
	 *         associated with this installer.
	 */
	public JPanel getPanel() {
		return container;
	}

	/**
	 * This is a refactored utility method to create informational labels at the top
	 * of this panel. This method was introduced to keep the code clutter in the
	 * constructor to a minimum. This method is called only once from the
	 * constructor.
	 */
	private void createOutputPane() {
		// Add the text pane to a scroll pane to handle large outputs
		JScrollPane jsp = new JScrollPane(pod.getTextPane());
		// Add the scroll pane to this tab as well.
		container.add(jsp, BorderLayout.CENTER);
	}

	/**
	 * This is a refactored utility method to create informational labels at the top
	 * of this panel. This method was introduced to keep the code clutter in the
	 * constructor to a minimum.
	 */
	private void createTopPanel() {
		// Create the fields and progress bar.
		JTextField hostName = new JTextField(server.getName());
		hostName.setBackground(Color.white);
		hostName.setBorder(new EtchedBorder());
		hostName.setEditable(false);
		hostName.setAlignmentX(0);

		JTextField instPath = new JTextField(server.getInstallPath());
		instPath.setEditable(false);
		instPath.setAlignmentX(0);
		instPath.setBackground(Color.white);
		instPath.setBorder(new EtchedBorder());

		// Pack the host name and path in a horizontal box
		Box topBox = Box.createHorizontalBox();
		topBox.add(Box.createHorizontalStrut(10));
		topBox.add(new JLabel("Server name: "));
		topBox.add(hostName);
		topBox.add(Box.createHorizontalStrut(20));
		topBox.add(new JLabel("Install path: "));
		topBox.add(instPath);
		topBox.add(Box.createHorizontalGlue());

		// Now create a progress bar to indicate we are doing some operation here.
		Box botBox = Box.createHorizontalBox();
		botBox.add(Box.createHorizontalStrut(10));
		botBox.add(new JLabel("Progress indicator: "));
		progressBar = new JProgressBar();
		botBox.add(progressBar);
		botBox.add(Box.createHorizontalGlue());

		// Now pack the botBox and topBox into another box.
		Box bigBox = Box.createVerticalBox();
		bigBox.add(Box.createVerticalStrut(5));
		bigBox.add(topBox);
		bigBox.add(Box.createVerticalStrut(5));
		bigBox.add(botBox);
		bigBox.add(Box.createVerticalStrut(5));

		// Add the topPanel to this main panel
		container.add(bigBox, BorderLayout.NORTH);
	}

	/**
	 * Helper method to add <i>informational</i> output.
	 * 
	 * This is a helper method that is called from various methods to show more
	 * detailed progress information to the user. This method is synonymous to
	 * calling <code>addLog("something", "info");</code>
	 * 
	 * @param entry The log entry to be appended.
	 */
	protected void addLog(String entry) {
		addLog(entry, "info");
	}

	/**
	 * Helper method to add string to output and scroll it.
	 * 
	 * This is a helper method that is called from various methods to show more
	 * detailed progress information to the user.
	 * 
	 * @param entry The log entry to be appended.
	 * @param style The style (that determines color) for the entry.
	 */
	protected void addLog(String entry, String style) {
		try {
			DefaultStyledDocument output = pod.getDoc();
			// Add the output in given style
			output.insertString(output.getLength(), entry, output.getStyle(style));
			// Scroll the output to the bottom is automatically done by the
			// ProcessOutputDisplay class.
		} catch (BadLocationException e) {
			ProgrammerLog.log(e);
		}
	}

	@Override
	protected Void doInBackground() throws Exception {
		try {
			// Setup install status for server.
			server.setStatus(Server.ServerStatusType.INSTALLING);

			// Start up the roving progress bar.
			progressBar.setIndeterminate(true);

			// Connect to the remote server.
			addLog("Attempting to connect to server " + server.getName() + "...\n");
			serverSession.connect();
			addLog("Successfully connected to server.\n");

			// Create the destination directory.
			serverSession.mkdir(server.getInstallPath());
			addLog("Successfully created install directory '" + server.getInstallPath() + "'\n");

			// Copy the installation files over.
			String instCmd = copyInstallFiles(serverSession);

			// Add optional command line option to try and install EAST
			if (instCmd != null && server.hasEASTInstalled()) {
				instCmd += " --buildEAST";
			}

			int exitCode = 0;
			if (instCmd != null) {
				// Run the install command and display outputs.
				addLog("\nRunning install command: " + instCmd + "\n");
				exitCode = serverSession.exec(instCmd, pod.getDoc());
			}

			// Create working directory on the server for storing job and est data files.
			if (exitCode == 0) {
				String jobDir = server.getInstallPath() + "/jobs";
				String estDir = server.getInstallPath() + "/estData";

				// Create the directories.
				addLog("Creating Job directory: " + jobDir + "...");
				serverSession.mkdir(jobDir);
				addLog("Done\n");

				addLog("Creating EST storage directory: " + estDir + "...");
				serverSession.mkdir(estDir);
				addLog("Done\n");
			}

			// Log the final result.
			addLog("\nFinal exit code: " + exitCode + "\n");

			if (exitCode != 0) {
				throw new Exception("One of the installation steps did not complete successfully.");
			}

			// The installation was successful. Update the server status.
			server.setStatus(Server.ServerStatusType.GOOD);
		} catch (Exception e) {
			String exception = Utilities.toString(e);
			exception += "\n*** INSTALLATION FAILED ***\n";
			addLog(exception, "stderr");
			server.setStatus(Server.ServerStatusType.INSTALL_FAILED);
		} finally {
			// Stop the roving progress bar.
			progressBar.setIndeterminate(false);

			// Close the connection.
			serverSession.disconnect();
			addLog("Connection closed.\n");
		}
		return null;
	}

	/**
	 * This is a helper method that is used to copy the install files to the server.
	 * Specifically this method copies "peace.tar.gz" and "install.sh".
	 * 
	 * @param ss The server session to be used for copying the file(s) to the
	 *           server.
	 * 
	 * @throws Exception Exception is generated if errors occur while copying the
	 *                   install files.
	 * 
	 * @return The install file to be run on the remote machine.
	 */
	private String copyInstallFiles(ServerSession ss) throws Exception {
		boolean isWindows = ss.getOSType().equals(Server.OSType.WINDOWS);
		String[] installFiles = (isWindows ? WINDOWS_FILES : LINUX_FILES);

		// Copy the install file to the remote machine.
		String SrcPath = (isWindows ? "installFiles/windows/" : "installFiles/linux/");

		for (String file : installFiles) {
			// Get the input stream to the install file (directly from the PEACE jar).
			addLog("Transferring install file: " + file + "...");
			InputStream is = Utilities.getStream(SrcPath + file);
			ss.copy(is, server.getInstallPath(), file, 0700);
			addLog("Done.\n");
		}

		// Build the final command if supplied
		String cmd = null;
		if (installFiles[installFiles.length - 1] != null) {
			// For now assuming only Linux/Unix
			cmd = "cd " + server.getInstallPath() + " && " + installFiles[installFiles.length - 1];
		}
		return cmd;
	}

	/**
	 * This is the server entry that is being currently installed via this wizard.
	 * This object is initially used to establish connection and later on to change
	 * status.
	 */
	private final Server server;

	/**
	 * This is the server session that is active in ServerWizard
	 */
	private final ServerSession serverSession;

	/**
	 * Helper class to display styled output from the installation process. This
	 * class contains an styled document along with a text pane to display the
	 * styled document.
	 */
	private final ProcessOutputDisplay pod;

	/**
	 * A roving progress bar to show the user we are doing some work. We really
	 * can't tell how much progress has happened but we can say we are doing
	 * something so that the user understands the GUI is not hanging.
	 */
	private JProgressBar progressBar;

	/**
	 * The top-level container that logically contains all the GUI elements
	 * associated with this installer. This container is created in the constructor.
	 */
	private final JPanel container;

	/**
	 * The list of files to be installed on a Linux machine. Note that the last
	 * entry is the command to be executed for completing the installation. If the
	 * installation does not require running a command, then make the last entry
	 * null.
	 */
	private static final String[] LINUX_FILES = { "peace.tar.gz", "install.sh", "listJobs.sh", "./install.sh" };

	/**
	 * The list of files to be installed on a Windows machine. Note that the last
	 * entry is the command to be executed for completing the installation. If the
	 * installation does not require running a command, then make the last entry
	 * null.
	 */
	private static final String[] WINDOWS_FILES = { "License.txt", "ReadMe.txt", "peace.exe", "launcher.exe",
			"listJobs.bat", null };
}
