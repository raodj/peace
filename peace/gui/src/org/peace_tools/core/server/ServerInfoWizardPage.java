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

package org.peace_tools.core.server;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;

import org.peace_tools.core.FileInfo;
import org.peace_tools.core.session.RemoteFileSystemView;
import org.peace_tools.core.session.RemoteServerSession;
import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;

/**
 * This class serves as the second interactive page in a ServerWizard.
 * This page permits the user to enter additional information about
 * the server on which PEACE is to be configured and installed.
 * When the "Next >" button is clicked this wizard verifies that the
 * directory for installing PEACE does not exist on the remote 
 * machine. 
 * 
 * <p><b>Note:</b>  This page does the same operation for both local 
 * and remote servers.</p>
 */
public class ServerInfoWizardPage extends GenericWizardPage 
implements Runnable, ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a text area for
	 * description, a text field for install directory, and a formatted
	 * field for polling time delay.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * @param server The in-memory server object that contains the
	 * core information about the server.
	 * @param lockPath Flag to indicate if this page must permit the
	 * user to edit the install path.
	 */
	public ServerInfoWizardPage(ServerWizard wizard, Server server,
			boolean lockPath) {
		this.wizard   = wizard;
		this.server   = server;
		this.lockPath = lockPath; 
		// Setup the title(s) for this page and border
		setTitle("Server Information", 
				"Set additional information about the server");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create panels with description, install folder, and
		// polling time.
		description = new JTextArea(3, 3);
		JScrollPane jsp = new JScrollPane(description);
		jsp.setMinimumSize(description.getPreferredSize());
		
		JComponent descBox = 
			Utilities.createLabeledComponents("Description for server:",
					"(This is for your reference & can be anything)", 0, false, 
					jsp);
		// Put the install path text file and browse button into a 
		// single horizontal panel.
		browse = Utilities.createButton(null, " Browse ", 
				"Browse", this, "Browse local file system", true);
		installPath = new JTextField(10);
		JPanel horizBox = new JPanel(new BorderLayout(0, 10));
		horizBox.add(installPath, BorderLayout.CENTER);
		horizBox.add(browse, BorderLayout.EAST);
		
		JComponent dirBox =
			Utilities.createLabeledComponents("Enter install directory (absolute path):",
					"(Directory must not exist for fresh installs)",
					0, false, horizBox);
		JComponent timeBox =
			Utilities.createLabeledComponents("Enter polling delay (seconds):",
					"(Delay between successive checks for job completion)", 4,
					false, (pollTime = new JSpinner(new SpinnerNumberModel(30, 10, 120, 1))));
		
		// Set up additional informational panel if the install
		// directory is editable for verification.
		JPanel labelPanel = null;
		if (!lockPath) {
			// Create the informational panel.
			labelPanel = createInfoPanel();
		} else {
			installPath.setEnabled(false);
		}

		// Pack the options along with a pretty icon and label into 
		// a sub-panel.
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
			descBox, Box.createVerticalStrut(10),
			dirBox, Box.createVerticalStrut(10),
			timeBox, Box.createVerticalStrut(10),
			labelPanel);
		// Set up the border
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * This is a refactored utility method to create informational
	 * labels. This method was introduced to keep the code clutter in
	 * the constructor to a minimum. 
	 * 
	 * @return A panel to which the labels have been added.
	 * 
	 */
	private JPanel createInfoPanel() {
		// Let the user know the remote directory will be validated when
		// they click the "Next>" button.
		fixedMsgs[0] = new JLabel("<html>Install path will be verified when " + 
				"the<br>Next button is clicked</html>", 
				Utilities.getIcon("images/16x16/Information.png"),
				JLabel.LEFT);
		// The fixedMsg to inform user to "wait" is a bit more involved.
		progressBar = new JProgressBar();
		JLabel waitMsg = new JLabel("Install path is being verified. " +
			"Please wait...", Utilities.getIcon("images/16x16/HourGlass.png"),
			JLabel.CENTER);
		JPanel waitMsgPanel = new JPanel();
		waitMsgPanel.setLayout(new BoxLayout(waitMsgPanel, BoxLayout.Y_AXIS));
		waitMsgPanel.add(waitMsg);
		waitMsgPanel.add(Box.createVerticalStrut(10));
		waitMsgPanel.add(progressBar);
		fixedMsgs[1] = waitMsgPanel;
		fixedMsgs[1].setBorder(new CompoundBorder(new EtchedBorder(EtchedBorder.LOWERED),
				new EmptyBorder(5, 30, 5, 30)));
		fixedMsgs[1].setVisible(false);
		// Ensure that the two messages take up the same space so that
		// the display switching is smooth.
		int height = fixedMsgs[1].getPreferredSize().height;
		height    -= fixedMsgs[0].getPreferredSize().height;
		fixedMsgs[0].setBorder(new EmptyBorder(height / 2, 5, height / 2, 5));
		// Pack the various Box'es into the main panel. For this
		// we need to set the x-alignment correctly.
		fixedMsgs[0].setAlignmentX(0);
		fixedMsgs[1].setAlignmentX(0);
		// Create and add the messages to the subPanel.
		JPanel subPanel = new JPanel(new BorderLayout(0, 0));
		subPanel.add(fixedMsgs[0], BorderLayout.NORTH);
		subPanel.add(fixedMsgs[1], BorderLayout.SOUTH);
		return subPanel;
	}
		
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the data being displayed in
	 * the GUI fields from the data stored in the in-memory
	 * Server object.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// First update the necessary information.
		description.setText(server.getDescription());
		installPath.setText(server.getInstallPath().trim());
		// Setup default install path for remote server.
		if ((installPath.getText().length() == 0) && 
			(server.isRemote())) {
			installPath.setText("/home/" + server.getUserID() + "/PEACE");
		}
		pollTime.setValue(new Integer((int) server.getPollTime()));
	}

	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		// Save the information for future reference.
		server.setDescription(description.getText());
		server.setPollTime(((Number) pollTime.getValue()).longValue());
		server.setInstallPath(installPath.getText().trim());
		if (lockPath) {
			// The install path was uneditable. So assume it is good.
			return true;
		}
		
		// Before we proceed further, we verify that the install path
		// that the user has specified is valid and unique. First do some
		// very basic checks on the path.
		String instPath = installPath.getText().trim();
		if (instPath.length() < 2) {
			JOptionPane.showMessageDialog(this, "The installation path " +
					"is either empty or too short.", "Invalid Install Path", 
					JOptionPane.ERROR_MESSAGE);
			return false;
		}
		// Check to ensure that the install path does not have spaces in it.
		if (instPath.contains(" ")) {
			// Spaces cause issues with scripts and batch files. So we
			// avoid using it.
			JOptionPane.showMessageDialog(this, NO_SPACE_MSG,
					"Spaces in install path", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		// Check if the path specified is absolute windows or Linux path.
		// We need a more reliable check than this..
		boolean absPath = ((":\\".equals(instPath.substring(1, 3)) || 
				instPath.charAt(0) == '/'));
		if (!absPath) {
			JOptionPane.showMessageDialog(this, "The installation path " +
					"must be an absolute Windows or Linux path.", 
					"Invalid Install Path",	JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		// Disable any changes to this page for now.
		Utilities.setEnabled(this, false);
		// But enable the label with message...
		fixedMsgs[0].setVisible(false);
		fixedMsgs[1].setVisible(true);
		Utilities.setEnabled(fixedMsgs[1], true);
		progressBar.setIndeterminate(true);
		validate();
		// Disable the buttons in in the wizard while this operation
		// is underway as we are going to spin off a separate thread
		// for this operation.
		dialog.setButtonStatus(0, 0, -1);

		// Create and establish a connection with a remote server
		// to help
		serverSession = SessionFactory.createSession(dialog, server);
		// Launch the connect call from a separate thread.
		// Start the connection thread.
		Thread connThread = new Thread(this);
		// Flag this thread as a daemon because some connections hang
		// around for a long time and we don't want this thread to.
		connThread.setDaemon(true);
		wizard.addThread(connThread);
		connThread.start();
		// Return false indicating that wizard must not change pages
		// for now. Later on once the connection thread has finished
		// its job, we will let the wizard change pages.
		return false;
	}

	/**
	 * This method checks to ensure that the install directory does not
	 * exists but can be created.  This method was primarily introduced
	 * to streamline the {@link #run()} method in this class. This 
	 * method is called just one from {@link #run()} method.
	 * 
	 * @throws Exception This method throws an exception if the install
	 * directory exists or if an error occurs when attempting to
	 * create it.
	 */
	private void checkInstallDir() throws Exception {
		final String instPath = server.getInstallPath();
		FileInfo pathInfo     = serverSession.fstat(instPath);
		if (pathInfo.exists()) {
			// Error when running remote command. Bail out.
			throw new Exception("The install path already exists on " +
				"the server.\nPlease choose a different install path.");				
		}
		try {
			// Ensure we can create the directory.
			serverSession.mkdir(instPath);
			// OK, we can create directory. For now delete it.
			serverSession.rmdir(instPath);
		} catch (Exception e) {
			Exception err = new Exception("Unable to create the specified install " +
			"path on the server.\nPlease choose a different install path.");
			err.initCause(e);
			throw err;
		}
	}
	
	/**
	 * A refactored utility method that checks to ensure a Linux/Unix machine
	 * has the basic set of development tools installed. This method is
	 * called only once from the {@link #run()} method. This method was
	 * introduced to streamline the {@link #run()} method.
	 * 
	 * @param streamsData The data returned by calling the command. This information
	 * is passed back to the caller for display.
	 * 
	 * @throws Exception This method throws an exception if the server does not
	 * have the necessary tools installed on it.
	 */
	private void checkBasicLinuxTools(String streamsData[]) throws Exception {
		final String cmd = "/usr/bin/which automake autoconf tar gzip make grep";
		if ((serverSession.exec(cmd, streamsData) != 0) || 
				(streamsData[0] == null) ||	(streamsData[1].length() > 0)) {
			throw new Exception("The server does not seem to have " +
					"the necessary software tools.\nPEACE runtime cannot be installed " +
					"without the necessary tools.\nThe software required are:\n" + cmd);
		}
	}

	/**
	 * A refactored utility method that checks to see if a Linux/Unix machine
	 * has MPI and compiler installed. After detection it informs the user and
	 * permits the user to abort installation (by throwing an exception) if
	 * the user does not want to proceed without MPI. This method is called 
	 * only once from the {@link #run()} method. This method was
	 * introduced to streamline the {@link #run()} method.
	 * 
	 * @param streamsData The data returned by calling the command. This information
	 * is passed back to the caller for display.
	 * 
	 * @throws Exception This method throws an exception if an error occurs
	 * or if the user decided to abort install (due to absence of MPI).
	 */
	private void checkMPICompiler(String streamsData[]) throws Exception {
		String cmd = "mpicxx --version";
		if ((serverSession.exec(cmd, streamsData) != 0) || (streamsData[0] == null) || 
				(streamsData[1].length() > 0)) {
			String execOutput = "Result of running command `mpicxx --version` on server returned:\n\n" +
				"ERROR STREAM:\n" + streamsData[1] + "\n\n" + "OUTPUT STREAM:\n" + streamsData[0];
			JPanel detailedMessage = Utilities.collapsedMessage(NO_MPI_MSG, execOutput); 
			int choice = JOptionPane.showConfirmDialog(wizard, 
					detailedMessage, "Unable to detect MPI", 
					JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
			if (choice == JOptionPane.NO_OPTION) {
				throw new Exception("The remote server does not seem to have " +
						"MPI installed on it.\nPEACE clustering engine install " +
				"cancelled by user.");
			}
		}
	}

	/**
	 * A refactored utility method that checks to see if a given server
	 * has GCC compiler installed (for building EAST). After detection 
	 * it informs the user and permits the user to abort installation 
	 * (by throwing an exception) if the user does not want to proceed
	 * without EAST. This method is called only once from the {@link #run()}
	 * method. This method was introduced to streamline the {@link #run()} method.
	 * 
	 * @param streamsData The data returned by calling the command. This information
	 * is passed back to the caller for display.
	 * 
	 * @throws Exception This method throws an exception if an error occurs
	 * or if the user decided to abort install (due to absence of MPI).
	 */
	protected void checkGCC(String streamsData[]) throws Exception {
		String cmd = "g++ --version";
		if ((serverSession.exec(cmd, streamsData) != 0) || (streamsData[0] == null) || 
				(streamsData[1].length() > 0)) {
			int choice = JOptionPane.showConfirmDialog(wizard, 
					NO_EAST_MSG, "Cannot install EAST (no g++)", 
					JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
			if (choice == JOptionPane.NO_OPTION) {
				throw new Exception("The server does not seem to have " +
						"g++ installed on it.\nPEACE clustering & assembler install " +
				"aborted.");
			}
		} else {
			// Found gcc. Set flag for future use
			wizard.setHaveGCC(true);
		}
	}

	/**
	 * A refactored utility method that checks to see if a Linux/Unix machine
	 * has AMOS tools installed on it. 
	 * 
	 * <p>AMOS (see http://sourceforge.net/apps/mediawiki/amos/) is a third-party 
	 * tool set that is used to convert ACE output file format generated by 
	 * EAST assembler to SAM file format. Similar to ACE, SAM is another widely 
	 * used text file format for storing assembly
	 * results (see http://samtools.sourceforge.net/). The process of converting
	 * an ACE file format to SAM is documented at: 
	 * http://sourceforge.net/apps/mediawiki/amos/index.php?title=Bank2contig#SAM_Conversion
	 * </p> 
	 * 
	 * <p><b>Note</b>We do not distribute AMOS tools with PEACE and we expect
	 * the user (or their system administrator) to install and setup
	 * AMOS tools for us. However, in this method we simply detect to see if
	 * the following three AMOS-tools that we use are available by running
	 * the following commands to simply return their version information:
	 * 
	 * <ul>
	 * <li>toAmos -v</li>
	 * <li>bank-transact -v</li>
	 * <li>bank2contig -v</li>
	 * </ul>
	 * 
	 * </p>    
	 * 
	 * <p>This method runs the above commands to detect the tools. If detection
	 * fails it informs the user and  permits the user to abort installation
	 * (by throwing an exception) if the user does not want to proceed without
	 * AMOS tools. This method is called only once from the {@link #run()} method. 
	 * This method was introduced to streamline the {@link #run()} method.</p>
	 * 
	 * <p><b>Note</b>: This method sets the AMOS tool availability in the server
	 * object owned by the wizard.</p>
	 * 
	 * @param streamsData The data returned by calling the command. This information
	 * is passed back to the caller for display.
	 * 
	 * @throws Exception This method throws an exception if an error occurs
	 * or if the user decided to abort install (due to absence of AMOS tools).
	 */
	private void checkAMOSTools(String streamsData[]) throws Exception {
		// Set output streams data to empty content to start with.
		streamsData[0] = "";
		streamsData[1] = "";
		final String commands[] = {"toAmos -v", "bank-transact -v", "bank2contig -v"};
		// Run the above three commands on the server to ensure we have those three
		// tools operational. Collate the outputs in the streamsData array. Stop
		// checking when a given command does not complete successfully.
		for(int i = 0; (i < commands.length); i++) {
			// Array to hold result from each command
			String cmdStreams[] = {"", ""};
			final int exitStatus = serverSession.exec(commands[i], cmdStreams);
			// Check to ensure command ran successfully. Note anomalous behavior that
			// AMOS tools write version to standard error rather than standard output!
			if ((exitStatus != 0) || (cmdStreams[1] == null) || (cmdStreams[0].length() > 0)) {
				// Command did not complete successfully! Create detailed error message
				String execOutput = "Result of running command " + 
					commands[i] + " on server returned:\n" +
					"EXIT CODE: " + exitStatus + "\n" +
					"ERROR STREAM:\n" + cmdStreams[1] + "\n" + 
					"OUTPUT STREAM:\n" + cmdStreams[0];
				JPanel detailedMessage = Utilities.collapsedMessage(NO_AMOS_MSG, execOutput);
				// Show it to the user prompting check to cancel install
				final int choice = JOptionPane.showConfirmDialog(wizard, 
						detailedMessage, "Unable to detect AMOS Tools", 
						JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
				if (choice == JOptionPane.NO_OPTION) {
					throw new Exception("The remote server does not seem to have " +
							"AMOS Tools installed on it.\nPEACE clustering engine install " +
							"cancelled by user.");
				} else {
					// The user want's to continue without AMOS tools.
					return;
				}
			} else {
				// Command completed successfully. Save outputs in main stream and
				// move onto the next command to check.
				String userInfo = commands[i] + "\n" + cmdStreams[1];
				streamsData[0] += userInfo;
			}
		}
		// When control drops here that means AMOS tools were readily
		// available. So set the flag in the server 
		wizard.getServer().setHasAMOSTools(true);
	}
	
	/**
	 * This method is invoked from a separate thread from the pageChanging()
	 * method. This method performs the task of attempting to connect to the
	 * remote server using the credentials supplied by the user and validating
	 * the install directory to ensure it does not exist.
	 * 
	 * <p><b>Note:</b>  The process of connecting to the remote machine is performed from a
	 *       separate thread because of the way Ganymede SSH has been
	 *       implemented. The call backs provided by Ganymede do not occur on
	 *       the Swing's EventDispatch thread (but from an generic thread).
	 *       Consequently, the GUI will be unresponsive (in most cases call
	 *       backs result in a blank window being displayed and the whole GUI
	 *       hangs) if the main EventDispatch thread is blocked. Consequently,
	 *       the connection is performed from a separate thread. When this
	 *       thread completes, it posts the necessary information back to the
	 *       main EventDispatch thread for updates.</p>
	 */
	public void run() {
		// Common exception to be reported to GUI (if any)
		Exception exp = null;
		// Array to hold result from remote command
		String streamsData[] = {"", ""};
		try {
			// Try to connect (with interactive prompts if
			// needed) to the remote server.
			serverSession.connect();
			// Now that we have connected, check if the directory does
			// not exist.
			checkInstallDir();
			// Ensure that the necessary tools are installed for Linux/
			// Unix family of servers.
			if (!Server.OSType.WINDOWS.equals(serverSession.getOSType())) {
				// Check if Unix/Linux server has basic development tools installed.
				checkBasicLinuxTools(streamsData);
				// OK, the basic stuff checks out. Next check if we have mpicc.
				String prevOutput = "STANDARD TOOLS DETAILS:\n" + streamsData[0];
				// Check to see if the server has MPI installed. This method
				// throws an exception if user decides to abort install.
				checkMPICompiler(streamsData);
				// Add current output with previous output.
				prevOutput = prevOutput + "\nMPI COMPILER DETAILS:\n" + streamsData[0]; 
				// Check to see if we have GCC for building EAST. The
				// gcc flag in the parent server wizard class is set by
				// the following method call.
				checkGCC(streamsData);
				prevOutput = prevOutput + "\nGCC (for EAST) DETAILS:\n" + streamsData[0];
				// Check to see if we have AMOS tools for converting ACE output
				// from EAST to SAM output file format.
				checkAMOSTools(streamsData);
				streamsData[0] = prevOutput + "\nAMOS TOOLS DETAILS:\n" + streamsData[0];
			} else {
				// Set a default detail so the user does not freak out.
				streamsData[0] = "No special tools are needed for windows servers.\n" +
						"PEACE GUI comes with pre-built runtime executable files.";
			}
			// The connection was successful! Credentials 
			// are valid and work just fine. Disconnect for now.
			serverSession.disconnect();
		} catch (Exception e) {
			// Save exception to be reported below.
			exp = e;
			// Close the connection (if any)
			serverSession.disconnect();
		}
		// Use a final variable so that we can access the exception
		// in the run() method below.
		final Exception result = exp;
		final JPanel    page   = this;
		final String    output = streamsData[0];
		// When control drops here we have either successfully
		// verified install path (no exceptions) or not.
		// Report the completion on the main AWTThread.
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				// First re-enable the buttons
				wizard.setButtonStatus(1, 1, -1);
				// In addition re-enable the input fields so user
				// can change it if desired.
				Utilities.setEnabled(page, true);
				// But disable the label with message...
				progressBar.setIndeterminate(false);
				fixedMsgs[1].setVisible(false);
				fixedMsgs[0].setVisible(true);
				validate();
				// Disable path control as configured.
				installPath.setEnabled(!lockPath);
				// Check and print any exceptions we may have had.
				if (result != null) {
					// The connection to the remote server failed.
					// Report it to the user.
					JPanel msg = Utilities.collapsedMessage(ErrorMsg, 
							Utilities.toString(result)); 
					JOptionPane.showMessageDialog(wizard, msg,
						"Invalid Information", JOptionPane.ERROR_MESSAGE);
				} else {
					JPanel msg = 
						Utilities.collapsedMessage("<html>" +
								"The install path does not exist (this is good) and will be<br>" +
								"created when PEACE runtime is installed. In addition, the<br>" +
								"server has all the necessary tools to build PEACE runtime.</html>",
								output);
					JOptionPane.showMessageDialog(wizard, msg,
						"Install Validation Success", JOptionPane.INFORMATION_MESSAGE);
					// In this case, simply let the user move onto the
					// next page.
					wizard.changePage(wizard.getCurrentPage(), 
							wizard.getCurrentPage() + 1);
				}
			}
		});
		// Remove the thread from the list of workers.
		wizard.removeThread(Thread.currentThread());
	}

	/**
	 * This method handles the tasks performed when the "Browse"
	 * button is clicked. It displays a JFileChooser and permits
	 * the user to select a local folder.
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		// For remote server we need to setup a suitable RemoteFileSystemView
		// for the JFileChooser to use.
		try {
			RemoteFileSystemView rfsv = null;
			if (server.isRemote()) {
				RemoteServerSession rss = (RemoteServerSession) SessionFactory.createSession(this, server);
				rfsv = new RemoteFileSystemView(rss);
			}
			// Setting file system view after constructor does not seem to work well..
			final JFileChooser jfc = (rfsv == null ? new JFileChooser() : new JFileChooser(rfsv));
			jfc.setDialogTitle("Choose install directory");
			jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			if (jfc.showDialog(this, "Install Here") == JFileChooser.APPROVE_OPTION) {
				// Copy the chosen directory to the work space combo box.
				String wsPath = jfc.getSelectedFile().getAbsolutePath();
				// Set the selected item in this path.
				installPath.setText(wsPath);
			}
		} catch (Exception e) {
			ProgrammerLog.log(e);
		}
	}

	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final ServerWizard wizard;
	
	/**
	 * Information about the actual server entry being edited. This
	 * reference is set when this class is instantiated and is never
	 * changed during the life time.
	 */
	private final Server server;
			
	/**
	 * Field to read/display the install path where PEACE is (or has
	 * been) installed on the given server.
	 */
	private JTextField installPath;
	
	/**
	 * Field to read and edit a brief description about the server.
	 * This information can be anything the user desires and is
	 * meaningful only to the user.
	 */
	private JTextArea description;
	
	/**
	 * Field to read/display the polling time.
	 */
	private JSpinner pollTime;
	
	/**
	 * These labels are created in the constructor but one of them
	 * remains hidden. The visibility of the labels is swapped
	 * just before and after the verification process. The idea is 
	 * let the user know that the wizard is doing some operation.
	 */
	private JComponent fixedMsgs[] = new JComponent[2];
	
	/**
	 * A roving progress bar to show the user we are doing 
	 * some work. We really can't tell how much progress has 
	 * happened but we can say we are doing something so that
	 * the user understands the GUI is not hanging.
	 */
	JProgressBar progressBar;
	
	/**
	 * If this flag is set to true, then this wizard permits the
	 * user only to change the install path on the server.
	 */
	private boolean lockPath;
	
	/**
	 * This instance variable is used to manage the connection to 
	 * the server just for testing purposes. Possibly this
	 * connection can be centralized into the ServerWizard for use
	 * throughout the wizard.
	 */
	private ServerSession serverSession;
	
	/**
	 * The browse button to be enabled only for local installs.
	 */
	private JButton browse;
	
	/**
	 * Just a simple error message that is displayed when install
	 * on a server/remote host could not be successfully verified.
	 */
	private static final String ErrorMsg = "<html>" + 
		"Unable to proceed further with the provided server information.<br/>" +
		"Either the installation path is not valid or the server does<br/>" +
		"not have the necessary software tools installed on it.<br/>" + 
		"<b>See details below for full details.</b><br>" +
		"Contact your system administrator to ensure all the necessary<br>" +
		"software tools are installed and the server is operational." +
		"</html>";

	/**
	 * A warning message that is displayed when mpicc is unavailable 
	 * on the remote machine on which PEACE is being installed.
	 */
	private static final String NO_MPI_MSG = "<html>" + 
		"<b>Could not successfully detect the MPI C++ compiler (mpicxx)</b><br>" +
		"The sources of this issue include:<ul>"+
		"<li>mpicxx was not found in the default path<br>" +
		"You need to set path in <b>.bashrc</b> for bash and in <b>.cshrc</b> for tcsh.</li>" +
		"<li>Possibly MPI is not installed on the machine or your module is<br>" +
		"not loaded. You need to contact your system adminstrator<br>" +
		"to determine if MPI is setup and if your module/path is not correct. Some<br>" +
		"clusters require MPI module(s) to be loaded manually, which your<br>" +
		"administrator can automate so that PEACE finds mpicc correctly.</li>" +
		"<li>There is a minor issue with MPI setup and you could proceed<br>" +
		"to see if further detection & configuration resolve this issue.</li>" +
		"<li>Note that PEACE uses a non-interactive shell to communicate<br>" +
		"with remote hosts. Therefore ensure your setup is consistent<br>" + 
		"in both interactive and non-interactive logins/shells.</li>" +
		"</ul><br>" +
		"You may still procceed further with the default non-MPI based C++<br>" +
		"compiler. <i>However, PEACE will most likely not run in parallel mode.</i><br><br>" +
		"Do you wish to proceed with default C++ compiler (if available)?</html>";

	/**
	 * A warning message that is displayed when AMOS tool detection
	 * failed on the machine on which PEACE is being installed. This message
	 * is displayed by the {@link #checkAMOSTools(String[])} method.
	 */
	private static final String NO_AMOS_MSG = "<html>" + 
		"<b>Could not successfully detect AMOS tools</b><br/>" +
		"AMOS tools (see <a href=\"http://sourceforge.net/apps/mediawiki/amos\">http://sourceforge.net/apps/mediawiki/amos</a>)<br/>" +
		"are needed only by EAST to convert ACE files to SAM file format.<br/>" +
		"This is not a required component and you may proceed further.<br/>" +
		"However, assembly jobs run on this server will not be able to<br/>"+
		"directly generate SAM files (but will generate FASTA &amp; ACE files).<br/>" +
		"<br/>" +
		"Here are a few troubleshooting tips:<ul>"+
		"<li>AMOS tools may not be installed on the server.<br/>" +
		"Please contact your system adminstrator for details.</li>" +
		"<li>AMOS module may not be loaded automatically or path is incorrect.<br/>"+
		"You may need to appropriately modify <tt>.bashrc</tt>  or <tt>.cshrc</tt><br/>"+
		"based on the login shell you use on the server.</li>"+
		"<li>Note that PEACE uses a non-interactive shell to communicate<br/>" +
		"with remote hosts. Therefore ensure your setup is consistent<br/>" + 
		"in both interactive and non-interactive logins/shells.</li>" +
		"</ul><br/>" +
		"You may still procceed further without SAM file generation feature.<br/>" +
		"Do you wish to proceed further?</html>";

	/**
	 * A simple message that is displayed to the user to indicate
	 * we don't accept spaces in installation path.
	 */
	private static final String NO_SPACE_MSG = "<html>" + 
		"The installation path has spaces in it. Spaces cause<br>" +
		"several compatibility problems and are hard to work around.<br>" +
		"Consequently, in the current version of PEACE spaces cannot<br>" +
		"be present in the installation path.<br><br>" +
		"<b>Please choose a path without spaces in it.<b></html>";
	
	/**
	 * A warning message that is displayed when gcc is unavailable 
	 * on a given machine/server and consequently EAST cannot be installed.
	 */
	private static final String NO_EAST_MSG = "<html>" + 
		"The C++ compiler from GNU Compiler Collection (GCC) (i.e., g++)<br>" +
		"was not found in the default path. Currently, GCC is needed for<br>" +
		"building EAST and consequently EAST will not be installed. <br><br>" +
		"However, you may still proceed further with installation<br>" +
		"of the PEACE clustering engine.</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
