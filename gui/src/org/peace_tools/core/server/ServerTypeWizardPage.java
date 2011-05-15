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
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JProgressBar;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Server.OSType;

/**
 * This class serves as the first interactive page in a ServerWizard.
 * This page permits the user to select the type of server entry
 * to be added. In addition, if the selected server is a remote
 * server then this wizard page permits the user to set the host name
 * and credentials for logging onto the remote server. When the 
 * "Next >" button is clicked then this wizard verifies that it
 * can connect with the given credentials. 
 */
public class ServerTypeWizardPage extends GenericWizardPage 
implements ActionListener, Runnable {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include a combo box to
	 * choose between local and remote servers. In addition, there is
	 * a separate panel that is setup for credentials associated
	 * with the remote host. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * @param server The in-memory server object that contains the
	 * core information about the server.
	 * @param changeCredentialsOnly Flag to indicate if this page
	 * must permit the user to just change the credentials.
	 */
	public ServerTypeWizardPage(WizardDialog wizard, Server server,
			boolean changeCredentialsOnly) {
		this.wizard = wizard;
		this.server = server;
		this.changeCredentialsOnly = changeCredentialsOnly; 
		// Setup the title(s) for this page and border
		setTitle("Server Type", "Set the type of server entry to add");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create a panel with the type combo box for the user to
		// choose from two different options.
		String[] options = {"Local Server (localhost)",
				"Remote Server (access via SSH)"};
		serverTypes = new JComboBox(options);
		serverTypes.setEditable(false);
		serverTypes.addActionListener(this);
		// Pack the options along with a pretty icon and label into 
		// a sub-panel.
		JPanel subPanel = new JPanel(new BorderLayout(2, 2));
		subPanel.setBorder(new EmptyBorder(10, 15, 10, 10));
		subPanel.add(new JLabel("Select type of server to add:"), BorderLayout.NORTH);
		subPanel.add(serverTypes, BorderLayout.SOUTH);
		// Add sub-panel to the main panel
		add(subPanel, BorderLayout.NORTH);
		
		// Create a sub panel with host name, port, login ID, and 
		// password for accessing a remote machine.
		JComponent srvrNameBox = 
			Utilities.createLabeledComponents("Enter name or IP Address of server:",
				null, 4, false, (hostName = new JTextField(10)));
		// Create the spinner to set non-default port number
		JComponent portBox =
			Utilities.createLabeledComponents("Port:", null, 4, false, 
					(port = new JSpinner(new SpinnerNumberModel(22, 1, 65535, 1))));
		Utilities.adjustDimension(portBox, -25, 0);
		// Create a box to hold the server name and port number in
		// a horizontal fashion.
		JComponent serverBox = Box.createHorizontalBox();
		serverBox.add(srvrNameBox);
		serverBox.add(Box.createHorizontalStrut(20));
		serverBox.add(portBox);

		// Utilities.adjustDimension(hostName, 200, 4);
		// Create the user-name and password fields
		JComponent userNameBox = 
			Utilities.createLabeledComponents("Login (user) ID:", null,
				4, false, (userName = new JTextField(10)));
		JComponent passwordBox = 
			Utilities.createLabeledComponents("Password:", null,
				4, false, (password = new JPasswordField(10)));
		// Wrap the credential inputs appropriately
		JComponent credentialBox = Box.createHorizontalBox();
		credentialBox.add(userNameBox);
		credentialBox.add(Box.createHorizontalStrut(30));
		credentialBox.add(passwordBox);
		// Let the user know the credentials will be validated when
		// they click the "Next>" button.
		fixedMsgs[0] = new JLabel("<html>Login credentials will be verified when " + 
				"<br>the Next button is clicked.<br></html>", 
				Utilities.getIcon("images/16x16/Information.png"),
				JLabel.LEFT);
		// The fixedMsg to inform user to wait is a bit more involved.
		progressBar = new JProgressBar();
		JLabel waitMsg = new JLabel("Credentials are being verified. " +
			"Please wait...", Utilities.getIcon("images/16x16/HourGlass.png"),
			JLabel.CENTER);
		JPanel waitMsgPanel = new JPanel(new BorderLayout(0, 3));
		waitMsgPanel.add(waitMsg, BorderLayout.NORTH);
		waitMsgPanel.add(progressBar, BorderLayout.SOUTH);
		fixedMsgs[1] = waitMsgPanel;
		fixedMsgs[1].setBorder(new CompoundBorder(new EtchedBorder(EtchedBorder.LOWERED),
				new EmptyBorder(5, 30, 5, 30)));
		Utilities.adjustFont(fixedMsgs[1], 0, 10, 1);
		fixedMsgs[1].setVisible(false);
		// Pack the various Box'es into a suitable panel.
		remoteSrvrInfo = Utilities.createLabeledComponents(null, null, 0, true,
				serverBox,
				Box.createVerticalStrut(10),
				credentialBox,
				Box.createVerticalStrut(10),
				fixedMsgs[0], fixedMsgs[1]);		
		remoteSrvrInfo.setBorder(new CompoundBorder(new TitledBorder(" Remote Server Data "),
				new EmptyBorder(10, 10, 10, 10)));
		// Add the remote server Information panel to the main page
		add(remoteSrvrInfo, BorderLayout.CENTER);
	}
		
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the data being displayed in
	 * the GUI fields from the data stored in the in-memory
	 * Server object.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Enable/disable the remoteSrvrInfo depending on the type of
		// server that has been currently selected. This done via the
		// action performed method by passing in an null event.
		this.actionPerformed(null);
	}

	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		// Save host name for local & remote machines
		server.setName(hostName.getText().trim());
		// Save port number 
		server.setPort(((Number) port.getValue()).intValue());
		// Before we proceed further to detect type of OS on the server, 
		// for remote server we validate the credentials supplied by the user --
 		// that is, verify if the credentials work to connect to the
		// remote machine. First update the credentials. 
		server.setUserID(userName.getText().trim());
		server.setPassword(password.getText());
		// Disable any changes to this page for now.
		Utilities.setEnabled(this, false);
		// But enable the label with message...
		fixedMsgs[0].setVisible(false);
		fixedMsgs[1].setVisible(true);
		Utilities.setEnabled(fixedMsgs[1], true);
		progressBar.setIndeterminate(true);
		validate();
		// Create and establish a connection with the server
		srvrSession = SessionFactory.createSession(dialog, server);
		// Disable the buttons in in the wizard while this operation
		// is underway as we are going to spin off a separate thread
		// for this operation (to handle remote server cases as well).
		dialog.setButtonStatus(0, 0, -1);
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
	 * This function checks to see if a given remote server is a Linux
	 * server. Currently, PEACE+GUI only support Linux as the remote
	 * server. It should be possible to extend the feature support for
	 * Unix/Linux remote server. However, for now we support only Linux.
	 * This check used to be part of the {@link #run()} method. However
	 * since the original method was getting cluttered, this functionality
	 * has been refactored out into its own method.
	 * 
	 * <p><b>Note:</b>This method is to be invoked only from the {@link #run()}
	 * method after a connection has been established for remote servers.</p>
	 * 
	 * @return The standard output and error streams that contain the 
	 * Linux version for display to user.
	 * 
	 * @exception Exception This method throws an exception if the remote 
	 * machine is not a Linux machine.
	 */
	private String[] checkRemoteServer() throws Exception {
		String streamsData[] = new String[2];
		if ((srvrSession.exec("uname -a", streamsData) != 0) ||
				(streamsData[0] == null) || 
				(streamsData[0].indexOf("Linux") == -1)) {
				// Error when running remote command. Bail out.
				throw new Exception("The remote machine does not seem " +
					"to be a Linux machine.\nCurrently PEACE only supports " +
					"Linux clusters. Sorry for the inconvenience.\nIf this " +
					"is a feature you need then please file a feature request.\n" +
					"[Standard output: " + streamsData[0] + "]\n" + 
					"[Standard error : " + streamsData[1] + "]\n");				
			}
		return streamsData;
	}
	
	/**
	 * This is a refactored helper method to display the results of credentials
	 * checks and operating system information. This method was introduced to keep the
	 * code clutter in the {@link #run()} method to a minimum. This method is
	 * called only once from the {@link #run()} method to display the 
	 * results of the checks performed by it. 
	 * 
	 * <p><b>Note:</b> This method resets the wizard page (including various
	 * controls) back to their nominal operating status. Consequently,
	 * calling this method has some additional side effects!</p>
	 * 
	 * @param result The exception (if any) that was generated from the various
	 * checks performed to verify credentials and detect OS type.
	 * 
	 * @param page The wizard page that logically owns the dialog box to be
	 * displayed by this method.
	 * 
	 * @param version The OS version, type, and any other information to be 
	 * displayed to the user. This is typically just a simple multi-line text
	 * string.
	 */
	private void showResults(final Exception result, final JPanel page,
	                         final String version) {
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
		// Disable other controls if only credentials are to change
		serverTypes.setEnabled(!changeCredentialsOnly);
		hostName.setEnabled(!changeCredentialsOnly);

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
				Utilities.collapsedMessage("<html>Successfully connected " +
						"to the Local/Remote server<br>" +
						"using the supplied information.</html>",
						"Server machine is: " + version);
			JOptionPane.showMessageDialog(wizard, msg,
				"Connection Success", JOptionPane.INFORMATION_MESSAGE);
			// Save information in the in-memory format
			server.setRemote(serverTypes.getSelectedIndex() == 1);
			server.setName(hostName.getText());
			server.setUserID(userName.getText());
			server.setPassword(password.getText());
			// In this case, simply let the user move onto the
			// next page.
			wizard.changePage(wizard.getCurrentPage(), 
					wizard.getCurrentPage() + 1);
		}
	}
	
	/**
	 * This method is invoked from a separate thread from the
	 * pageChanging() method to detect OS type for server. 
	 * This method performs the following tasks:
	 * 
	 * <ol>
	 * <li>First it attempts to connect to the remote server using the 
	 * credentials supplied by the user. The process of 
	 * connecting to the remote machine is performed from a separate
	 * thread because of the way Ganymede SSH has been implemented.
	 * The call backs provided by Ganymede do not occur on the 
	 * Swing's EventDispatch thread (but from an generic thread). 
	 * Consequently, the GUI will be unresponsive (in most cases
	 * call backs result in a blank window being displayed and the
	 * whole GUI hangs) if the main EventDispatch thread is blocked.
	 * Consequently, the connection is performed from a separate</li>
	 * thread. <li>
	 * 
	 * <li>It detects and sets the OS type in the server entry to match
	 * the OS type of the target server.</li>
	 * </ol>
	 * 
	 * When this thread completes, it posts the necessary
	 * information back to the main EventDispatch thread for updates. 
	 */
	public void run() {
		// Common exception to be reported to GUI (if any)
		Exception exp = null;
		// Array to hold result from remote command
		String streamsData[] = new String[2];
		try {
			// Try to connect (with interactive prompts if
			// needed to a remote server).
			srvrSession.connect();
			// Now that we have connected, check if the remote
			// host is a Linux machine.
			if (server.isRemote()) {
				streamsData = checkRemoteServer();
			}
			// Setup the OS type for future reference.
			server.setOSType(srvrSession.getOSType());
			// The connection was successful! Credentials 
			// are valid and work just fine. Disconnect for now.
			srvrSession.disconnect();
		} catch (Exception e) {
			System.out.println(e);
			e.printStackTrace();
			// Save exception to be reported below.
			exp = e;
			// Close the connection (if any)
			srvrSession.disconnect();
		}
		// Use a final variable to that we can access the exception
		// in the run() method below.
		final Exception result = exp;
		final JPanel    page   = this;
		String    osInfo       = (server.isRemote() ? streamsData[0] + "\n" : "") +
			"Operating System (OS) = ";
			 
		// Add some disambiguation information for Unix type hosts
		if (OSType.UNIX.equals(server.getOSType())) {
			osInfo += "Unix/Unix-compatible OS\n" + 
				"(such as: OS-X, Solaris, etc.)";
		} else {
			osInfo += server.getOSType();
		}
		// When control drops here we have either successfully
		// connected (no exceptions) or we failed connecting.
		// Report the completion on the main AWTThread.
		final String version = osInfo;
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				// Do the operations in a refactored method.
				showResults(result, page, version);
			}
		});
		// Remove the thread from the list of workers.
		wizard.removeThread(Thread.currentThread());
	}
	
	@Override
	public void actionPerformed(ActionEvent arg0) {
		// This method is called whenever the user changes the
		// option between local and remote server. This method
		// enables/disables the remote server information.
		server.setRemote(serverTypes.getSelectedIndex() == 1);
		Utilities.setEnabled(remoteSrvrInfo, serverTypes.getSelectedIndex() == 1);
		// If local host is selected stuff name of host as server name.
		if (serverTypes.getSelectedIndex() == 0) {
			hostName.setText("localhost");
			userName.setText(System.getProperty("user.name"));
			password.setText("");
		} else {
			hostName.setText(server.getName());
			final int portNumber = (server.getPort() == -1 ? 22 : server.getPort());
			port.setValue(new Integer(portNumber));
			userName.setText(server.getUserID());
			password.setText(server.getPassword());
		}
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final WizardDialog wizard;
	
	/**
	 * Information about the actual server entry being edited. This
	 * reference is set when this class is instantiated and is never
	 * changed during the life time.
	 */
	private final Server server;
	
	/**
	 * The combo box that permits the user to choose between a local
	 * and remote server setting.
	 */
	private JComboBox  serverTypes;
	
	/**
	 * The panel that contains the data entry fields related to a
	 * remote server entry. This panel contains all the fields so
	 * that they can be enabled and disabled in one big swoop whenever
	 * the user changes between local and remote server settings. 
	 */
	private JPanel remoteSrvrInfo;
	
	/**
	 * Field to read/display the host name to the user. The host name
	 * maybe the FQN or the IP address of the host.
	 */
	private JTextField hostName;
	
	/**
	 * Field to read/display the port number set for remote servers. 
	 * This entry is enabled only for remote servers. The default
	 * port number is 22 (for ssh).
	 */
	private JSpinner port;
	
	/**
	 * Field to read/display the user/login ID. This information is
	 * meaningful for remote servers.
	 */
	private JTextField userName;
	
	/**
	 * Field to read/display the password ID. This information is
	 * meaningful for remote servers.
	 */
	private JTextField password;
	
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
	 * user only to change the credentials associated with a remote
	 * server.
	 */
	private boolean changeCredentialsOnly;
	
	/**
	 * This instance variable is used to manage the connection to 
	 * the remote server just for testing purposes. Possibly this
	 * connection can be centralized into the ServerWizard for use
	 * throughout the wizard.
	 */
	private ServerSession srvrSession;
	
	/**
	 * Just a simple error message that is displayed when connection
	 * to the remote host could not be successfully established.
	 */
	private static final String ErrorMsg = "<html>" + 
		"Failed to connect to the remote server. Either your<br>" +
		"supplied host name and credentials are incorrect or<br>" +
		"the remote host is not a Linux server compatible with PEACE.<br>" +
		"<br>Please verify that the information you provided is<br>" +
		"correct and try again." +
		"</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
