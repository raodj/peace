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

package org.peace_tools.core.session;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.InetAddress;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JProgressBar;
import javax.swing.JTextField;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.FileInfo;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

import ch.ethz.ssh2.Connection;
import ch.ethz.ssh2.KnownHosts;
import ch.ethz.ssh2.SFTPException;
import ch.ethz.ssh2.SFTPv3Client;
import ch.ethz.ssh2.SFTPv3FileAttributes;
import ch.ethz.ssh2.SFTPv3FileHandle;
import ch.ethz.ssh2.ServerHostKeyVerifier;
import ch.ethz.ssh2.Session;
import ch.ethz.ssh2.StreamGobbler;

/**
 * A remote server session based on the secure shell (SSH) protocol.
 * 
 * <p>This class provides an implementation of the ServerSession API.
 * Specifically, this class provides a session that can be used to
 * interact with a remote host via the secure shell (SSH) protocol.
 * The secure shell protocol is a defacto standard for interacting
 * with remote servers via the Internet today. It provides all the
 * necessary security features to safely interact with remote 
 * hosts and almost all super computing clusters mandate the use of
 * SSH for interactions. </p>
 * 
 * <p>This class uses the Ganymede SSH 
 * implementation for establishing SSH connections. Ganymede SSH
 * supports only ssh-2 protocol. Please refer to Genymede SSH website
 * for further details: <A HREF="http://www.ganymed.ethz.ch/ssh2/">
 * http://www.ganymed.ethz.ch/ssh2/</A>. PEACE distributes Ganymede
 * license file as per Ganymede licensing requirements.</p>
 * 
 */
public class RemoteServerSession extends ServerSession 
implements ServerHostKeyVerifier {
	/**
	 * The constructor merely passes parameters to the base class that
	 * initializes the instance variables. No special operations are
	 * performed in the constructor.
	 * 
	 * @param server The server data entry that provides the necessary
	 * information to connect to the server.
	 * @param parent The parent component that should be used to 
	 * create GUI elements that may be needed for any interactive
	 * operations.
	 */
	protected RemoteServerSession(Component parent, Server server) {
		super(server, parent);
		connection = null;
		osType     = null;
		purpose    = null;
	}

	/**
	 * Connect to the server in order to perform various operations. 
	 * 
	 * This method must be used to establish a connection to a server
	 * before performing any tasks. This method is overridden in the
	 * derived class to perform the following operations:
	 * 
	 * <ol>
	 * 
	 * <li>It first initializes the list of known hosts (servers
	 * we have connected to before) from the ".KnownHosts" file.</li>
	 * 
	 * <li>It ensures that the host name for the remote host can 
	 * be resolved, thereby establishing its basic validity.</li>
	 * 
	 * <li>It connects to the remote server by creating a 
	 * Ganymede connection object.</li>
	 * 
	 * <li>It then authenticates with the user by providing the
	 * user name and password stored in the Server object supplied
	 * when this class was insantiated.</li>
	 * 
	 * </ol>
	 * 
	 * @note The process of establishing a connection can be a 
	 * time consuming task. In some cases incorrect host names can
	 * cause long connection times (until the connection times out
	 * which can be in minutes). Consequently, it is best to call
	 * this method from a separate <b>daemon</b> thread.
	 * 
	 * @throws IOException This method throws IO exceptions in 
	 * the case of errors. If an error occurs, then a connection
	 * was not established and the caller will have to try again
	 * to establish a connection.
	 */
	@Override
	public void connect() throws IOException {
		// Load known hosts information first.
		loadKnownHosts();
		// Check and ensure that the host name is valid. The following
		// call will generate an exception if host is invalid.
		InetAddress.getByName(server.getName());
		// Create connection to the remote host. Use a temporary 
		// variable so that if exceptions are thrown our instance 
		// variable continues to remain valid.
		Connection connection = new Connection(server.getName());
		// Setup preferred key algorithms to be used for this host if the
		// host is in the known hosts list.
		String[] hostkeyAlgos = 
			knownHosts.getPreferredServerHostkeyAlgorithmOrder(server.getName());
		if (hostkeyAlgos != null) {
			connection.setServerHostKeyAlgorithms(hostkeyAlgos);
		}
		// Establish the connection. The following call will call the
		// verifyServerHostKey to check if the host key is present in the
		// known host database.
		connection.connect(this);
		// Next authenticate the user using the supplied credentials
		if (!connection.isAuthMethodAvailable(server.getUserID(), "password")) {
			throw new IOException("Non-interactive authentication using " +
					"password is not supported by remote server.\n" +
					"That is rather strange. You need to contact the " +
					"system adminstrator to enable support for\n" +
			"non-interactive logins using password.");
		}
		int retryCount = 3;
		do {
			// If password is null get a new password.
			getPassword();
			// OK, try to connect to the remote server by authenticating
			// using the user's ID and password
			if (!connection.authenticateWithPassword(server.getUserID(), 
					server.getPassword())) {
				// Authentication failed. Either the username or password
				// is incorrect. Generate exception.
				if (retryCount == 0) {
					throw new IOException("Authentication failed.\n" +
						"The user name or password is incorrect.");
				}
				retryCount--;
				server.setPassword(null); // reset password
			} else {
				break;
			}
		}  while (retryCount > 0);
		// OK, the connection was established successfully. Set our
		// instance variable to the local value.
		this.connection = connection;
	}

	/**
	 * Helper method to check and get password from the user.
	 */
	private void getPassword() throws IOException {
		if (server.getPassword() != null) {
			// We have valid password to work with.
			return;
		}
		// Create text fields for user name and password.
		JTextField userID = new JTextField(10);
		userID.setEditable(false);
		userID.setText(server.getUserID());
		JPasswordField password = new JPasswordField(10);

		// Create components by laying them out appropriately
		JPanel credPanel = new JPanel(new GridLayout(2, 2, 0, 3));
		credPanel.add(new JLabel("User id:"));
		credPanel.add(userID);
		Utilities.adjustDimension(userID, 0, 6);
		credPanel.add(new JLabel("Password:"));
		Utilities.adjustDimension(password, 0, 6);
		credPanel.add(password);
		// Another panel to control the size o the grid layout 
		// to ensure it looks decent
		JPanel msgPanel = new JPanel(new BorderLayout(0, 5));
		msgPanel.add(credPanel, BorderLayout.SOUTH);
		// Add a label indicating server information.
		JLabel subInfo = new JLabel("<html>Enter login credentials for <b>" + 
				server.getName() + "</b></html>");
		msgPanel.add(subInfo, BorderLayout.NORTH);
		// If purpose has been given add purpose information into another
		// panel.
		if (purpose != null) {
			JLabel info = new JLabel(purpose, Utilities.getIcon("images/32x32/Information.png"), 
									 JLabel.LEFT);
			 JPanel outer = new JPanel(new BorderLayout(0, 0));
			 outer.add(info, BorderLayout.CENTER);
			 outer.add(msgPanel, BorderLayout.SOUTH);
			 // Set message panel to be the outer most one now.
			 msgPanel = outer;
		}
		// Pack all the elements into an array
		Object items[] = { msgPanel };
		int result = 
		JOptionPane.showConfirmDialog(null, items, "Enter Password", 
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);
		if (result == JOptionPane.OK_OPTION) {
			server.setUserID(userID.getText());
			server.setPassword(new String(password.getPassword()));
		} else {
			throw new IOException("User cancelled credential input");
		}
	}
	
	/**
	 * Method to disconnect from a remote server. 
	 * 
	 * This method disconnects from the remote server if it is 
	 * connected. All current sessions will be terminated.
	 */
	@Override
	public void disconnect() {
		if (connection != null) {
			connection.close();
			connection = null;
			osType = null;
		}
	}

	/**
	 * This method can be used to run a brief command that 
	 * produces succinct output.
	 * 
	 * This method provides a convenient API to execute a 
	 * command on a remote server, buffer the output, and
	 * return the resulting output as a string. Since this
	 * method buffers and returns the output, it must be used
	 * only for jobs that return small volumes of data. In 
	 * addition, since the method returns the results only
	 * after the command has fully executed, it does not
	 * provide any interactive features to the user. 
	 * Consequently, it should be used only for jobs that
	 * run for a short period of time.
	 * 
	 * @note The connection to the remote server must have
	 * been established successfully via a call to connect method.
	 * 
	 * @param command The command line to be executed. This
	 * command must be compatible with the target machine's 
	 * OS. Otherwise this method will generate an exception.
	 * 
	 * @param outputs The buffered results from the standard
	 * output and standard error streams of the remote process.
	 * 
	 * @return The exit code from the remote command that was
	 * run.
	 * 
	 * @throws Exception If the execution produces an error then
	 * this method throws an exception.
	 */
	@Override
	public int exec(final String command, String[] outputs) throws Exception {
		// The following code was developed based on the API example
		// from Ganymede SSH.
		if (connection == null) {
			throw new IOException("Not connected to remote server.");
		}
		// Create a session via the connection
		Session session = connection.openSession();
		// Setup the command for execution on remote machine.
		session.execCommand(command);
		// Process the output streams. Stream gobbler is used to 
		// buffer the data from standard out (on a different thread)
		// while we read standard error in this thread.
		StreamGobbler stdout = new StreamGobbler(session.getStdout());
		// Read all the data into a single string from standard error
		// stream into a single string.
		outputs[1] = Utilities.readFullStream(session.getStderr());
		// Convert the standard output to a string as well.
		outputs[0] = Utilities.readFullStream(stdout);
		// Save exit status.
		Integer exitStatus = session.getExitStatus();
		session.close();
		// Some times the exit status can be null! Don't know why but it is
		// Have to file it with Ganymede SSH implementation. If exitStatus is
		// not null then use it. Otherwise use the length of standard error
		// as a measure of the exit status.
		int exitCode = (exitStatus != null) ? exitStatus : 
			outputs[1].length();
		// Return the exit code from the remote process
		return exitCode;
	}

	/**
	 * This method can be used to run a long running command that 
	 * may produce verbose output.
	 * 
	 * This method provides a convenient API to execute a 
	 * command on a remote server, stream the output, and
	 * return the exit code from the remote command. The
	 * outputs are streamed to a given StyledDocument. This
	 * method uses styles named "stdout", "stderr", and 
	 * "warning" (if available in the given output document)
	 * 
	 * @note The connection to the remote server must have
	 * been established successfully via a call to connect method.
	 * 
	 * @param command The command line to be executed. This
	 * command must be compatible with the target machine's 
	 * OS. Otherwise this method will generate an exception.
	 * 
	 * @param output The styled document to which the standard
	 * output and standard error streams are to be written.
	 * 
	 * @return The exit code from the remote command that was
	 * run.
	 * 
	 * @throws Exception If the execution produces an error then
	 * this method throws an exception.
	 */
	@Override
	public int exec(final String command, DefaultStyledDocument output) throws Exception {
		// Ensure we have a valid connection first.
		if (connection == null) {
			throw new IOException("Not connected to remote server.");
		}
		// Create a session via the connection
		Session session = connection.openSession();
		// Setup the command for execution on remote machine.
		session.execCommand(command);
		// Process the output streams. 
		StreamGobbler stderr = new StreamGobbler(session.getStderr());
		// Read stdout one line at a time and add it to the output
		BufferedReader stdin = new BufferedReader(new InputStreamReader(session.getStdout()));
		String line = null;
		while ((line = stdin.readLine()) != null) {
			// Got another line of standard output. Display it.
			output.insertString(output.getLength(), line + "\n", 
					output.getStyle("stdout"));
		}
		// Flush out any pending data on the standard error stream
		String stdErrData = Utilities.readFullStream(stderr);
		output.insertString(output.getLength(), stdErrData, 
				output.getStyle("stderr"));
		// Save exit status.
		Integer exitStatus = session.getExitStatus();
		session.close();
		// Some times the exit status can be null! Don't know why but it is
		// Have to file it with Ganymede SSH implementation. If exitStatus is
		// not null then use it. Otherwise use the length of standard error
		// as a measure of the exit status.
		int exitCode = (exitStatus != null) ? exitStatus : stdErrData.length();
		// Return the exit code from the remote process
		return exitCode;
	}

	/**
	 * Interactive verifier for use with Ganymede SSH callback.
	 * 
	 * This method is called by the Ganymede SSH layer once it has
	 * established initial communication with the remote server. This
	 * method is invoked to verify if the SSH client should proceed
	 * with the connection, given the server's credentials. This
	 * method checks the credentials of the server against the values
	 * in the KnownHosts file. If the value is present then this 
	 * method proceeds with the connection. Otherwise it provides
	 * the necessary information to the user to determine if the
	 * user wants to proceed. If the user wants to proceed then this
	 * method adds the host information to the KnownHosts file and
	 * returns true.
	 * 
	 * @param hostname The host name to be added to the list.
	 * @param port The port associated with the remote connection. 
	 * Currently this is not used.
	 * @param serverHostKeyAlgorithm The string name for the certificate
	 * encryption algorithm (such as: ssh-rsa2 etc.) 
	 * @param serverHostKey The host key (actual digest/figerprint) 
	 * 
	 *  @return This method returns true to indicate that the 
	 *  connection should proceed further.
	 */
	public boolean verifyServerHostKey(String hostname, int port, 
			String serverHostKeyAlgorithm, byte[] serverHostKey) throws Exception {
		// Look up the knownHosts list to see if this host's keys match up with
		// the values that we have seen before.
		int result;
		synchronized (knownHostsLock) {
			result = knownHosts.verifyHostkey(hostname, 
					serverHostKeyAlgorithm, serverHostKey);
		}
		if (result == KnownHosts.HOSTKEY_IS_OK) {
			// This server is in the data base and checks out OK. 
			// Go ahead and complete handshake with the server.
			return true;
		}

		// This remote server is not know to us or its fingerprint has changed!
		// The typicall norm here is to let the user know about it and check to
		// ensure that the user want's to actually connect to the machine. Below
		// we build a message to be displayed to the user.
		String message;
		if (result == KnownHosts.HOSTKEY_IS_NEW) {
			// Possibly the first time connection but must verify with user.
			message = "This is the first time PEACE is connecting to:\n";
		} else {
			// The host key has changed. This warrants a warning.
			message = "The host key for the server has changed! <br>";
		}
		// Add additional information and finger prints to the message:
		message += "     Host: " + hostname;
		InetAddress serverAddress = InetAddress.getByName(hostname);
		message += " (IP: " + serverAddress.getHostAddress() + ")\n";
		// Obtain and add finger prints
		String hexFingerprint = KnownHosts.createHexFingerprint(serverHostKeyAlgorithm, serverHostKey);
		String babbleFingerprint = KnownHosts.createBubblebabbleFingerprint(serverHostKeyAlgorithm,
				serverHostKey);
		message += "     Server's " + serverHostKeyAlgorithm + " Fingerprint: " + 
		hexFingerprint + "\n     Server's Babble Fingerprint: " + 
		babbleFingerprint + "\n";
		// Add the trailing question.
		message += "Do you want to procceed with the connection?";

		// Verify if the user wants to proceed
		int choice = JOptionPane.showConfirmDialog(parent, message, 
				"Connect?", JOptionPane.YES_NO_OPTION);
		if (choice == JOptionPane.YES_OPTION) {
			// Use helper method to add known host
			addKnownHost(hostname, port, serverHostKeyAlgorithm, 
					serverHostKey);
			return true;
		}
		return false;
	}

	/**
	 * Method to add a new host entry to the list of known hosts.
	 * 
	 * This method provides an MT-safe approach to add a new remote
	 * server/host entry to the list of known hosts. This method adds
	 * a new entry to both the in-memory knownHosts list and to the
	 * persistent file containing the list of known hosts.
	 * 
	 * @note This method checks to ensure that the host being added
	 * is indeed unique. This method consumes any exceptions generated
	 * by Ganymede (but logs it in the programmer log).
	 * 
	 * @param hostname The host name to be added to the list.
	 * @param port The port associated with the remote connection. 
	 * Currently this is not used.
	 * @param serverHostKeyAlgorithm The string name for the certificate
	 * encryption algorithm (such as: ssh-rsa2 etc.) 
	 * @param serverHostKey The host key (actual digest/figerprint) 
	 * as reported by the remote host.
	 */
	private void addKnownHost(String hostname, int port, 
			String serverHostKeyAlgorithm, byte[] serverHostKey) {
		// Ganymede SSH example suggests a hashed host name.
		String hashedHostname = KnownHosts.createHashedHostname(hostname);
		synchronized (knownHostsLock) {
			try	{
				if (knownHosts.verifyHostkey(hostname, serverHostKeyAlgorithm, 
						serverHostKey) == KnownHosts.HOSTKEY_IS_OK) {
					// This is a duplicate entry. Ignore it.
					return;
				}
				// Add the host key to the in-memory database as indicated in
				// Ganymede example.
				knownHosts.addHostkey(new String[] { hashedHostname }, 
						serverHostKeyAlgorithm, serverHostKey);
				// Also try to add the key to a known host file
				KnownHosts.addHostkeyToFile(new File(KnownHostPath), 
						new String[] { hashedHostname },
						serverHostKeyAlgorithm, serverHostKey);
			} catch (IOException exp) {
				// Log error in programmer log
				ProgrammerLog.log(exp);
			}
		}
	}

	/**
	 * This is a helper method that is invoked just before a remote session
	 * establishes connection with a server. This method attempts to load
	 * the list of known servers from the "KnownHosts" file. This file is
	 * stored in the main PEACE folder in user's home directory.
	 * 
	 * @note This method loads the known hosts file only once, the first
	 * time it is called in the GUI process. Subsequent calls to this method
	 * simply return. Therefore calling this method frequently is OK. 
	 */
	private void loadKnownHosts() {
		synchronized(knownHostsLock) {
			if (knownHosts != null) {
				// We have an known hosts database. Nothing else to be done.
				return;
			}
			// Create the known hosts class for further use
			knownHosts = new KnownHosts();
			// Try and load the known hosts data base for use.
			File knownHostFile   = new File(KnownHostPath);
			if (knownHostFile.exists()) {
				try	{
					// Load the data from the known hosts file.
					knownHosts.addHostkeys(knownHostFile);
				} catch (IOException e)	{
					// Log error in programmer log.
					ProgrammerLog.log(e);
					UserLog.log(UserLog.LogLevel.WARNING, "SSH", 
					"Error loading data from Known Hosts file.");
				}
			} else {
				// Try and create a default empty file so that we can
				// append to it at a later date.
				try {
					knownHostFile.createNewFile();
					// Cut a user log.
					UserLog.log(UserLog.LogLevel.NOTICE, "SSH",
					"Created initial empty known hosts file.");
				} catch (IOException e) {
					ProgrammerLog.log(e);
					UserLog.log(UserLog.LogLevel.WARNING, "SSH", 
					"Error creating initial empty Known Hosts file.");
				}
			}
		}
	}

	/**
	 * Determine the type of OS that this session is connected to.
	 * 
	 * This method can be used to determine the type of the OS that
	 * this session is connected to. 
	 * 
	 * @note The session must be connected prior to this call.
	 * 
	 * @return The type of OS this session is associated with.
	 * 
	 * @throws Exception If the detection of OS fails or a connection
	 * does not exist then this method throws an exception.
	 */
	@Override
	public OSType getOSType() throws Exception {
		if (connection == null) {
			throw new IOException("Remote session not connected.");
		}
		if (osType != null) {
			// We already know the OS type os just return it.
			return osType;
		}

		// Now that the session is connected try to run uname -a to
		// determine type of the remote operating system.
		String streamsData[] = {"", ""};
		if ((exec("uname", streamsData) != 0) ||
				(streamsData[0].length() == 0) ||
				(streamsData[1].length() != 0)) {
			throw new Exception("Unable to determine the remote " + 
					"machine's OS type. Possibly it is not a Linux " +
			"or an Unix machine.");
		}
		// Determine OS type based on the response string
		osType = (streamsData[0].indexOf("Linux") != -1) ? 
				OSType.LINUX : OSType.UNIX;
		return osType;
	}

	/**
	 * Copy given data from an input stream to a given file on the
	 * remote machine.
	 *
	 * @note This method uses SFTP to copy the data. The connection
	 * to the remote host must have already been established before
	 * invoking this method.
	 * 
	 * @param srcData The source stream that provides the data to
	 * be copied.
	 * 
	 * @param destDirectory The destination directory to which the
	 * data is to be copied. This method assumes that the remote
	 * directory has already been created.
	 * 
	 * @param destFileName The name of the destination file to which
	 * the data is to be copied.
	 * 
	 * @param mode The POSIX compliant mode string (such as: "0600"
	 * or "0777" to be used as the mode for the target file.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	public void copy(InputStream srcData, String destDirectory, 
			String destFileName, String mode) throws IOException {
		if (connection == null) {
			throw new IOException("Remote session not connected.");
		}
		// Create file attributes for passing to server
		SFTPv3FileAttributes attr = new SFTPv3FileAttributes();
		attr.permissions = Integer.parseInt(mode, 8);
		// Compute the final path and target file name.
		if (destDirectory == null) {
			destDirectory = ".";
		}
		final String targetFile = destDirectory + "/" + destFileName; 
		// Create an SFTP client.
		SFTPv3Client sftp = null;
		SFTPv3FileHandle destFile = null;
		// A try..finally block to ensure sftp and ftp connections
		// get closed
		try {
			sftp = new SFTPv3Client(connection);
			// Create a file for writing (truncate any existing files)
			destFile = sftp.createFileTruncate(targetFile, attr);
			// Read data from input stream and write data to the destFile
			byte buffer[]  = new byte[8096];
			int  bytesRead = 0; long destOffset = 0;
			while ((bytesRead = srcData.read(buffer, 0, buffer.length)) != -1) {
				sftp.write(destFile, destOffset, buffer, 0, bytesRead);
				destOffset += bytesRead;
			}
		} finally {
			if (destFile != null) {
				sftp.closeFile(destFile);
			}
			if (sftp != null) {
				sftp.close();
			}
		}
	}

	/**
	 * Copy file from a remote machine to a given output stream.
	 *
	 * @note This method uses SFTP to copy the data. The connection
	 * to the remote host must have already been established before
	 * invoking this method.
	 * 
	 * @param destData The destination stream to which the data is
	 * to be written.
	 * 
	 * @param srcDirectory The source directory from where the
	 * file is to be copied. 
	 * 
	 * @param srcFileName The name of the source file from where
	 * the data is to be copied.
	 * 
	 * @param progBar The progress bar to be used indicate the
	 * file copy progress.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	@Override
	public void copy(OutputStream destData, String srcDirectory, 
			String srcFileName, JProgressBar progBar) throws IOException {
		if (connection == null) {
			throw new IOException("Remote session not connected.");
		}
		// Compute the final path and target file name.
		if (srcDirectory == null) {
			srcDirectory = ".";
		}
		final String sourceFile = srcDirectory + "/" + srcFileName; 
		// Create an SFTP client.
		SFTPv3Client sftp        = null;
		SFTPv3FileHandle srcFile = null;
		// A try..finally block to ensure sftp and ftp connections
		// get closed
		try {
			// Make the progress bar indeterminate for now.
			if (progBar != null) {
				progBar.setIndeterminate(true);
			}
			sftp = new SFTPv3Client(connection);
			// Create a file for reading (or get exception)
			srcFile = sftp.openFileRO(sourceFile);
			// Now determine the file size and set progress bar size
			if (progBar != null) {
				SFTPv3FileAttributes attribs = sftp.lstat(sourceFile);
				progBar.setMaximum((int) attribs.size.intValue());
				progBar.setIndeterminate(false);
			}
			// Read data from input stream and write data to the destFile
			byte buffer[]  = new byte[4096];
			int  bytesRead = 0; 
			long srcOffset = 0;
			while ((bytesRead = sftp.read(srcFile, srcOffset, buffer, 0, buffer.length)) != -1) {
				destData.write(buffer, 0, bytesRead);
				srcOffset += bytesRead;
				if (progBar != null) {
					progBar.setValue((int) srcOffset);
				}
			}
		} finally {
			// Ensure that the file is closed correctly
			if (srcFile != null) {
				sftp.closeFile(srcFile);
			}
			if (sftp != null) {
				sftp.close();
			}
		}
	}
	
	/**
	 * Obtain information about a given path on the remote machine.
	 * 
	 * @note This method uses SFTP to copy the data. The connection
	 * to the remote host must have already been established before
	 * invoking this method.
	 * 
	 * @param path The path (absolute or relative to home directory)
	 * of the file whose meta data is to be retrieved.
	 * 
	 * @return A FileInfo object containing the information about
	 * the path.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	public FileInfo fstat(String path) throws IOException {
		if (connection == null) {
			throw new IOException("Remote session not connected.");
		}
		// Create an SFTP client.
		SFTPv3Client sftp = null;
		FileInfo info = null;
		try {
			sftp = new SFTPv3Client(connection);
			SFTPv3FileAttributes attribs = sftp.lstat(path);
			if (attribs != null) {
				// Translate SFTP file attributes to FileInfo style
				// attributes.
				int attributes = 0;
				attributes |= (attribs.isDirectory()   ? FileInfo.DIR_ATTRIB  : 0);
				attributes |= (attribs.isRegularFile() ? FileInfo.FILE_ATTRIB : 0);
				attributes |= ((attribs.uid & 00400) != 0) ? FileInfo.READ_ATTRIB : 0;
				attributes |= ((attribs.uid & 00200) != 0) ? FileInfo.WRITE_ATTRIB : 0;
				attributes |= ((attribs.uid & 00100) != 0) ? FileInfo.EXEC_ATTRIB : 0;
				// Now create the file info object.
				info = new FileInfo(path, attribs.atime.longValue() * 1000L, 
						attribs.size.longValue(), attributes);
			} else {
				// The file does not exist? Create a dummy info.
				info = new FileInfo(path, -1, -1, 0);
			}
		} catch (SFTPException exp) {
			if (exp.getMessage().startsWith("No such file")) {
				// This is an expected exception if file does not exist
				info = new FileInfo(path, -1, -1, 0);
			} else {
				// Unexcepted exception
				throw exp;
			}
		} finally {
			if (sftp != null) {
				sftp.close();
			}
		}
		// Return file information to caller.
		return info;
	}

	/**
	 * Creates a directory on the server.
	 * 
	 * This method must be used to create a directory entry on the
	 * server.
	 * 
	 * @note The connection to the server must have already been 
	 * successfully established via a call to the connect method 
	 * before invoking this method.
	 * 
	 * @param directory The fully path to the directory to be created.
	 * 
	 * @throws IOException This method throws an exception if the 
	 * directory could not be created.
	 */
	public void mkdir(String directory) throws Exception {
		if (connection == null) {
			throw new IOException("Remote session not connected.");
		}
		// Now that the session is connected try to run mkdir to
		// create the specified directory on the remote machine.
		SFTPv3Client sftp = null;
		try {
			sftp = new SFTPv3Client(connection);
			sftp.mkdir(directory, 0700);
		} finally {
			if (sftp != null) {
				sftp.close();
			}
		}
	}

	/**
	 * Deletes an <i>empty</i> directory on the server.
	 * 
	 * This method must be used to delete a directory entry on the
	 * server. The directory must be empty.
	 * 
	 * @note The connection to the server must have already been 
	 * successfully established via a call to the connect method 
	 * before invoking this method.
	 * 
	 * @param directory The fully path to the directory to be deleted.
	 * 
	 * @throws IOException This method throws an exception if the 
	 * directory could not be deleted.
	 */
	public void rmdir(String directory) throws Exception {
		if (connection == null) {
			throw new IOException("Remote session not connected.");
		}
		// Now that the session is connected try to run rmdir to
		// delete the specified directory on the remote machine
		SFTPv3Client sftp = null;
		try {
			sftp = new SFTPv3Client(connection);
			sftp.rmdir(directory);
		} finally {
			if (sftp != null) {
				sftp.close();
			}
		}
	}
	
	/**
	 * A simple method to set a purpose message for this session.
	 * 
	 * This method can be used to set up a purpose message for a
	 * server session. The purpose message is displayed
	 * to the user when prompting for inputs from the user for 
	 * credentials. The message serves the purpose of appraising the
	 * user about the purpose of the session.
	 *  
	 * @param text This string is used to create a label (possibly
	 * with an icon on it). So it can be plain text or HTML. If the
	 * message is long, then ensure it is properly broken into 
	 * multiple lines so that dialog boxes don't get too large.
	 */
	public void setPurpose(String text) {
		purpose = text;
	}
	
	/**
	 * The connection to the remote server via which the remote server can
	 * be accessed for performing various operations. The connection is
	 * created via the connect()  method.
	 */
	private Connection connection;

	/**
	 * The last known OS type for the remote server. This value is
	 * set after the getOSType() method is called. If a connection
	 * is lost or closed, then this value is reset.
	 */
	private OSType osType;

	/**
	 * This is a convenience class that is provided by Ganymede SSH to 
	 * store information about servers/hosts that we have already
	 * connected to. This data is loaded once initially and is shared
	 * by all RemoteServerSession instances. This enables us to load the
	 * data once and keep updating and using it rather than loading it
	 * each time (which should save some CPU cycles).
	 */
	private static KnownHosts knownHosts = null;

	/**
	 * This object is merely used to arbitrate access to the knownHosts
	 * object so that operations on knownHosts is MT-Safe.
	 */
	private static Object knownHostsLock = new Boolean(false);

	/**
	 * A simple textual information indicating the purpose for this
	 * session. This string is more meanigful to the user and is 
	 * merely used to provide additional information when prompting
	 * for inputs from the user.
	 */
	private String purpose;
	
	/**
	 * The OS-specific path where the list of known hosts (that is
	 * the servers to which we have connected before) is stored. This
	 * list is applicable only for remote hosts.
	 */
	private static final String KnownHostPath = 
		Utilities.getDefaultDirectory() + "/.KnownHosts";
}
