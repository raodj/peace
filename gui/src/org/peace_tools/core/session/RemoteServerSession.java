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

import java.awt.Component;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.InetAddress;
import java.net.ServerSocket;

import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.FileInfo;
import org.peace_tools.generic.Log.LogLevel;
import org.peace_tools.generic.PasswordDialog;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

import com.jcraft.jsch.ChannelExec;
import com.jcraft.jsch.ChannelSftp;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.SftpATTRS;
import com.jcraft.jsch.SftpException;
import com.jcraft.jsch.UIKeyboardInteractive;
import com.jcraft.jsch.UserInfo;

/**
 * A remote server session based on the secure shell (SSH) protocol.
 * 
 * <p>
 * This class provides an implementation of the ServerSession API. Specifically,
 * this class provides a session that can be used to interact with a remote host
 * via the secure shell (SSH) protocol. The secure shell protocol is a de facto
 * standard for interacting with remote servers via the Internet today. It
 * provides all the necessary security features to safely interact with remote
 * hosts and almost all super computing clusters mandate the use of SSH for
 * interactions.
 * </p>
 * 
 * <p>
 * This class uses the JSch SSH implementation for establishing SSH connections.
 * JSch SSH supports the recent SSH2 protocol. Please refer to JSch website for
 * further details: <a href="http://www.jcraft.com/jsch">JSch</a>
 * </p>
 * 
 * <p>
 * <b>Note:</b>The remote server session automatically turns-on maximum data
 * compression automatically.
 * </p>
 * 
 * <p>
 * PEACE redistributes JSch and includes their BSD style license file as per
 * JSch licensing requirements. Here is a link to JSch license file:
 * <a href="http://www.jcraft.com/jsch/LICENSE.txt">License</a>
 * </p>
 */
public class RemoteServerSession extends ServerSession implements UserInfo, UIKeyboardInteractive {
	/**
	 * The constructor merely passes parameters to the base class that initializes
	 * the instance variables. No special operations are performed in the
	 * constructor.
	 * 
	 * @param server The server data entry that provides the necessary information
	 *               to connect to the server.
	 *
	 * @param parent The parent component that should be used to create GUI elements
	 *               that may be needed for any interactive operations.
	 */
	protected RemoteServerSession(Component parent, Server server) {
		super(server, parent);
	}

	/**
	 * Connect to the server in order to perform various operations.
	 * 
	 * This method must be used to establish a connection to a server before
	 * performing any tasks. This method is overridden in the derived class to
	 * perform the following operations:
	 * 
	 * <ol>
	 * <li>If a connection has already been established then this method returns
	 * immediately.</li>
	 * 
	 * <li>It first initializes the list of known hosts (servers we have connected
	 * to before) from the ".KnownHosts" file.</li>
	 * 
	 * <li>It ensures that the host name for the remote host can be resolved,
	 * thereby establishing its basic validity.</li>
	 * 
	 * <li>It connects to the remote server by creating a JSch session object. The
	 * process of creating a session requires authentication by providing the user
	 * name and password stored in the Server object supplied when this class was
	 * instantiated. If a valid password is not available, then this method prompts
	 * the user to obtain the password. It is important to set the purpose for this
	 * session to inform the user about the need to connect.</li>
	 * </ol>
	 * 
	 * <p>
	 * <b>Note:</b> The process of establishing a connection can be a time consuming
	 * task. In some cases incorrect host names can cause long connection times
	 * (until the connection times out which can be in minutes).
	 * </p>
	 * 
	 * @throws IOException This method throws IO exceptions in the case of errors.
	 *                     If an error occurs, then a connection was not established
	 *                     and the caller will have to try again to establish a
	 *                     connection.
	 */
	@Override
	public void connect() throws IOException {
		// If connection is already set, then do nothing.
		if (session != null) {
			return;
		}

		// Setup the known hosts file as needed
		setKnownHosts();

		// Check and ensure that the host name is valid. The following call will
		// generate an exception if host is invalid.
		InetAddress.getByName(server.getName());

		// Create connection to the remote host. Use a temporary variable so that if
		// exceptions are thrown our instance variable continues to remain valid.
		Session connection = null;

		int retryCount = 3;
		do {
			try {
				// Try and establish a connection. If an exception occurs then we can't do
				// anything further in this step.
				connection = sshClientLibrary.getSession(server.getUserID(), server.getName(), server.getPort());
				connection.setUserInfo(this);
			} catch (JSchException e) {
				// This type of exception is not recoverable!
				throw new IOException(e);
			}
			try {
				// Request server to enable level-9 compression for this session.
				connection.setConfig("compression.s2c", "zlib@openssh.com,zlib,none");
				connection.setConfig("compression.c2s", "zlib@openssh.com,zlib,none");
				connection.setConfig("compression_level", "9");

				// Establish the connection. The following call will prompt for adding entries
				// to ".KnownHosts" file and for the user's password (if the server does not
				// have a password already set).
				connection.connect();

				// When control drops here and the connection is not connected, that indicates a
				// routine flow in which the user decided to cancel the connection
				if (!connection.isConnected()) {
					throw new IOException("User interrupted loggin to server");
				}
			} catch (JSchException e) {
				// UserLog.log(LogLevel.NOTICE, "RemoteServerSession", exp.getMessage());
				ProgrammerLog.log(e);

				// An exception occurred. If we are not out of retries, then we will prompt the
				// user for the password again. Reset password.
				server.setPassword(null);
				connection.disconnect();
			} catch (RuntimeException e) {
				if (USER_INTERRUPTED_EXP_MSG.equals(e.getMessage())) {
					// This is an exception that we expect to get
					throw new IOException(e);
				}
				// This is not an exception we expect to get. Rethrow it.
				throw e;
			} finally {
				// Always decrement retry count to ensure we never get into an infinite loop.
				retryCount--;
			}
		} while (retryCount > 0 && !connection.isConnected());

		if (retryCount > 0 && connection.isConnected()) {
			// The connection was established successfully. Set our instance variable to
			// the local value.
			this.session = connection;
		} else {
			// We could not connect even after 3-retries. Throw an exception and bail out
			// from here.
			throw new IOException("Authentication failed.\nThe user name or password is incorrect.");
		}
	}

	/**
	 * A helper method to setup the known hosts file (and information).
	 *
	 * This is a refactored helper method that is used to setup the known hosts file
	 * for use by JSch (SSH client library). The known hosts file caches SSH
	 * key/signatures of servers with which the user has established connections in
	 * previous runs. This file is automatically managed by JSch. Here, we just
	 * provide the file name for use (if we have not already done so). This method
	 * is invoked from the {@link #connect()} method.
	 */
	private void setKnownHosts() {
		// Setup the known hosts repository if one is not already set
		synchronized (sshClientLibrary) {
			if (sshClientLibrary.getHostKeyRepository().getKnownHostsRepositoryID() == null) {
				try {
					sshClientLibrary.setKnownHosts(KNOWN_HOST_PATH);
				} catch (JSchException exp) {
					// Let user know that known hosts saving is disabled!
					ProgrammerLog.log(exp);
					JPanel msg = Utilities.collapsedMessage(KNOWN_HOSTS_ERROR, Utilities.toString(exp));
					JOptionPane.showMessageDialog(parent, msg, "Unable to use Known Hosts File",
							JOptionPane.WARNING_MESSAGE);
				}
			}
		}
	}

	/**
	 * Method to obtain the current secure-shell connection to the server.
	 * 
	 * This method can be used to obtain the current secure shell connection being
	 * used by this server session. Note that the {@link #connect()} method must be
	 * invoked to establish a valid session. If a valid connection (AKA JSch
	 * session) is returned by this method then it can be readily used to create a
	 * SSH-channel.
	 * 
	 * @return The current secure-shell connection to the server. If a connection
	 *         has not been established then this method returns null.
	 */
	public Session getConnection() {
		return session;
	}

	/**
	 * Implementation for corresponding method in JSch's UserInfo class.
	 * 
	 * This method provides a simple implementation for the corresponding method in
	 * the JSch's UserInfo class. This method checks to see if the server (to which
	 * the session is being established) already has a password setup. If a valid
	 * password is already available this method simply returns that password value.
	 * Otherwise this method calls the promptPassword() method to obtain a password
	 * from the user.
	 * 
	 * @return This method returns the current password set for the server. If a
	 *         valid password is not set then this method returns null.
	 */
	@Override
	public String getPassword() {
		if (server.getPassword() == null) {
			promptPassword(null);
		}
		return server.getPassword();
	}

	/**
	 * Implementation for JSch's UserInfo interface to check and get password from
	 * the user.
	 * 
	 * This method implements the corresponding method in JSch's UnserInfo
	 * interface. This method is invoked from within JSch to prompt the user for a
	 * password if an SSH handshake requires a password to be entered. This method
	 * prompts the user for a password only if {@code server.getPassword()} returns
	 * null. This method displays a properly layed-out modal dialog box to obtain
	 * the password. The purpose (if set) is displayed when prompting for password.
	 *
	 * @param message The message to be displayed to the user. This message is
	 *                passed in by JSch. Currently, this value is ignored by this
	 *                method (and therefore it can be null) giving higher preference
	 *                to the purpose set for the session.
	 * 
	 * @return This method returns true if the user clicks the "OK" button to
	 *         proceed with SSH handshake. If the user clicks the "Cancel" button,
	 *         then this method returns false indicating that further SSH processing
	 *         should not occur.
	 * 
	 * @see setPurpose
	 */
	@Override
	public boolean promptPassword(String message) {
		if (server.getPassword() != null) {
			// We have valid password to work with.
			return true;
		}
		PasswordDialog pwdDialog = new PasswordDialog(server.getUserID(), false, server.getName(), purpose);
		if (pwdDialog.prompt(parent)) {
			server.setUserID(pwdDialog.getUserID());
			server.setPassword(pwdDialog.getPassword());
			return true;
		}
		// When control drops here that indicates that the user did not click the "OK"
		// button. So stop SSH operations
		return false;
	}

	/**
	 * Dummy implementation for corresponding method in JSch's UserInfo class.
	 * 
	 * This method is expected to return a pass phrase to be used with SSH-digests.
	 * Currently, PEACE does not support this feature and consequently, this method
	 * merely provides a dummy implementation.
	 * 
	 * @return Currently, this method always returns null as SSH digests are not
	 *         fully implemented.
	 */
	@Override
	public String getPassphrase() {
		return null;
	}

	/**
	 * Dummy implementation for corresponding method in JSch's UserInfo class.
	 * 
	 * This method is expected to prompt the user for a pass phrase to be used with
	 * SSH-digests. Currently, PEACE does not support this feature and consequently,
	 * this method merely provides a dummy implementation.
	 * 
	 * @return Currently, this method always returns false indicating that further
	 *         SSH handshake should not occur using pass phrases
	 */
	@Override
	public boolean promptPassphrase(String message) {
		return false;
	}

	/**
	 * Simple implementation for corresponding method in JSch's UserInfo class.
	 * 
	 * This method is expected to display messages to the user. The messages are
	 * typically initial license or usage requirements that are sent by the server
	 * as a part of the SSH-protocol hand shake. This method merely displays the
	 * information to the user.
	 * 
	 * @param message The message to be displayed to the user.
	 */
	@Override
	public void showMessage(String message) {
		// Create a text area to display the message and place it inside a scroll pane
		// to permit display of large messages using a decent sized GUI window.
		JTextArea msgDisplay = new JTextArea(message, 60, 5);
		JScrollPane jsp = new JScrollPane(msgDisplay);

		// Show the message to the user.
		JOptionPane.showMessageDialog(parent, jsp, "Secure Shell (SSH) Message", JOptionPane.INFORMATION_MESSAGE);
	}

	/**
	 * Implementation for corresponding method in JSch's UserInfo class.
	 * 
	 * This method provides a simple implementation for the corresponding method in
	 * the JSch's UserInfo class. This method is typically invoked by JSch to check
	 * if a user would like to connect to a server whose SSH-signature is not in the
	 * ".KnownHosts" file. Essentially this method displays a dialog (with suitable
	 * message) and prompts the user for a Yes/No type answer.
	 * 
	 * @param message The message to be displayed to the user. This message is
	 *                provided by JSch and we don't have much control on it.
	 * 
	 * @return This method returns true if the user clicks on the "Yes" button. If
	 *         the user click the "No" button then this method returns false.
	 */
	@Override
	public boolean promptYesNo(String message) {
		// The message object to be filled in further below.
		JComponent msgDisplay = null;
		int msgType = JOptionPane.INFORMATION_MESSAGE;
		boolean expOnNo = false;

		// JSch messages are OK for folks who are familiar with SSH and such. However,
		// for a common user, here try and create a more meaningful message in
		// situations where the user is connecting to the server for the first time.
		if (message.startsWith("The authenticity of host") && message.contains("can't be established")) {
			// This is a situation where the user is typically connecting for the first
			// time. Display a custom message.
			String rsaTag = "RSA key fingerprint is ";
			int rsaKeyPos = message.indexOf(rsaTag) + rsaTag.length() + 1;
			String rsaKey = message.substring(rsaKeyPos, message.indexOf('.', rsaKeyPos));

			// Format and display the message to the user.
			String custMsg = String.format(FIRST_SSH_CONNECTION_MSG, server.getName(), rsaKey);
			msgDisplay = Utilities.collapsedMessage(custMsg, message, false);
			expOnNo = true;
		} else if (message.startsWith("WARNING: ")) {
			// This is a situation when the RSA key for the server has changed. Here we
			// create a custom and more informative message than the one displayed by JSch.
			msgDisplay = Utilities.collapsedMessage(SSH_HOST_CHANGE_MSG, message, false);
			msgType = JOptionPane.WARNING_MESSAGE;
			expOnNo = true;
		} else {
			// Create a text area to display the message from JSch in other cases
			JTextArea jta = new JTextArea(message, 60, 5);
			msgDisplay = new JScrollPane(jta);
		}
		// Display message to user and get user's Yes/No choice
		int choice = JOptionPane.showConfirmDialog(parent, msgDisplay, "Secure Shell (SSH) Protocol Interaction",
				JOptionPane.YES_NO_OPTION, msgType);

		// If the user clicked "No" for a "Warning" message, we throw an exception here
		// to communicate a serious issue.
		if (choice == JOptionPane.NO_OPTION && expOnNo) {
			// A serious issue from which we should not retry connections etc.
			throw new RuntimeException(USER_INTERRUPTED_EXP_MSG);
		}

		// Return true if the user clicked on the "Yes" button.
		return (choice == JOptionPane.YES_OPTION);
	}

	/**
	 * Implementation for corresponding method in JSch's UIKeyboardInteractive
	 * class.
	 * 
	 * This method is typically invoked by JSch to get user's input when additional
	 * information is needed to authenticate, e.g 2nd factor authentication. It
	 * displays an input dialog and get the user's response.
	 *
	 * @return the user's response
	 */
	@Override
	public String[] promptKeyboardInteractive(String destination, String name, String instruction, String[] prompt,
			boolean[] echo) {
		String[] str = new String[1];
		if (prompt[0].startsWith("Password")) {
			str[0] = getPassword();
		} else {
			String result = JOptionPane.showInputDialog(parent, String.join("\n", prompt));
			str[0] = result;
		}
		return str;
	}

	/**
	 * Method to disconnect from a remote server.
	 * 
	 * This method disconnects from the remote server if it is connected. All
	 * current sessions will be terminated.
	 */
	@Override
	public void disconnect() {
		if (session != null) {
			session.disconnect();
			session = null;
			osType = null;
		}
	}

	/**
	 * This method can be used to run a brief command that produces succinct output.
	 * 
	 * This method provides a convenient API to execute a command on a remote
	 * server, buffer the output, and return the resulting output as a string. Since
	 * this method buffers and returns the output, it must be used only for jobs
	 * that return small volumes of data. In addition, since the method returns the
	 * results only after the command has fully executed, it does not provide any
	 * interactive features to the user. Consequently, it should be used only for
	 * jobs that run for a short period of time.
	 * 
	 * <p>
	 * <b>Note:</b> The connection to the remote server must have been established
	 * successfully via a call to connect method.
	 * </p>
	 * 
	 * @param command The command line to be executed. This command must be
	 *                compatible with the target machine's OS. Otherwise this method
	 *                will generate an exception.
	 * 
	 * @param outputs The buffered results from the standard output and standard
	 *                error streams of the remote process. Specifically, outputs[0]
	 *                will contain the standard output data while outputs[1] will
	 *                contain the standard error stream.
	 * 
	 * @return The exit code from the remote command that was run.
	 * 
	 * @throws Exception If the execution produces an error then this method throws
	 *                   an exception.
	 */
	@Override
	public int exec(final String command, String[] outputs) throws Exception {
		if (session == null) {
			throw new IOException("Not connected to remote server.");
		}

		// Create a channel via the session
		ChannelExec channel = (ChannelExec) session.openChannel("exec");

		// Setup the command for execution on remote machine.
		channel.setCommand(command);

		// Process the output streams. A suitable stream is used to buffer the data from
		// standard out while we read standard error explicitly.
		ByteArrayOutputStream stdout = new ByteArrayOutputStream(8192);
		// Setup the standard output stream to which data is to be written when the the
		// command runs.
		channel.setOutputStream(stdout);
		// No buffers for error stream.
		channel.setErrStream(null);

		// Now run the command on the remote server
		channel.connect();

		// Read all the data into a single string from standard error stream into a
		// single string.
		outputs[1] = Utilities.readFullStream(channel.getErrStream());
		// Convert the standard output to a string as well.
		outputs[0] = stdout.toString();

		// Save exit status.
		int exitCode = channel.getExitStatus();

		channel.disconnect();

		// Return the exit code from the remote process
		return exitCode;
	}

	/**
	 * This method can be used to run a long running command that may produce
	 * verbose output.
	 * 
	 * This method provides a convenient API to execute a command on a remote
	 * server, stream the output, and return the exit code from the remote command.
	 * The outputs are streamed to a given StyledDocument. This method uses styles
	 * named {@code "stdout"}, {@code "stderr"}, and {@code "warning"} (if available
	 * in the given output document)
	 * 
	 * <p>
	 * <b>Note:</b> The connection to the remote server must have been established
	 * successfully via a call to connect method.
	 * <p>
	 * 
	 * @param command The command line to be executed. This command must be
	 *                compatible with the target machine's OS. Otherwise this method
	 *                will generate an exception.
	 * 
	 * @param output  The styled document to which the standard output and standard
	 *                error streams are to be written.
	 * 
	 * @return The exit code from the remote command that was run.
	 * 
	 * @throws Exception If the execution produces an error then this method throws
	 *                   an exception.
	 */
	@Override
	public int exec(final String command, DefaultStyledDocument output) throws Exception {
		if (session == null) {
			throw new IOException("Not connected to remote server.");
		}

		// Create a session via the connection
		ChannelExec channel = (ChannelExec) session.openChannel("exec");

		// Setup the command for execution on remote machine.
		channel.setCommand(command);

		// Process the output streams. The following output stream buffers the data from
		// standard error (on a different thread) while we read standard output in this
		// thread.
		ByteArrayOutputStream stderr = new ByteArrayOutputStream(8192);
		channel.setErrStream(stderr);

		// Read stdout one line at a time and add it to the output
		BufferedReader stdin = new BufferedReader(new InputStreamReader(channel.getInputStream()));

		// Now run the command on the remote server
		channel.connect();

		String line = null;
		while ((line = stdin.readLine()) != null) {
			// Got another line of standard output. Display it.
			output.insertString(output.getLength(), line + "\n", output.getStyle("stdout"));
		}

		// Flush out any pending data on the standard error stream
		String stdErrData = stderr.toString();
		output.insertString(output.getLength(), stdErrData, output.getStyle("stderr"));

		// Save exit status.
		int exitCode = channel.getExitStatus();

		channel.disconnect();

		// Return the exit code from the remote process
		return exitCode;
	}

	/**
	 * Determine the type of OS that this session is connected to.
	 * 
	 * This method can be used to determine the type of the OS that this session is
	 * connected to.
	 * 
	 * <p>
	 * <b>Note:</b> The session must be connected prior to this call.
	 * </p>
	 * 
	 * @return The type of OS this session is associated with.
	 * 
	 * @throws Exception If the detection of OS fails or a connection does not exist
	 *                   then this method throws an exception.
	 */
	@Override
	public Server.OSType getOSType() throws Exception {
		if (osType != null) {
			return osType;
		}
		if (session == null) {
			throw new IOException("Remote session not connected.");
		}
		// Now that the session is connected try to run uname -a to determine type of
		// the remote operating system.
		String[] streamsData = new String[2];
		if (exec("uname", streamsData) != 0 || streamsData[0].length() == 0 || streamsData[1].length() != 0) {
			throw new Exception("Unable to determine the remote machine's OS type.\n"
					+ "Possibly it is not a Linux or an Unix machine.");
		}
		// Determine OS type based on the response string
		osType = (streamsData[0].indexOf("Linux") != -1) ? Server.OSType.LINUX : Server.OSType.UNIX;
		return osType;
	}

	/**
	 * Copy given data from an input stream to a given file on the remote machine.
	 *
	 * <p>
	 * <b>Note:</b> This method uses SFTP to copy the data. The connection to the
	 * remote host must have already been established before invoking this method.
	 * </p>
	 * 
	 * @param srcData       The source stream that provides the data to be copied.
	 * 
	 * @param destDirectory The destination directory to which the data is to be
	 *                      copied. This method assumes that the remote directory
	 *                      has already been created.
	 * 
	 * @param destFileName  The name of the destination file to which the data is to
	 *                      be copied.
	 * 
	 * @param mode          The POSIX compatible mode string (such as: 0600 or 0777)
	 *                      to be used as the mode for the target file. Note the
	 *                      leading zero! It is important to ensure that these
	 *                      numbers are Octal values.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	@Override
	public void copy(InputStream srcData, String destDirectory, String destFileName, int mode) throws IOException {
		if (session == null) {
			throw new IOException("Remote session not connected.");
		}

		// Compute the final path and target file name.
		if (destDirectory == null) {
			destDirectory = ".";
		}

		String targetFile = destDirectory + "/" + destFileName;

		// Create an SFTP client.
		ChannelSftp sftp = null;
		OutputStream destFile = null;

		// A try..finally block to ensure sftp and ftp connections get closed
		try {
			// Create a sftp channel and connect the channel.
			sftp = (ChannelSftp) session.openChannel("sftp");
			sftp.connect();

			// Create a file for writing (truncate any existing files)
			destFile = sftp.put(targetFile, ChannelSftp.OVERWRITE);

			// Read data from input stream and write data to the destFile
			byte[] buffer = new byte[8096];
			int bytesRead = 0;
			while ((bytesRead = srcData.read(buffer, 0, buffer.length)) != -1) {
				destFile.write(buffer, 0, bytesRead);
			}
			
			// Done writing data. Close stream and change permissions
            destFile.close();
            destFile = null;
			
            // Get current status information about the file.
			SftpATTRS attribs = sftp.stat(targetFile);
			attribs.setPERMISSIONS(mode);
			sftp.setStat(targetFile, attribs);
		} catch (JSchException | SftpException e) {
			throw new IOException(e);
		} finally {
			if (destFile != null) {
				destFile.close();
			}
			if (sftp != null) {
				sftp.disconnect();
			}
		}
	}

	/**
	 * Copy file from a remote machine to a given output stream.
	 *
	 * <p>
	 * <b>Note:</b> This method uses SFTP to copy the data. The connection to the
	 * remote host must have already been established before invoking this method.
	 * </p>
	 * 
	 * @param destData     The destination stream to which the data is to be
	 *                     written.
	 * 
	 * @param srcDirectory The source directory from where the file is to be copied.
	 * 
	 * @param srcFileName  The name of the source file from where the data is to be
	 *                     copied.
	 * 
	 * @param progBar      The progress bar to be used indicate the file copy
	 *                     progress.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	@Override
	public void copy(OutputStream destData, String srcDirectory, String srcFileName, JProgressBar progBar)
			throws IOException {
		if (session == null) {
			throw new IOException("Remote session not connected.");
		}

		// Compute the final path and target file name.
		if (srcDirectory == null) {
			srcDirectory = ".";
		}

		String sourceFile = srcDirectory + "/" + srcFileName;

		// Create an SFTP client.
		ChannelSftp sftp = null;
		InputStream srcFile = null;

		// A try..finally block to ensure sftp channels get closed
		try {
			// Make the progress bar indeterminate for now.
			if (progBar != null) {
				progBar.setIndeterminate(true);
			}

			// Open a sftp channel and connect to server to secure-FTP data
			sftp = (ChannelSftp) session.openChannel("sftp");
			sftp.connect();

			// Obtain information about the source file to be copied
			SftpATTRS attribs = sftp.stat(sourceFile);
			if (attribs == null) {
				throw new IOException("The source file " + sourceFile + " was not found.");
			}

			// Open source file for reading (or get exception)
			srcFile = sftp.get(sourceFile);

			// Now using the file size and set progress bar size
			if (progBar != null) {
				progBar.setMaximum((int) attribs.getSize());
				progBar.setIndeterminate(false);
			}

			// Read data from input stream and write data to the destFile
			byte[] buffer = new byte[4096];
			int bytesRead = 0;
			long srcOffset = 0;
			while ((bytesRead = srcFile.read(buffer, 0, buffer.length)) != -1) {
				destData.write(buffer, 0, bytesRead);
				srcOffset += bytesRead;
				if (progBar != null) {
					progBar.setValue((int) srcOffset);
				}
			}
		} catch (JSchException | SftpException e) {
			throw new IOException(e);
		} finally {
			// Ensure that the file is closed correctly
			if (srcFile != null) {
				srcFile.close();
			}
			if (sftp != null) {
				sftp.disconnect();
			}
		}
	}

	/**
	 * Obtain information about a given path on the remote machine.
	 * 
	 * <p>
	 * <b>Note:</b> This method uses SFTP to obtain the data. The connection to the
	 * remote host must have already been established before invoking this method.
	 * </p>
	 * 
	 * @param path The path (absolute or relative to home directory) of the file
	 *             whose meta data is to be retrieved.
	 * 
	 * @return A FileInfo object containing the information about the path.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	@Override
	public FileInfo fstat(String path) throws IOException {
		if (session == null) {
			throw new IOException("Remote session not connected.");
		}

		FileInfo info = null;

		try {
			RemoteFile rf = new RemoteFile(path, session);
			SftpATTRS attribs = rf.getAttributes();

			if (attribs != null) {
				String logEntry = String.format("SFTP gave 0x%s as attributes for %s\n",
						Integer.toHexString(attribs.getPermissions()), path);
				ProgrammerLog.log(logEntry);

				// Translate SFTP file attributes to FileInfo style attributes.
				int attributes = 0;
				attributes |= (rf.isDirectory() ? FileInfo.DIR_ATTRIB : 0);
				attributes |= (rf.isFile() ? FileInfo.FILE_ATTRIB : 0);
				attributes |= (rf.canRead() ? FileInfo.READ_ATTRIB : 0);
				attributes |= (rf.canWrite() ? FileInfo.WRITE_ATTRIB : 0);
				attributes |= (rf.canExecute() ? FileInfo.EXEC_ATTRIB : 0);

				// Now create the file info object.
				info = new FileInfo(path, attribs.getATime() * 1000L, attribs.getSize(), attributes);
			} else {
				// The file does not exist? Create a dummy info.
				ProgrammerLog.log("SFTP says " + path + " does not exit (but no exception).\n");
				info = new FileInfo(path, -1, -1, 0);
			}
		} catch (Exception e) {
			if (e.getMessage().indexOf("No such file") != -1) {
				// This is an expected exception if file does not exist
				ProgrammerLog.log("SFTP says " + path + " does not exit.\n");
				info = new FileInfo(path, -1, -1, 0);
			} else {
				// Unexpected exception
				throw new IOException(e);
			}
		}
		// Log to track translation.
		ProgrammerLog.log("FileInfo for " + path + " = " + info + "\n");
		// Return file information to caller.
		return info;
	}

	/**
	 * Creates a directory on the server.
	 * 
	 * This method must be used to create a directory entry on the server.
	 * 
	 * <p>
	 * <b>Note:</b> The connection to the server must have already been successfully
	 * established via a call to the connect method before invoking this method.
	 * </p>
	 * 
	 * @param directory The fully path to the directory to be created.
	 * 
	 * @throws IOException This method throws an exception if the directory could
	 *                     not be created.
	 */
	@Override
	public void mkdir(String directory) throws Exception {
		if (session == null) {
			throw new IOException("Remote session not connected.");
		}
		// Now that the session is connected try to run mkdir to create the specified
		// directory on the remote machine.
		ChannelSftp sftp = null;
		try {
			sftp = (ChannelSftp) session.openChannel("sftp");
			sftp.connect();
			sftp.mkdir(directory);
		} finally {
			if (sftp != null) {
				sftp.disconnect();
			}
		}
	}

	/**
	 * Deletes an <i>empty</i> directory on the server.
	 * 
	 * This method must be used to delete a directory entry on the server. The
	 * directory must be empty.
	 * 
	 * <p>
	 * <b>Note:</b> The connection to the server must have already been successfully
	 * established via a call to the connect method before invoking this method.
	 * </p>
	 * 
	 * @param directory The fully path to the directory to be deleted.
	 * 
	 * @throws IOException This method throws an exception if the directory could
	 *                     not be deleted.
	 */
	@Override
	public void rmdir(String directory) throws Exception {
		if (session == null) {
			throw new IOException("Remote session not connected.");
		}
		// Now that the session is connected try to run rmdir to delete the specified
		// directory on the remote machine
		ChannelSftp sftp = null;
		try {
			sftp = (ChannelSftp) session.openChannel("sftp");
			sftp.connect();
			sftp.rmdir(directory);
		} finally {
			if (sftp != null) {
				sftp.disconnect();
			}
		}
	}

	/**
	 * A simple method to set a purpose message for this session.
	 * 
	 * This method can be used to set up a purpose message for a server session. The
	 * purpose message is displayed to the user when prompting for inputs from the
	 * user for credentials. The message serves the purpose of appraising the user
	 * about the purpose of the session.
	 * 
	 * @param text This string is used to create a label (possibly with an icon on
	 *             it). So it can be plain text or HTML. If the message is long,
	 *             then ensure it is properly broken into multiple lines so that
	 *             dialog boxes don't get too large.
	 */
	@Override
	public void setPurpose(String text) {
		purpose = text;
	}

	/**
	 * Method to forward a local port to a remote host and port number.
	 * 
	 * <p>
	 * This method is a convenience method that is meaningful only for remote server
	 * sessions. This method uses SSH port forwarding support in JSch to forward
	 * connections from a local port to a given remote host and port.
	 * </p>
	 * 
	 * @param localPort  The local port on the local machine to be forwarded to a
	 *                   given remote host. If this value is -1, then a free port is
	 *                   detected by this method and used.
	 * 
	 * @param remoteHost The host or IP address of the remote machine to which a
	 *                   connection is to be forwarded.
	 * 
	 * @param remotePort The remote port number to which the connection is to be
	 *                   forwarded.
	 * 
	 * @return This method returns the local port that has been forwarded.
	 * @throws IOException
	 */
	@Override
	public int forwardPort(int localPort, String remoteHost, int remotePort) throws IOException {
		try {
			if (localPort == -1) {
				// Find out a local socket that is free.
				ServerSocket tempSocket = new ServerSocket(0);
				localPort = tempSocket.getLocalPort();
				tempSocket.close();
			}
			// Get JSch to forward the port for us.
			return session.setPortForwardingL(localPort, remoteHost, remotePort);
		} catch (Exception e) {
			ProgrammerLog.log(e);
			UserLog.log(LogLevel.WARNING, "RemoteServerSession", e.getMessage());
			throw new IOException(e);
		}
	}

	/**
	 * The session to the remote server via which the remote server can be accessed
	 * for performing various operations. The connection is created via the
	 * connect() method. A session consists of multiple, independent channels. Each
	 * channel performs different types of tasks.
	 */
	private Session session;

	/**
	 * The last known OS type for the remote server. This value is set after the
	 * getOSType() method is called. If a session is lost or closed, then this value
	 * is reset.
	 */
	private Server.OSType osType;

	/**
	 * A simple textual information indicating the purpose for this session. This
	 * string is more meaningful to the user and is merely used to provide
	 * additional information when prompting for inputs from the user.
	 */
	private String purpose;

	/**
	 * The top-level JSch object that is used to create sessions to multiple
	 * servers.
	 * 
	 * This is a top-level, shared JSch object that represents the top-level
	 * component in the JSch SSH-Client library. There is only one instance of this
	 * class that is used to create multiple, independent sessions to a server. Each
	 * session can have multiple, independent channels that are used to perform
	 * various tasks.
	 */
	private static final JSch sshClientLibrary = new JSch();

	/**
	 * The OS-specific path where the list of known hosts (that is the servers to
	 * which we have connected before) is stored. This list is applicable only for
	 * remote hosts.
	 */
	private static final String KNOWN_HOST_PATH = Utilities.getDefaultDirectory() + "/.KnownHosts";

	/**
	 * A generic informational message that is displayed to the user when any error
	 * occurs when trying to set the known hosts file to be used by JSch.
	 */
	private static final String KNOWN_HOSTS_ERROR = "<html>Error occured when attempting to use the known hosts<br>"
			+ "file: '" + KNOWN_HOST_PATH + "'.<br>"
			+ "Please rectify this issue appropriately. You can still continute<br>"
			+ "to use PEACE. But caching of known hosts for secure shell connection<br>will be disabled.";

	/**
	 * A first-time connection message that is formatted and displayed to the user.
	 * 
	 * This string contains an HTML message that is suitably formatted (to fill-in
	 * additional information) and displayed to the user. This message is used when
	 * the user connects to an Server for the first time and the server entry is not
	 * in the Known hosts file.
	 */
	private static final String FIRST_SSH_CONNECTION_MSG = "<html>The authenticity of host %s<br/>"
			+ "cannot be verified as it is not a \"known host\".<br/>The RSA fingerprint key is: %s<br/><br/>"
			+ "This is most likely because this is the first time you<br/>"
			+ "connecting to this host via PEACE and this message is<br/>"
			+ "normal when connecting via secure shell (SSH) protocol.<br/><br/>"
			+ "<b>Would you like to add this server to the \"known hosts\"<br/>"
			+ "and proceed with the connection?</b></html>";

	/**
	 * Message that is formatted and displayed to the user to warn about change in
	 * RSA finger print.
	 * 
	 * This string contains an HTML message that is suitably formatted (to fill-in
	 * additional information) and displayed to the user. This message is used when
	 * the user connects to an Server but the server's RSA finger print key has
	 * changed.
	 */
	private static final String SSH_HOST_CHANGE_MSG = "<html><b>The server's identification has changed!</b><br/>"
			+ "<b>It is possible that someone is doing something nasty</b><br/>"
			+ "(someone could be eavesdropping via man-in-the-middle type attack).<br/>"
			+ "It is aslo possible that the RSA host key for the server has changed.<br/><ul>"
			+ "<li>If this is a server maintained by your department or university<br/>"
			+ "it is normally safe to proceed with using the server.</li>"
			+ "<li>If not please contact your server administrators to verify that<br/>"
			+ "the change in finger print is expected prior to using the server.</li></ul>"
			+ "<b>Would you like to update the server's entry in \"known hosts\"<br/>"
			+ "and proceed with the connection?</b></html>";

	/**
	 * A static message that is included as a part of the RuntimeException generated
	 * by some of the methods in this class.
	 */
	private static final String USER_INTERRUPTED_EXP_MSG = "User has interrupted SSH connection to server";
}
