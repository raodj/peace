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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;

import javax.swing.JProgressBar;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.FileInfo;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

/**
 * A local server session to run jobs on the local host.
 * 
 * <p>This class provides an implementation of the ServerSession API.
 * Specifically, this class provides a session that can be used to
 * interact with the local PC using the same API as that used for
 * remote hosts. The consistent API eases development of the core GUI
 * modules.</p> 
 * 
 */
public class LocalServerSession extends ServerSession {
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
	protected LocalServerSession(Component parent, Server server) {
		super(server, parent);
		purpose = null;
	}

	/**
	 * Connect to the server in order to perform various operations. 
	 * 
	 * This method is really meaningful only for a remote server session
	 * that requires connections to be established. For a local server
	 * there is absolutely nothing to do and this method is empty.
	 * 
	 * @throws IOException The exception signature is only for API
	 * compliance. This method does not really throw any exceptions
	 */
	@Override
	public void connect() throws IOException {
	}
	
	/**
	 * Method to disconnect from a remote server. 
	 * 
	 * Similar to the connect method, this method is blank because there
	 * isn't a physical connection to tear down.
	 */
	@Override
	public void disconnect() {
	}
	
	/**
	 * This method can be used to run a brief command that 
	 * produces succinct output.
	 * 
	 * This method provides a convenient API to execute a 
	 * command on the local machine, buffer the output, and
	 * return the resulting output as a string. Since this
	 * method buffers and returns the output, it must be used
	 * only for jobs that return small volumes of data. In 
	 * addition, since the method returns the results only
	 * after the command has fully executed, it does not
	 * provide any interactive features to the user. 
	 * Consequently, it should be used only for jobs that
	 * run for a short period of time.
	 * 
	 * <p><b>Note:</b>  The connection to the remote server must have
	 * been established successfully via a call to connect method.</p>
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
		Process process = startProcess(command);
		// Read all the data into a single string from standard output
		// assuming standard error does not fill up.
		outputs[0] = Utilities.readFullStream(process.getInputStream());
		// Convert the standard output to a string as well.
		outputs[1] = Utilities.readFullStream(process.getErrorStream());
		// Return exit code for the process ensuring it has actually finished.
		return process.waitFor();
	}
	
	/**
	 * Helper method to start process to execute command via OS-specific
	 * command processor.
	 * 
	 * @param command The command to be run.
	 * @return The process object created by running the command via 
	 * OS-specific command processor.
	 */
	private Process startProcess(String command) throws Exception {
		// Start up the process via OS-specific command processor
		String cmdHandler = "/bin/bash";
		String parameter  = "-c";
		if (Server.OSType.WINDOWS.equals(getOSType())) {
			cmdHandler = "cmd.exe";
			parameter  = "/c";
		}
		Process process = new ProcessBuilder(cmdHandler, parameter, command).start();
		return process;
	}
	
	/**
	 * This method can be used to run a long running command that 
	 * may produce verbose output.
	 * 
	 * This method provides a convenient API to execute a 
	 * command on a local server, stream the output, and
	 * return the exit code from the process. The
	 * outputs are streamed to a given StyledDocument. This
	 * method uses styles named "stdout", "stderr", and 
	 * "warning" (if available in the given output document)
	 * 
	 * <p><b>Note:</b>  This method waits for the process to complete before
	 * returning control. Consequently, this method must be called
	 * from a separate thread if the GUI should not block.<p>
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
		// Start up the process.
		Process process = startProcess(command);
		// Process the output stream (assuming standard error does not fill up)
		// Read stdout one line at a time and add it to the output
		BufferedReader stdin = new BufferedReader(new InputStreamReader(process.getInputStream()));
		String line = null;
		while ((line = stdin.readLine()) != null) {
			// Got another line of standard output. Display it.
			output.insertString(output.getLength(), line + "\n", 
					output.getStyle("stdout"));
		}
		// Read any pending data on the standard error stream
		String stdErrData = Utilities.readFullStream(process.getErrorStream());
		output.insertString(output.getLength(), stdErrData, 
				output.getStyle("stderr"));
		// Return the exit code from the process
		return process.waitFor();
	}
	
	/**
	 * Determine the type of OS that this session is connected to.
	 * 
	 * This method can be used to determine the type of the OS that
	 * this session is connected to. 
	 * 
	 * @return The type of OS this session is associated with.
	 * 
	 * @throws Exception If the detection of OS fails or a connection
	 * does not exist then this method throws an exception.
	 */
	@Override
	public Server.OSType getOSType() throws Exception {
		String osName = System.getProperty("os.name");
		// Determine OS type based on the response string
		Server.OSType osType = Server.OSType.UNIX;
		if (osName.indexOf("Windows") != -1) {
			osType = Server.OSType.WINDOWS;
		}
		if (osName.indexOf("Linux") != -1) {
			osType = Server.OSType.LINUX;
		}
		return osType;
	}
	
	/**
	 * Copy given data from an input stream to a given file on the
	 * local machine.
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
	 * @param mode The POSIX compatible mode number (such as: 0600
	 * or 0777 to be used as the mode for the target file.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	public void copy(InputStream srcData, String destDirectory, 
			String destFileName, int mode) throws IOException {
		// Compute the final path and target file name.
		if (destDirectory == null) {
			destDirectory = ".";
		}
		final String targetFile = destDirectory + File.separator + destFileName;
		// Create a file to update attributes.
		File destFile = new File(targetFile);
		FileOutputStream os = new FileOutputStream(destFile);
		// Read data from input stream and write data to the destFile
		byte buffer[]  = new byte[8096];
		int  bytesRead = 0;
		while ((bytesRead = srcData.read(buffer, 0, buffer.length)) != -1) {
			os.write(buffer, 0, bytesRead);
		}
		// Close destination file
		os.close();
		// Update attributes
		setPerms(destFile, mode & 0700, true);  // owner
		// setPerms(destFile, mode.charAt(3), false); // others
	}
	
	/**
	 * Copy file from one directory to a given output stream.
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
		// Compute the final path and target file name.
		if (srcDirectory == null) {
			srcDirectory = ".";
		}
		final String sourceFile = srcDirectory + "/" + srcFileName; 
		// A try..finally block to ensure streams get closed
		FileInputStream srcFile = null;
		try {
			// Make the progress bar indeterminate for now.
			if (progBar != null) {
				progBar.setIndeterminate(true);
			}
			// Create an input stream
			File file = new File(sourceFile);
			srcFile = new FileInputStream(sourceFile);
			// Now determine the file size and set progress bar size
			if (progBar != null) {
				progBar.setMaximum((int) file.length());
				progBar.setIndeterminate(false);
			}
			// Read data from input stream and write data to the destFile
			byte buffer[]  = new byte[4096];
			int  bytesRead = 0; 
			long srcOffset = 0;
			while ((bytesRead = srcFile.read(buffer, 0, buffer.length)) != -1) {
				destData.write(buffer, 0, bytesRead);
				srcOffset += bytesRead;
				if (progBar != null) {
					progBar.setValue((int) srcOffset);
				}
			}
		} finally {
			// Ensure that the file is closed correctly
			if (srcFile != null) {
				srcFile.close();
			}
		}
	}
	
	/**
	 * Creates a directory on the server.
	 * 
	 * This method must be used to create a directory entry on the
	 * server.
	 * 
	 * @param directory The fully path to the directory to be created.
	 * 
	 * @throws IOException This method throws an exception if the 
	 * directory could not be created.
	 */
	public void mkdir(String directory) throws Exception {
		File dir = new File(directory);
		if (!dir.mkdir()) {
			String msg = "Unable to create directory '" + directory + "'";
			throw new IOException(msg);
		}
	}
	
	/**
	 * Deletes an <i>empty</i> directory on the server.
	 * 
	 * This method must be used to delete a directory entry on the
	 * server. The directory must be empty.
	 * 
	 * <p><b>Note:</b>The connection to the server must have already been 
	 * successfully established via a call to the connect method 
	 * before invoking this method.</p>
	 * 
	 * @param directory The fully path to the directory to be deleted.
	 * 
	 * @throws IOException This method throws an exception if the 
	 * directory could not be deleted.
	 */
	public void rmdir(String directory) throws Exception {
		File dir = new File(directory);
		if (!dir.delete()) {
			String msg = "Unable to create directory '" + directory + "'";
			throw new IOException(msg);
		}
	}

    /**
     * Obtain information about a given path on the server.
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
		File file = new File(path);
		if (!file.exists()) {
			// File does not exist. Return default information.
			return new FileInfo(path, -1, -1, 0);
		} 
		// File exists create complete file info
		int attributes = 0;
		attributes |= (file.isDirectory() ? FileInfo.DIR_ATTRIB   : 0);
		attributes |= (file.isFile()      ? FileInfo.FILE_ATTRIB  : 0);
		attributes |= (file.canRead()     ? FileInfo.READ_ATTRIB  : 0);
		attributes |= (file.canWrite()    ? FileInfo.WRITE_ATTRIB : 0);
		attributes |= (file.canExecute()  ? FileInfo.EXEC_ATTRIB  : 0);
		// Now create the file info object.
		return new FileInfo(path, file.lastModified(), file.length(), attributes);
	}

	/**
	 * Set up the permissions for a given file.
	 * 
	 * This is a helper method that uses a given digit in the permission
	 * string to setup the read, write, and execute information for a given
	 * file.
	 * 
	 * @param file The file whose permissions are to be set.
	 * 
	 * @param perm A POSIX compatible digit (0-7) that indicates the flags
	 * for read (4), write (2), and execute (1) values.
	 * 
	 * @param owner Flag to indicate if the status is for the owner (true) or 
	 * others (false). 
	 */
	private void setPerms(File file, int perm, boolean owner) {
		file.setReadable((perm & 0x4) != 0, owner);
		file.setWritable((perm & 0x2) != 0, owner);
		file.setExecutable((perm & 0x1) != 0, owner);
	}
	
	/**
	 * A simple method to set a purpose message for this session.
	 * 
	 * This method can be used to set up a purpose message for a
	 * server session. The purpose message is displayed
	 * to the user when prompting for inputs from the user.
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
	 * Obtain the purpose set for this session.
	 * 
	 * @return The purpose set for this session. If a purpose is
	 * not set then this method returns null.
	 */
	public String getPurpose() {
		return purpose;
	}
	
	/**
	 * A simple textual information indicating the purpose for this
	 * session. This string is more meanigful to the user and is 
	 * merely used to provide additional information when prompting
	 * for inputs from the user.
	 */
	private String purpose;
}
