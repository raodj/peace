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
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import javax.swing.JProgressBar;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.FileInfo;
import org.peace_tools.workspace.Server;

/**
 * A common base class for local and remote connections.
 * 
 * This is an abstract (non-instantiable) base class that serves as
 * a common interface for both local and remote server sessions. 
 * Sub-systems that exclusively operate with the ServerSession are
 * guaranteed to be agnostic to local and remote servers and
 * they would operate in a consistent manner. Local server sessions
 * typically directly run the commands on the local machine using
 * suitable JVM API calls. On the other hand remote sessions run
 * the commands on a remote machine via a secure shell (SSH) 
 * protocol.
 * 
 *  @see LocalServerSession
 *  @see RemoteServerSession
 */
public abstract class ServerSession {
	/**
	 * The constructor merely initializes the instance variables to
	 * the values supplied as parameters.
	 * 
	 * @param server The server data entry that provides the necessary
	 * information to connect to the server.
	 * @param parent The parent component that should be used to 
	 * create GUI elements that may be needed for any interactive
	 * operations.
	 */
	public ServerSession(Server server, Component parent) {
		this.server = server;
		this.parent = parent;
	}
	
	/**
	 * Connect to the server in order to perform various operations. 
	 * 
	 * This method must be used to establish a connection to a server
	 * before performing any tasks. This method is overridden in the
	 * derived class to perform the necessary connection operations.
	 * In the case of a remote session an actual network connection
	 * is established with the remote server which entails creating
	 * an SSH session and authenticating with the server. 
	 * 
	 * <p>In a single session, a connection needs to be established
	 * only once. After a connection has been established multiple
	 * commands and copy operations can be performed without requiring
	 * to connect and disconnect from the session. </p>
	 * 
	 * <p><b>Note:</b>  The process of establishing a connection can be a 
	 * time consuming task. In some cases incorrect host names can
	 * cause long connection times (until the connection times out
	 * which can be in minutes). Consequently, it is best to call
	 * this method from a separate daemon thread.</p>
	 * 
	 * @throws IOException This method throws IO exceptions in 
	 * the case of errors. If an error occurs, then a connection
	 * was not established and the caller will have to try again
	 * to establish a connection.
	 */
	public abstract void connect() throws IOException;
	
	/**
	 * Disconnect and close connection to a server.
	 * 
	 * This method must be used to disconnect an existing connection
	 * to a server. Invoking this method on session that is not
	 * connected has no side effects (and derived session classes
	 * do no special operations if a connection does not exist).
	 * 
	 * <p><b>Note:</b>  Disconnecting from a server that may be performing a 
	 * long running operation may have undesired side effects of 
	 * the commands aborting and possibly leaving data files in
	 * an inconsistent state.</p>
	 */
	public abstract void disconnect();
	
	/**
	 * This method can be used to run a brief command that 
	 * produces succinct output.
	 * 
	 * This method provides a convenient API to execute a 
	 * command on a server, buffer the output, and
	 * return the resulting output(s) as a string(s). Since this
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
	public abstract int exec(final String command, 
			String[] outputs) throws Exception;
	
	/**
	 * This method can be used to run a long running command that 
	 * may produce verbose output.
	 * 
	 * This method provides a convenient API to execute a 
	 * command on a server, stream the output, and
	 * return the exit code from the remote command. The
	 * outputs are streamed to a given StyledDocument. This
	 * method uses styles named "stdout", "stderr", and 
	 * "warning" (if available in the given output document)
	 * 
	 * <p><b>Note:</b>  The connection to the remote server must have
	 * been established successfully via a call to connect method.</p>
	 * 
	 * @param command The command line to be executed. This
	 * command must be compatible with the target machine's 
	 * OS. Otherwise this method will generate an exception.
	 * 
	 * @param output The styled document to which the standard
	 * output and standard error streams are to be written. This
	 * parameter cannot be null.
	 * 
	 * @return The exit code from the remote command that was
	 * run.
	 * 
	 * @throws Exception If the execution produces an error then
	 * this method throws an exception.
	 */
	public abstract int exec(final String command, 
			DefaultStyledDocument output) throws Exception;
	
	/**
	 * Determine the type of OS that this session is connected to.
	 * 
	 * This method can be used to determine the type of the OS that
	 * this session is connected to. 
	 * 
	 * <p><b>Note:</b>  The session must be connected prior to this call. 
	 * Implementations must try to make this method as quick as
	 * possible, by caching the OS type, if necessary.</p>
	 * 
	 * @return The type of OS this session is associated with.
	 * 
	 * @throws Exception If the detection of OS fails or a connection
	 * does not exist then this method throws an exception.
	 */
	public abstract Server.OSType getOSType() throws Exception;
	
	/**
	 * Copy given data from an input stream to a given file on the
	 * server.
	 *
	 * <p><b>Note:</b>  The connection to the server must have already been 
	 * successfully established via a call to the connect method 
	 * before invoking this method.</p>
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
	 * @param mode The POSIX compatible mode string (such as: 0600
	 * or 0777 to be used as the mode for the target file. Note the leading
	 * zero -- is is absolutely important and is needed to ensure that
	 * the number is in octal representation.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	public abstract void copy(InputStream srcData, 
			String destDirectory, String destFileName, 
			int mode) throws IOException;
	
	/**
	 * Copy file from a remote machine to a given output stream.
	 *
	 * <p><b>Note:</b>   The connection to the remote host must have already
	 * been established before invoking this method.</p>
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
	 * file copy progress. This parameter can be null.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	public abstract void copy(OutputStream destData, 
			String srcDirectory, String srcFileName, 
			JProgressBar progBar) throws IOException;
	
	/**
	 * Creates a directory on the server.
	 * 
	 * This method must be used to create a directory entry on the
	 * server.
	 * 
	 * <p><b>Note:</b>  The connection to the server must have already been 
	 * successfully established via a call to the connect method 
	 * before invoking this method.</p>
	 * 
	 * @param directory The fully path to the directory to be created.
	 * 
	 * @throws IOException This method throws an exception if the 
	 * directory could not be created.
	 */
	public abstract void mkdir(String directory) throws Exception;
	
	/**
	 * Remove an <i>empty</i> directory on the server.
	 * 
	 * This method must be used to remove a directory entry on the
	 * server. Note that the directory must be empty in order for
	 * this operation to succeed.
	 * 
	 * <p><b>Note:</b>  The connection to the server must have already been 
	 * successfully established via a call to the connect method 
	 * before invoking this method.</p>
	 * 
	 * @param directory The fully path to the directory to be deleted.
	 * 
	 * @throws IOException This method throws an exception if the 
	 * directory could not be deleted.
	 */
	public abstract void rmdir(String directory) throws Exception;

	/**
	 * Obtain information about a given path on the server.
	 * 
	 * <p><b>Note:</b>  This method uses SFTP to copy the data. The connection
	 * to the remote host must have already been established before
	 * invoking this method.</p>
	 * 
	 * @param path The path (absolute or relative to home directory)
	 * of the file whose meta data is to be retrieved.
	 * 
	 * @return A FileInfo object containing the information about
	 * the path.
	 * 
	 * @throws IOException This method throws exceptions on errors.
	 */
	abstract public FileInfo fstat(String path) throws IOException;
	
	/**
	 * Method to forward a local port to a remote host and port number.
	 * 
	 * <p>This method is a convenience method that is meaningful only for
	 * remote server sessions. This method  uses SSH port forwarding 
	 * support in JSch to forward connections from a local port to a 
	 * given remote host and port.</p>  
	 * 
	 * @param localPort The local port on the local machine to be forwarded
	 * to a given remote host. If this value is -1, then a free port is
	 * detected by this method and used.
	 * 
	 * @param remoteHost The host or IP address of the remote machine to
	 * which a connection is to be forwarded.
	 * 
	 * @param remotePort The remote port number to which the connection
	 * is to be forwarded.
	 * 
	 * @return This method returns the local port that has been forwarded.
	 * @throws IOException 
	 */
	public abstract int forwardPort(int localPort, String remoteHost, 
			int remotePort) throws IOException;
	
	/**
	 * A simple method to set a purpose message for this session.
	 * 
	 * This method can be used to set up a purpose message for a
	 * server session. The purpose message is typically displayed
	 * to the user when prompting for inputs from the user for 
	 * credentials. The message serves the purpose of appraising the
	 * user about the purpose of the session.
	 *  
	 * @param text This string is used to create a label (possibly
	 * with an icon on it). So it can be plain text or HTML. If the
	 * message is long, then ensure it is properly broken into 
	 * multiple lines so that dialog boxes don't get too large.
	 */
	abstract public void setPurpose(String text);
	
	/**
	 * The server data entry that provides the necessary information to 
	 * connect to the server.
	 */
	protected final Server server;
	
	/**
	 * The parent component that visually owns this server session. This
	 * parent component is used to display dialog boxes used for prompting
	 * the user for password etc (as needed) by the derived sessions.
	 */
	protected final Component parent;
}
