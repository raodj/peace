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
import java.awt.Color;
import java.awt.Component;
import java.awt.Dialog;

import javax.swing.BorderFactory;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.WindowConstants;

import org.peace_tools.generic.CustomBorder;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;


/**
 * A helper class to create a remote server session in background while
 * displaying a suitable message on the foreground.
 * 
 * This is a helper class that is used to establish connections with a
 * remote server in the background. Connecting to a remote server can
 * sometimes be a time consuming task that makes the GUI unresponsive
 * causing confusion to the user. This class ensures that the GUI remains
 * responsive while a connection is being established to a server in 
 * the background. This class can be used in the following manner:
 * 
 * <code>
 *   RemoteServerSessionMaker rssm = 
 * 					new RemoteServerSessionMaker(getWindow(), getServer(),
 * 							"Verifying & validating SSH-tunnel server");
 * 	 rssm.start();
 * 	 ServerSession session = rssm.getSession();
 * </code>
 * 
 * <p>Note that this class creates a dialog box and draws all the graphics
 * for the dialog box make it look something like a dialog box. However,
 * the actual dialog box is not drawn because the Close button on dialog
 * boxes cannot be disabled in Java.</p>
 */
public class RemoteServerSessionMaker extends SwingWorker<RemoteServerSession, RemoteServerSession> {
	/**
	 * The only constructor for this class.
	 * 
	 * The constructor sets up the dialog box that is displayed by this
	 * class. The dialog is actually made visible by the caller by calling
	 * the {@link #start()} method.
	 * 
	 * @param parent The parent window for the temporary dialog box
	 * displayed by this method.
	 * 
	 * @param srvr The server to which a connection is to be established.
	 * This should ideally be a remote server (otherwise connections to
	 * local servers are fast and this class is not really needed).
	 * 
	 * @param purpose The reason why the connection is being established.
	 * This reason is displayed to the user to provide additional 
	 * feedback.
	 */
	public RemoteServerSessionMaker(Component parent, Server srvr,
			String purpose) {
		this.srvr    = srvr;
		this.purpose = purpose;
		this.parent  = parent;
		// Create the progress bar and set its properties.
		JProgressBar progressBar = new JProgressBar();
		progressBar.setIndeterminate(true);
		progressBar.setStringPainted(true);
		progressBar.setString("Please wait...");
		// Create an informative label
		final String ConMsg = String.format(CONNECTING_MSG, 
				srvr.getName(), srvr.getPort());
		JLabel infoLabel = new JLabel(ConMsg, 
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		// Put the progress bar and infoLabel into a panel.
		JPanel infoPanel = new JPanel(new BorderLayout(10, 10));
		infoPanel.setBorder(BorderFactory.createEmptyBorder(10, 20, 10, 20));
		infoPanel.add(infoLabel, BorderLayout.NORTH);
		infoPanel.add(progressBar, BorderLayout.SOUTH);
		// Set up an overall top-level container panel.
		waitPanel = new JDialog(SwingUtilities.windowForComponent(parent), 
				"Connecting to Server...", Dialog.ModalityType.APPLICATION_MODAL);
		waitPanel.setUndecorated(true);
		waitPanel.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
		waitPanel.setResizable(false);
		// Use purpose to create a top-label to make it look nice
		JLabel title = new JLabel(purpose, JLabel.CENTER);
		final Color bgCol = progressBar.getForeground();
		title.setBackground(bgCol);
		title.setBorder(BorderFactory.createCompoundBorder(
				new CustomBorder(bgCol, bgCol, Color.BLACK, bgCol),
				BorderFactory.createEmptyBorder(3, 10, 3, 10)));
		title.setOpaque(true);
		// Create a top-level content pane to contain various GUI components.
		JPanel contentPane = new JPanel(new BorderLayout(10, 10));
		// Make a prettier border for the undecorated dialog
		//contentPane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(Color.BLACK),
		//		BorderFactory.createRaisedBevelBorder()));
		contentPane.setBorder(BorderFactory.createLineBorder(Color.BLACK));
				
		// Setup sub-elements
		contentPane.add(title, BorderLayout.NORTH);
		contentPane.add(infoPanel, BorderLayout.CENTER);
		// Add content pane to the dialog
		waitPanel.setContentPane(contentPane);
	}

	/**
	 * Determine if connection was successfully established and obtain
	 * connection.
	 * 
	 * <b>Note: This method blocks until the {@link #doInBackground()}
	 * method has completed its task.
	 *  
	 * @return This method returns a connected server session on success.
	 * On errors this method returns null.
	 */
	public RemoteServerSession getSession() {
		RemoteServerSession result = null;
		try {
			result = get();
		} catch (Exception e) {
			ProgrammerLog.log(e);
		}
		return result;
	}
	
	/**
	 * Method to start the actual work of creating a remote server session.
	 * 
	 * This method is a convenience method that must be used to 
	 * display the dialog (indicating work is underway) and
	 * actually start connecting to a server in the background.
	 */
	public void start() {
		// Finish layout.
		waitPanel.pack();
		waitPanel.setLocationRelativeTo(parent); 
		// Start background task.
		execute();
		// Make the wait dialog visible, blocking the caller.
		waitPanel.setVisible(true);
	}
	
	@Override
	protected RemoteServerSession doInBackground() throws Exception {
		// Try and connect to the remote server in the background.
		RemoteServerSession result = null;
		try {
			RemoteServerSession rss = (RemoteServerSession) SessionFactory.createSession(waitPanel, srvr);
			rss.setPurpose(purpose);
			rss.connect();
			result = rss;
		} catch (Exception e) {
			// Save exception for processing in the done() method
			exp = e;
		}
		return result;
	}

	/** Method to report results at the end of the test.
	 * 
	 * This method is automatically invoked from the main Swing thread 
	 * once the {@link #doInBackground()} method has completed. This
	 * method hides the {@link #waitPanel} dialog and displays the
	 * results from running MetaSim.
	 */
	@Override
	protected void done() {
		// Display failure information if any.
		if (exp != null) {
			final String errMsg = String.format(CONNECTING_ERROR_MSG,
					srvr.getName(), srvr.getPort());
			JPanel fullInfo = Utilities.collapsedMessage(errMsg, 
					Utilities.toString(exp));
			JOptionPane.showMessageDialog(waitPanel, fullInfo,
					"Error connecting to Server", JOptionPane.ERROR_MESSAGE);
		}
		// Hide the waitPanel.
		waitPanel.setVisible(false);
		waitPanel.dispose();
		super.done();
	}
	
	/**
	 * A simple dialog box requesting the user to wait while MetaSim is run.
	 * This dialog box provides the user with brief information about the
	 * fact that MetaSim is being run.
	 */
	private JDialog waitPanel;
	
	/**
	 * The exception that was generated (if any) when attempting to
	 * connect to the remote server.
	 * This is initialized to null and is later on
	 * changed if an exception occurs in the {@link #doInBackground()}
	 * method. It is displayed by the {@link #done()} method.
	 */
	private Exception exp = null;
	
	/**
	 * The server to which a remote session is to be established.
	 * This instance variable is initialized in the constructor and
	 * is never changed during the lifetime of this object.
	 */
	private final Server srvr;

	/**
	 * A simple textual information indicating the purpose for this
	 * session. This string is more meaningful to the user and is 
	 * merely used to provide additional information when prompting
	 * for inputs from the user.
	 */
	private final String purpose;
	
	/**
	 * The parent component that logically owns an instance of this
	 * class. This object is used to center the progress panel that
	 * is displayed by this class.
	 */
	private final Component parent;
	
	/**
	 * A text message that is suitably formatted and shown in the 
	 * dialog displayed by this class.
	 */
	private static final String CONNECTING_MSG =
		"<html>" +
		"Attempting to connect to server <br>" +
		"%s (Port: %d)<br>" +
		"This may take a few moments. Please wait..." +
		"</html>";
	
	/**
	 * A text message that is suitably formatted and shown to the 
	 * user if an error occurs when establishing connection.
	 */
	private static final String CONNECTING_ERROR_MSG =
		"<html>" +
		"Error occured when connecting to the following Server<br>" +
		"%s (Port: %d)<br>" +
		"Please verify host name, port, and login cerdentials." +
		"</html>";
}