package org.peace_tools.views;

import java.awt.BorderLayout;
import java.util.Date;

import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

/** Class to display list of jobs/processes running on a server.
 * 
 * This class extends the DetailView base class and provide the 
 * necessary functionality to display details about the job. 
 * This class utilizes most of the infrastructure provided
 * by the base class and the core work is done in the
 * run method.
 */
public class ServerJobsView extends DetailView {
	/**
	 * The constructor.
	 * 
	 * The constructor initializes the various instance variables and
	 * creates the different controls.
	 *  
	 * @param server The server whose jobs/processes are to be listed.
	 * 
	 * @param userID The ID of the user whose jobs are to be listed.
	 * This parameter is optional and can be null.
	 */
	public ServerJobsView(Server server, String userID) {
		super("User ID: ", (userID != null ? userID : "<n/a>"), server);
		this.userID = userID;
		// Setup the panels
		DefaultStyledDocument outputs[] = new DefaultStyledDocument[1];
		String TabNames[] = {"Jobs/Processes"};
		// The base class methods creates a JTabbedPane. However, we
		// have only a single tab. So we remove the component out
		// of the tab pane.
		JTabbedPane tabPane = createOutputDocs(outputs, TabNames);
		add(tabPane.getComponent(0), BorderLayout.CENTER);
		jobInfo = outputs[0];
		// Start the refresh process.
		actionPerformed(null);
	}
	
	/**
	 * The thread method that retrieves and populates job information.
	 * 
	 * This method is a long running method that performs the task
	 * of interacting with the server to obtain job information.
	 * This method must be called from a separate thread. If any
	 * error occurs this method displays it in a dialog box. At the 
	 * end of the update process this method re-enables the "Refresh" 
	 * button.
	 */
	@Override
	public void run() {
		try {
			// Setup suitable range for progress bar
			progressBar.setMinimum(0); progressBar.setMaximum(4);
			// Start up initial progress
			progressBar.setEnabled(true);
			progressBar.setValue(0);
			// Reset the dynamic information.
			resetDynamicInfo();
			// First try and connect to the server.
			ServerSession session = SessionFactory.createSession(this, server);
			// Setup purpose so the user knows why we are connecting
			session.setPurpose("<html>Attempting to obtain/update " +
					"jobs/process information<br>" +
					"running on server " + server.getName() + "</html>");
			// Connect to the server.
			session.connect();
			progressBar.setValue(1);
			// Obtain list of jobs running on the server
			String[] streams = {"", ""};
			final String script = server.getInstallPath() + "/listJobs";
			runCommand(session, script, (userID != null ? userID : ""), streams);
			// Tick progress
			progressBar.setValue(2);
			// Add outputs to the appropriate styled document.
			updateDocument(jobInfo, "Job/Process Information", 
					streams[0], "stdout");
			progressBar.setValue(3);
			// Check for empty outputs and let user know things are fine.
			if ((userID != null) && (streams[0].trim().isEmpty())) {
				// Empty output for user jobs. This can happen if the
				// user does not have any jobs pending.
				updateDocument(jobInfo, "Job/Process Information", 
						NO_JOB_INFO, "stdout");				
			}
			// Wind up the session.
			session.disconnect();
			progressBar.setValue(4);
		} catch (Exception exp) {
			String errMsg = String.format(INFO_ERR_MSG, server.getName());
			JPanel fullMsg = Utilities.collapsedMessage(errMsg, 
					Utilities.toString(exp));
			JOptionPane.showMessageDialog(this, fullMsg, 
					"Unable to obtain full job information", 
					JOptionPane.ERROR_MESSAGE);
		} finally {
			refreshButton.setEnabled(true);
			progressBar.setEnabled(false);
		}
	}
		
	/**
	 * Helper method to reset outputs, errors, and scripts to "unavailable"
	 * status.
	 * 
	 * This method is a helper method that is used to set/reset the 
	 * information in the dynamically loaded tabs to their initial
	 * state indicating the information is unavailable. 
	 */
	private void resetDynamicInfo() throws BadLocationException {
		Date now = new Date();
		String message = "\nThis information is currently unavailable.\n" +
			"Click the Refresh button to obtain/update information.\n"    +
			"\nTimestamp: " + now + "\n";
		// Update time stamp when information was last updated.
		timestamp.setText(now.toString());
		// Reset information in the three dynamic tags
		updateDocument(jobInfo, "Job/Process Information", 
				message, "subtitle");
	}
	
	/**
	 * The styled documents that is used to display the list of 
	 * jobs/processes running on a process.
	 */
	private DefaultStyledDocument jobInfo;
	
	/**
	 * The ID of the user whose jobs/processes are to be listed
	 * in this view.
	 */
	final String userID;
	
	/**
	 * A formattable string that is formatted and displayed to the
	 * user as an error message.
	 */
	private static final String INFO_ERR_MSG = 
		"<html>Unable to determine all the information about jobs<br>" +
		"or processes from server %s<br>" +
		"because an error occured. Only partial information was<br>" +
		"obtained. See details below for more information.<html>";
	
	/**
	 * A simple message that is displayed to the user if there is
	 * no job information available.
	 */
	private static final String NO_JOB_INFO = "\nIt appears that none" +
			"of your jobs/process are currently\n" +
			"running on this server. Therefore there is no job/process\n" +
			"information to display at this time.";
	
	/**
	 * A generated serialization UID to keep compiler happy.
	 */
	private static final long serialVersionUID = -20535421438589052L;
}
