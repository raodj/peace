package org.peace_tools.views;

import java.awt.BorderLayout;
import java.util.Date;

import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.MSTClusterData;
import org.peace_tools.workspace.MSTData;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/** Class to display details about a job.
 * 
 * This class extends the DetailView base class and provide the 
 * necessary functionality to display details about the job. 
 * This class utilizes most of the infrastructure provided
 * by the base class and the core work is done in the
 * run method.
 */
public class JobDetailsView extends DetailView {
	/**
	 * The constructor.
	 * 
	 * The constructor initializes the various instance variables and
	 * creates the different tabs and panels.
	 *  
	 * @param job The job entry that is being displayed in this panel. 
	 * 
	 * @param server The server on which the job information resides.
	 */
	public JobDetailsView(Job job, Server server) {
		super("Job ID: ", job.getJobID(), server);
		this.job    = job;
		// Setup the panels
		outputs = new DefaultStyledDocument[4];
		String TabNames[] = {"Overview", "Outputs", "Errors", "Scripts"};
		add(createOutputDocs(outputs, TabNames), BorderLayout.CENTER);
		// Populate the first overview tab with static information.
		setJobOverview();
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
			progressBar.setMinimum(0);
			progressBar.setMaximum(7);
			// Start up initial progress
			progressBar.setEnabled(true);
			progressBar.setValue(0);
			// Reset the dynamic information.
			resetDynamicInfo();
			// First try and connect to the server.
			ServerSession session = SessionFactory.createSession(this, server);
			// Setup purpose so the user knows why we are connecting
			session.setPurpose("<html>Attempting to obtain/update " +
					"information for job " + job.getJobID() + "<br>" +
					"running on server " + server.getName() + "</html>");
			// Connect to the server.
			session.connect();
			progressBar.setValue(1);
			// Obtain stdout information first.
			String[] streams = {"", ""};
			final String script = job.getPath() + "/jobRunner";
			runCommand(session, script, "output", streams);
			// Tick progress
			progressBar.setValue(2);
			// Add outputs to the appropriate styled document.
			updateDocument(outputs[1], DocTitles[1], streams[0], "stdout");			
			// Next obtain stderr information.
			runCommand(session, script, "error", streams);
			progressBar.setValue(3);
			// Add error info to the appropriate styled document.
			updateDocument(outputs[2], DocTitles[2], streams[0], "stderr");
			progressBar.setValue(4);
			// Next get script information.
			runCommand(session, script, "scripts", streams);
			progressBar.setValue(5);
			updateDocument(outputs[3], DocTitles[3], streams[0], "scripts");
			progressBar.setValue(6);
			// Wind up the session.
			session.disconnect();
			progressBar.setValue(7);
		} catch (Exception exp) {
			String errMsg = String.format(INFO_ERR_MSG, job.getJobID(),
					server.getName());
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
		updateDocument(outputs[1], DocTitles[1], message, "subtitle");
		updateDocument(outputs[2], DocTitles[2], message, "subtitle");
		updateDocument(outputs[3], DocTitles[3], message, "subtitle");
	}

	/**
	 * Helper method to populate first tab in job details.
	 * 
	 * This method is called from the constructor to populate the basic
	 * information about a job entry in the job tab. The necessary information
	 * is directly obtained from the work space information. This method was
	 * primarily introduced to streamline the code in the constructor.
	 * 
	 * <p>
	 * <b>Note:</b> This method assumes that
	 * {@link #createOutputDocs(DefaultStyledDocument[], String[])} method has
	 * already been called to create the tabs where the information is to be
	 * displayed.
	 * </p>
	 */
	private void setJobOverview() {
		try {
			assert( outputs[0] != null );
			// Add the output in given style
			outputs[0].insertString(outputs[0].getLength(), 
					"Overview\n", 
					outputs[0].getStyle("title"));
			// Build the overview information
			// Append information about the server and nodes.
			outputs[0].insertString(outputs[0].getLength(), "\nServer Information:\n", 
					outputs[0].getStyle("subtitle"));
			String serverInfo = String.format(SERVER_INFO, 
					server.getName(), job.getNodes(), job.getCPUsPerNode(),
					job.getMemory(), job.getMaxRunTime());
			outputs[0].insertString(outputs[0].getLength(), serverInfo, 
					outputs[0].getStyle("overview"));
			outputs[0].insertString(outputs[0].getLength(), "\nMST Data Summary:\n", 
					outputs[0].getStyle("subtitle"));
			// The MST summary information.
			MSTData mst = Workspace.get().getMSTData(job.getJobID());
			if (mst != null) {
				outputs[0].insertString(outputs[0].getLength(), mst.getSummary("\t"), 
						outputs[0].getStyle("overview"));
			}
			// Append clustering information/parameters.
			outputs[0].insertString(outputs[0].getLength(), "\nCluster Data Summary:\n", 
					outputs[0].getStyle("subtitle"));
			MSTClusterData cluster = Workspace.get().getClusterData(job.getJobID());
			if (cluster != null) {
				outputs[0].insertString(outputs[0].getLength(), cluster.getSummary("\t"), 
					outputs[0].getStyle("overview"));
			}
		} catch (BadLocationException e) {
			ProgrammerLog.log(e);
		}
	}
	
	/**
	 * The styled documents that is used to display the output from
	 * the job information.
	 */
	private DefaultStyledDocument outputs[];
	
	/**
	 * The job whose details is being displayed by this panel.
	 * This value is set in the constructor and is never changed
	 * during the lifetime of this class. 
	 */
	final private Job job;
	
	/**
	 * A format string to format the information on a Job/server.
	 * This string is used when populating information about
	 * the job in the setJobOverview() method.
	 */
	private static final String SERVER_INFO = 
		"\tServer: %s\n" + 
		"\tNodes: %d\n" +
		"\tCPUs per Node: %d\n" +
		"\tMax Memory (GB): %d\n" + 
		"\tMax run time (hours): %d\n";
	
	/**
	 * The title strings associated with the 4 documents. The strings
	 * are used to reset the document information and to pouplate
	 * information in them.
	 */
	private static final String[] DocTitles = 
	{"Job Overview/Summary", "Output Stream (stdout)", 
		"Error Stream (stderr)", "Scripts & Files"	};
	
	/**
	 * A formattable string that is formatted and displayed to the
	 * user as an error message.
	 */
	private static final String INFO_ERR_MSG = 
		"<html>Unable to determine all the information about job %s<br>" +
		"from server %s<br>" +
		"because an error occured. Only partial information was obtained.<br>" +
		"See details below for more information.<html>";
	
	/**
	 * A generated serialization UID to keep compiler happy.
	 */
	private static final long serialVersionUID = -20535421438589052L;
}
