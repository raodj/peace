package org.peace_tools.core;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.MSTClusterData;
import org.peace_tools.workspace.MSTData;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.ServerList;
import org.peace_tools.workspace.Workspace;

public class DeleteDialog extends JDialog 
implements ActionListener, Runnable {
	/**
	 * The constructor.
	 * 
	 * The constructor merely passes the parameters to the base class and
	 * sets up the various components to be displayed in this dialog.
	 *  
	 * @param mainFrame The main frame that ultimately owns all views and
	 * this dialog box. This information is required to correctly create a
	 * modal dialog that does not steal focus from other applications.
	 * 
	 * @param wsEntry The work space entry to be deleted. This object must
	 * be a valid work space entry such as: Job, Server, MSTClusterData,
	 * MSTData, or DataSet.
	 * 
	 */
	public DeleteDialog(MainFrame mainFrame, Object wsEntry) {
		super(mainFrame, "Verify entry deletion", true);
		// Setup final references to the parameters
		this.mainFrame = mainFrame;
		this.wsEntry   = wsEntry;
		// Ensure wsEntry is valid
		final int typeCode = getTypeCode(this.wsEntry); 
		assert(typeCode != -1);
		// Set the content pane for the dialog.
		JPanel contentPane = new JPanel(new BorderLayout(15, 5));
		setContentPane(contentPane);
		contentPane.setBorder(new EmptyBorder(0, 20, 10, 20));
		// Setup the top-level icon
		icon = new JLabel(UIManager.getIcon("OptionPane.warningIcon"));
		contentPane.add(icon, BorderLayout.WEST);
		
		// Setup the content panel with the necessary message
		JPanel msgPanel = new JPanel(new BorderLayout(0, 0));
		// Create entry specific custom message.
		String msg = String.format(Message[typeCode], wsEntry.toString());
		msg        = "<html>" + msg + "</html>";
		mainMessage= new JLabel(msg);
		msgPanel.add(mainMessage, BorderLayout.CENTER);
		// Create the check box
		deleteFiles = new JCheckBox("Delete physical file or directory as well?", true);
		msgPanel.add(deleteFiles, BorderLayout.SOUTH);
		// Add message panel to the dialog box.
		contentPane.add(msgPanel, BorderLayout.CENTER);
		
		// Create the button panel 
		JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 30, 10));
		buttonPanel.add(Utilities.createButton(null, "  OK  ", "ok", this, 
				"Proceed with the deletion operation", true));
		buttonPanel.add(Utilities.createButton(null, "Cancel", "cancel", this, 
				"Don't remove any entries.", true));
		
		// Create the indeterminate progress bar and associated message.
		JPanel progPanel = new JPanel(new BorderLayout(20, 5));
		progPanel.setBorder(new EmptyBorder(20, 5, 0, 20));
		delPath = new JLabel("Removing: ", JLabel.CENTER); 
		progPanel.add(delPath, BorderLayout.CENTER);
		progPanel.add(new JLabel(Utilities.getIcon("images/16x16/HourGlass.png")), BorderLayout.WEST);
		progressBar = new JProgressBar(JProgressBar.HORIZONTAL);
		progPanel.add(progressBar, BorderLayout.SOUTH);
		
		// Finally put the buttonPanel and progPanel in a card layout
		// so that they can switched easily later on.
		bottomPanel = new JPanel(new CardLayout());
		bottomPanel.add(buttonPanel, "buttons");
		bottomPanel.add(progPanel, "progbar");
		// Add bottom panel to the south of the dialog
		contentPane.add(bottomPanel, BorderLayout.SOUTH);
		
		// Now pack the dialog preparing it for display
		pack();
		// Lock in the size for this dialog.
		this.setResizable(false);
		// Set default close operation
		this.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		// Center the dialog on the main frame.
		setLocationRelativeTo(this.mainFrame);
	}

	/**
	 * An internal method to returns a type code for a given work space
	 * entry.
	 * 
	 * This is a helper method that is used to determine a type code
	 * for a work space entry object. The type code is simply an integer
	 * value for a given work space object. Integer values ease performing
	 * various tasks in this dialog. The integer values are assigned
	 * as per the following table:
	 * 
	 * <table>
	 * 	<th><td>Object Type</td>       <td>Code</td></tr>
	 *  <tr><td>Job</td>               <td>0</td></tr>
	 *  <tr><td>Server</td>            <td>1</td></tr>
	 *  <tr><td>MSTClusterData</td>    <td>2</td></tr>
	 *  <tr><td>MSTData</td>           <td>3</td></tr>
	 *  <tr><td>DataSet</td>           <td>4</td></tr>
	 * </table>
	 * 
	 * @param wsObject The wsObject whose type code is to be returned.
	 * 
	 * @return The type code for the given wsObject. If the object is 
	 * not a valid work space entry, then this method returns -1.
	 */
	private int getTypeCode(Object wsObject) {
		if (wsObject instanceof Job) {
			return 0;
		} else if (wsObject instanceof Server) {
			return 1;
		} else if (wsObject instanceof MSTClusterData) {
			return 2;
		} else if (wsObject instanceof MSTData) {
			return 3;
		} else if (wsObject instanceof DataSet) {
			return 4;
		}
		return -1;
	}

	/**
	 * Handle clicking of OK or Cancel buttons.
	 * 
	 * This method intercepts the click from OK or Cancel button
	 * and appropriately performs the necessary operation.
	 * 
	 * @param event The action event to be processed by this method.
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		if ("cancel".equals(event.getActionCommand())) {
			// The user canceled the operation. Simply close dialog
			disposeDialog();
			return;
		}
		// When control drops here the user really wants to delete
		// entries. Setup dialog to indicate the action is in progress
		CardLayout cl = (CardLayout) bottomPanel.getLayout();
		cl.last(bottomPanel);
		progressBar.setIndeterminate(true);
		icon.setIcon(UIManager.getIcon("OptionPane.informationIcon"));
		this.deleteFiles.setEnabled(false);
		this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		// Update main message to look grey'd out.
		String msg = mainMessage.getText();
		mainMessage.setText(Utilities.enableHTMLText(msg, false));
		// Finally start the actual deletion process from a background
		// thread so that the GUI does not appear to hang.
		Thread delThread = new Thread(this);
		delThread.start();
	}

	/**
	 * This method is called when the user clicks the OK button to
	 * delete entry(s). This method is invoked from a separate thread
	 * so that the GUI does not appear to hang (while files are being
	 * deleted on remote machines which could be a time consuming
	 * operation).
	 */
	@Override
	public void run() {
		final int entryType = getTypeCode(wsEntry);
		// First try and delete files if requested.
		if (this.deleteFiles.isSelected()) {
			try {
				// Different helper methods are used depending on
				// the entry type.
				switch (entryType) {
				case 0: // job entry
					delJobFiles((Job) wsEntry);
					break;
				case 1: // server entry
					Server server = (Server) wsEntry;
					deleteRemoteDir(server, server.getInstallPath());
					break;
				case 2: // MSTClusterData
					MSTClusterData cluster = (MSTClusterData) wsEntry;
					deleteLocalFile(cluster.getPath());
					break;
				case 3: // MSTData
					MSTData mst = (MSTData) wsEntry;
					deleteLocalFile(mst.getPath());
					break;
				}
			} catch (Exception e) {
				String msg = String.format(FILE_DEL_ERROR, wsEntry.toString());
				JPanel info = Utilities.collapsedMessage(msg, Utilities.toString(e));
				int choice = JOptionPane.showConfirmDialog(this,
						info, "Proceed further?", JOptionPane.WARNING_MESSAGE, 
						JOptionPane.YES_NO_OPTION);
				if (choice == JOptionPane.NO_OPTION) {
					disposeDialog();
					return;
				}				
			}
		}
		// Next, remove the actual entry from the work space.
		switch (entryType) {
		case 0: // job entry
			Workspace.get().getJobList().remove((Job) wsEntry);
			break;
		case 1: // server entry
			Workspace.get().getServerList().remove((Server) wsEntry);
			break;
		case 2: // MSTClusterData
			MSTClusterData cluster = (MSTClusterData) wsEntry;
			cluster.getDataSet().remove(cluster);
			break;
		case 3: // MSTData
			MSTData mst = (MSTData) wsEntry;
			mst.getDataSet().remove(mst);
			break;
		case 4: // Data set
			removeDataSet((DataSet) wsEntry);
		}
		// Save the updated work space information
		mainFrame.saveDelayedWorkspace();
		// All done.
		disposeDialog();
	}

	/**
	 * Helper method to hide and dispose this dialog.
	 * 
	 * This method simply hides and disposes the dialog from the
	 * main swing thread.
	 */
	private void disposeDialog() {
		// Now reset and dispose the GUI as removal is all done.
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				DeleteDialog.this.setVisible(false);
				DeleteDialog.this.dispose();
			}			
		});
	}
	
	/**
	 * Helper method to delete a data set along with all underlying 
	 * entries and files. This is a helper method that is invoked 
	 * from the run() method to delete  a complete data set. This
	 * method first deletes all the cluster entries (deleting
	 * files if indicated) followed by all MST entries (deleting
	 * MST files if indicated). If all the sub-entries were deleted
	 * successfully then this method also removes the data set 
	 * entry. 
	 * 
	 * @param dataSet The data set to be removed from the work space.
	 */
	private void removeDataSet(DataSet dataSet) {
		// First try and remove all cluster entries.
		ArrayList<MSTClusterData> clusters = new ArrayList<MSTClusterData>(dataSet.getClusterList());
		for (MSTClusterData clsEntry: clusters) {
			try {
				// Delete physical file for the cluster if requested.
				if (deleteFiles.isSelected()) {
					deleteLocalFile(clsEntry.getPath());
				}
				// Remove the cluster entry from the data set.
				dataSet.remove(clsEntry);
			} catch (Exception e) {
				// Error occurred. See if user want's to continue.
				String msg = String.format(DATASET_ERROR, clsEntry.toString());
				JPanel info = Utilities.collapsedMessage(msg, Utilities.toString(e));
				int choice = JOptionPane.showConfirmDialog(this,
						info, "Proceed further?", JOptionPane.QUESTION_MESSAGE, 
						JOptionPane.YES_NO_OPTION);
				if (choice == JOptionPane.NO_OPTION) {
					return;
				}
			}
		}
		// Next remove all MST entries from the data set.
		ArrayList<MSTData> mstList = new ArrayList<MSTData>(dataSet.getMSTList());
		for (MSTData mstEntry: mstList) {
			try {
				// Delete physical file for the MST if requested.
				if (deleteFiles.isSelected()) {
					deleteLocalFile(mstEntry.getPath());
				}
				// Remove the cluster entry from the data set.
				dataSet.remove(mstEntry);
			} catch (Exception e) {
				// Error occurred. See if user want's to continue.
				String msg = String.format(DATASET_ERROR, mstEntry.toString());
				JPanel info = Utilities.collapsedMessage(msg, Utilities.toString(e));
				int choice = JOptionPane.showConfirmDialog(this,
						info, "Proceed further?", JOptionPane.QUESTION_MESSAGE, 
						JOptionPane.YES_NO_OPTION);
				if (choice == JOptionPane.NO_OPTION) {
					return;
				}
			}
		}
		// Now finally remove the data set entry itself.
		Workspace.get().removeDataSet(dataSet);
	}
	
	/**
	 * Helper method to delete files for a given job.
	 * 
	 * This helper method was introduced to streamline the code in
	 * the run method. Several operations need to be performed to
	 * remove a job entry which includes locating a suitable server
	 * entry.
	 * 
	 * @param job The job entry to be removed.
	 */
	private void delJobFiles(Job job) throws Exception {
		// First locate a server entry in the workspace for this job
		ServerList srvrList = Workspace.get().getServerList();
		Server server = srvrList.getServer(job.getServerID());
		if (server == null) {
			// Could not locate a server entry. Can't really delete
			// files. Ask the user what to do.
			String msg = String.format(CANT_FIND_SERVER, job.getJobID());				
			int choice = JOptionPane.showConfirmDialog(this, msg, 
					"Can't find server", JOptionPane.YES_NO_OPTION, 
					JOptionPane.QUESTION_MESSAGE);
			if (choice == JOptionPane.NO_OPTION) {
				throw new Exception("Unable to delete files corresponding to job " +
						job.getJobID() + " as server entry was not found.");
			}
			// Return without doing anything further
			return;
		}
		// Now that we have a valid server entry, try and delete the
		// job directory on the server.
		deleteRemoteDir(server, job.getPath());
	}
	
	/**
	 * Helper method to delete a local file.
	 * 
	 * This is a simple method that deletes a local file on the
	 * local machine.
	 * 
	 * @param fileName The full path to the file to be deleted.
	 * @throws Exception This method throws an exception on errors.
	 */
	private void deleteLocalFile(String fileName) throws Exception {
		// Set information about file being deleted.
		delPath.setText("Deleting: " + fileName + " ...");
		File file = new File(fileName);
		if (!file.exists()) {
			// Nothing to delete if file is not there
			return;
		}
		if (!file.delete()) {
			throw new Exception("Unable to delete file:\n" + fileName);
		}
	}
	
	/**
	 * Helper method to recursively delete file on server.
	 * 
	 * This is a helper method that is used to recursively delete
	 * directories on a local or remote server. This method executes
	 * different commands (depending on server OS type) to recursively
	 * delete directories.
	 *
	 * @param server The server entry which the file is to be deleted.
	 * 
	 * @param directory The path to the directory to be deleted.
	 *
	 * @throws Exception This method throws an exception on errors.
	 */
	private void deleteRemoteDir(Server server, String directory) throws Exception {
		delPath.setText("Connecting to " + server.getName());
		// Create and establish a connection with the server
		ServerSession session = SessionFactory.createSession(this, server);
		// Setup purpose for the session if user is prompted
		session.setPurpose("<html>Attempting to delete the directory<br>" +
				"<i>" + directory + "</i><br>" +
				"as part of the removing the workspace entry.");
		// Try to connect (with interactive prompts if needed) to the server.
		session.connect();
		// OK, try and delete the file
		delPath.setText("Deleting: " + directory + " ...");
		// Setup the command to recursively delete directory based on
		// the OS of the server. For this we first determine OS type. 
		// Needless to say It is either windows or Linux
		String rmDirCmd = "rm -rf " + directory;
		ServerSession.OSType osType = session.getOSType();
		if (ServerSession.OSType.WINDOWS.equals(osType)) {
			rmDirCmd = "cmd /Q /C rmdir /s /q " + directory;
		}
		// Execute the command to recursively delete directories
		// Array to hold result from remote command
		String streamsData[] = {"", ""};
		if ((session.exec(rmDirCmd, streamsData) != 0) ||
				(streamsData[0] == null) || 
				(streamsData[1].length() > 0)) {
			// Unable to delete the directory. This is bad.
			throw new Exception("Unable to delete directory: " + 
					directory + "\non server: " + server.getName() + "\n" +
					"The server reported that:\n" + streamsData[1]);
		}
		// The files were successfully deleted.
		session.disconnect();
	}
	
	/**
	 * A convenient reference to the main frame that logically owns
	 * this dialog.
	 */
	private final MainFrame mainFrame;

	/**
	 * The actual work space entry to be removed by this dialog. This
	 * object is a valid work space entry such as: Job, Server, 
	 * MSTClusterData, MSTData, or DataSet. However, different actions
	 * and messages are generated depending on the entry being deleted.
	 */
	private final Object wsEntry;

	/**
	 * A check box to determine if the user wants the physical file to
	 * be deleted as well. If this check box is checked, the physical
	 * files associated with the entry are deleted first before the
	 * entry is deleted.
	 */
	private JCheckBox deleteFiles;

	/**
	 * The icon that is displayed in this dialog box to the left.
	 * The icon is stored separately to provide the option of 
	 * changing it at a later date.
	 */
	private JLabel icon;

	/**
	 * The progress bar that is placed in indeterminate mode when
	 * deleting remote entries and directories.
	 */
	private JProgressBar progressBar;

	/**
	 * The bottom panel in the dialog box that contains both the
	 * OK/Cancel buttons and the progress bar. Only one of them
	 * is shown as this panel is set to use a CardLayout manager.
	 */
	private JPanel bottomPanel;

	/**
	 * This label is used to display the actual path being deleted
	 * when physical file deletion is requested. Initially it is
	 * set to being empty.
	 */
	private JLabel delPath;

	/**
	 * This label contains the main message displayed in this dialog.
	 * This label's message is updated when the deletion process starts
	 * to make it appear disabled.
	 */
	private JLabel mainMessage;
	
	/**
	 * A set of predefined static messages to be displayed to the
	 * user indicating the entry to be removed and the implications
	 * of removing the entry. Begin and end HTML tags are added after
	 * the message is formatted in the constructor.
	 */
	private static final String Message[] = {
		"Are you sure you would like to remove<br>" +
		"the job entry with <i>Job ID: %s</i>?<br><br>" +
		"The job entry will be removed from the workspace.<br>" +
		"If you choose to delete the physical files then all<br>" +
		"the job information will be permanently lost.",
		"Are you sure you would like to remove the entry for<br>" +
		"server: <i>%s</i>?<br><br>" +
		"All jobs and job information running on this server will be lost.<br>" +
		"(<i>ensure there are no pending jobs on this server</i>)<br><br>" +
		"If you choose to delete the physical files then PEACE<br>" +
		"will be uninstalled from the server.",
		"Are you sure you would like to remove the entry<br>" +
		"<i>%s</i>?<br>",
		"Are you sure you would like to remove the entry<br>" +
		"<i>%s</i>?<br>",
		"Are you sure you would like to remove the data set entry<br>" +
		"<i>%s</i>?<br><br>" +
		"All underlying MST and clustering entries will also be removed.<br>" +
		"If you choose to delete the physical files then the MST and <br>" +
		"clustering data files will also be deleted."
	};

	/**
	 * A message that is formatted and displayed to the user when an 
	 * error occurs when deleting files or directories corresponding
	 * to a given entry in the work space.
	 */
	private static final String FILE_DEL_ERROR = 
		"<html>" +
		"An error occured when attempting to delete the file(s) for<br>" +
		"workspace entry: <i>%s</i><br>" +
		"Click on details below for additional information.<br><br>" +
		"Do you want to proceed with removing this entry from the workspace?" +
		"</html>";
	
	/**
	 * A message that is formatted and displayed to the user when a 
	 * server corresponding to a given job entry could not be located
	 * in the work space.
	 */
	private static final String CANT_FIND_SERVER = 
		"<html>The server entry for job <i>%s</i> could not be located.<br>" +
		"The files associated with this job cannot be deleted.<br>" +
		"(there is no point in keeping this job entry in the workspace)<br>" +
		"Do you still wish to proceed with removing this job entry?</html>";
	
	/**
	 * Message to be displayed to the user (after formatting) to report
	 * error when deleting a workspace entry.
	 */
	private static final String DATASET_ERROR = 
		"<html>Error occured when attempting to remove entry:<br>" +
		"%s<br><i>See details below for more information.<br><br>" +
		"Do you want to continue removing the data set?</html>";

	/**
	 * A generated serialization UID to keep compiler happy. 
	 */
	private static final long serialVersionUID = -4658096968134590995L;
}
