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

package org.peace_tools.core;

import java.awt.BorderLayout;
import java.awt.Frame;
import java.awt.GridLayout;
import java.io.File;
import java.io.FileOutputStream;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.MSTClusterData;
import org.peace_tools.workspace.MSTData;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * This is a dialog class for copying generated files.
 * 
 * This is a convenience dialog that is used to interactively copy 
 * files from a given server to a local file. This dialog is typically
 * launched whenever a job successfully completes.  
 */
public class FileCopyDialog extends JDialog implements Runnable {
	/**
	 * The constructor for creating and configuring this dialog.
	 * 
	 * The constructor organizes the various components in this dialog
	 * box and lays out the dialog in preparation for the copy operation.
	 * The copy operation is jump started as a background thread when the 
	 * dialog is shown. 
	 * 
	 * @param parent The parent window to which this wizard dialog logically
	 * belongs.
	 */
	public FileCopyDialog(Frame parent, Job job) {
		// Initialize the JDialog
		super(parent, "", true);
		// Save reference to the job
		this.job = job;
		// Obtain and save the MSTData and ClusterDataEntries that
		// contain the name of the destination file for this job.
		Workspace workspace    = Workspace.get();
		MSTData mstEntry       = workspace.getMSTData(job.getJobID());
		MSTClusterData cluster = workspace.getClusterData(job.getJobID());
		// Save the target file names in an array for later reference.
		targetFileNames    = new String[2];
		targetFileNames[0] = mstEntry.getPath();
		targetFileNames[1] = cluster.getPath();
		// Setup the labels using a helper method in a panel.
		fileInfo    = new JLabel[2];
		fileInfo[0] = createLabel("MST", targetFileNames[0]); 
		fileInfo[1] = createLabel("Cluster", targetFileNames[1]);
		// Create a panel to contain the labels.
		JPanel labelPanel = new JPanel(new GridLayout(2, 1));
		labelPanel.add(fileInfo[0]);
		labelPanel.add(fileInfo[1]);
		// Create a panel with information and labels 
		JPanel infoPanel = new JPanel(new BorderLayout(5, 0));
		infoPanel.add(new JLabel("<html><b>Please wait while the following File(s) are<br>" + 
			"being transferred into the workspace:</b></html>"), BorderLayout.NORTH);
		infoPanel.add(labelPanel, BorderLayout.SOUTH);
		
		// Now create the panel with the label and progress bar
		// to display while copying.
		JPanel progPanel = new JPanel(new BorderLayout(5, 2));
		progPanel.add((progressInfo = new JLabel(" ")), BorderLayout.NORTH);
		progPanel.add((progressBar  = new JProgressBar()), BorderLayout.SOUTH);
		// Set the font on the progress info to be a bit smaller
		Utilities.adjustFont(progressInfo, -2, 8, 1);
		
		// Finally pack all the components into a main panel.
		JPanel mainPanel = new JPanel(new BorderLayout(5, 5));
		mainPanel.setBorder(new EmptyBorder(10, 10, 10, 10));
		final String IconFile = "images/32x32/Download.png";
		mainPanel.add(new JLabel(Utilities.getIcon(IconFile)), BorderLayout.WEST);
		mainPanel.add(infoPanel, BorderLayout.CENTER);
		mainPanel.add(progPanel, BorderLayout.SOUTH);
		// Finally, add all the components and sub-panels
		setLayout(new BorderLayout(5, 5));
		add(mainPanel, BorderLayout.CENTER);
		setResizable(false);
	}
	
	/**
	 * Helper method to create a label with succinct file name.
	 * 
	 * @param type The type of the file to be displayed in the label.
	 * @param fileName The name of the target file to be copied.
	 * @return A JLabel containing the file information for display.
	 */
	private JLabel createLabel(String type, String fileName) {
		File temp    = new File(fileName);
		String msg   = type + " data file " + temp.getName();
		JLabel label = new JLabel(msg,
				Utilities.getIcon("images/16x16/Box.png"),
				JLabel.LEFT);
		label.setBorder(new EmptyBorder(5, 5, 5, 5));
		return label;
	}
	
	/**
	 * The main method that lays out the components in this wizard
	 * and displays the wizard. 
	 */
	public void showDialog() {
		// Pack the wizard forcing layout to happen
		pack();
		// Center the window on the frame.
		Utilities.centerPanel((Frame) getParent(), this);
		// Start the background copy thread.
		Thread copyThread = new Thread(this);
		copyThread.start();
		// Finally show the wizard.
		setVisible(true);
	}
	
	/**
	 * Thread method to do the actual copying operation.
	 * 
	 * This method is invoked from a separate background thread to
	 * perform the actual copy operation. This method uses the
	 * server session created below to copy the files.
	 */
	@Override
	public void run() {
		int i = 0;
		try {
			Thread.sleep(250);
			Server server = Workspace.get().getServerList().getServer(job.getServerID());
			session = SessionFactory.createSession(this, server);
			session.connect();
			// Copy the files, one file at a time.
			for(i = 0; (i < targetFileNames.length); i++) {
				// Set the progress bar into indeterminate mode.
				progressBar.setValue(0);
				progressBar.setIndeterminate(true);
				// Indicate that we are copying a given file.
				this.fileInfo[i].setIcon(Utilities.getIcon("images/16x16/Download.png"));
				// First create a output stream to copy data to
				File destFile = new File(targetFileNames[i]);
				progressInfo.setText("Transferring: " + destFile.getName());
				// Create an output stream.
				FileOutputStream fos = new FileOutputStream(targetFileNames[i]);
				session.copy(fos, job.getPath(), destFile.getName(), progressBar);
				// Close the destination file.
				fos.close();
				// Update the success status.
				this.fileInfo[i].setIcon(Utilities.getIcon("images/16x16/CheckedBox.png"));
			}
			// Set status of the job to success
			job.setStatus(JobBase.JobStatusType.SUCCESS);
		} catch (Exception e) {
			String errMsg = "<html>Error copying the following file to the GUI<br>" +
					"File: " + targetFileNames[i] + "<br>" +
					"<b>You may retry the copy operation again.<b></html>";
			JPanel msg = Utilities.collapsedMessage(errMsg, Utilities.toString(e));
			JOptionPane.showMessageDialog(this, msg, 
					"Error copying file", JOptionPane.ERROR_MESSAGE);
		} finally {
			// Ensure that the wizard shuts down.
			final JDialog dialog = this;
			SwingUtilities.invokeLater(new Runnable() {
				@Override
				public void run() {
					dialog.setVisible(false);
					dialog.dispose();
				}
			});
		}
	}
	
	/**
	 * The job entry for which the files are to be copied from the
	 * server (on which the job was run) to the local machine on
	 * which the GUI is running.
	 */
	private final Job job;
	
	/**
	 * These labels are created in the constructor for each file to
	 * be copied. The labels provide visual information about the
	 * files to the user.
	 */
	private JLabel[] fileInfo;
	
	/**
	 * Label displayed on top of the progress bar indicating the
	 * file being currently copied.
	 */
	private JLabel progressInfo;
	
	/**
	 * The progress bar that is displayed at the bottom of the 
	 * dialog box to indicate the progress made in copying
	 * files.
	 */
	private JProgressBar progressBar;
	
	/**
	 * The server session to be used for copying the necessary files.
	 * This session is created in the run() method once the dialog
	 * is displayed.
	 */
	private ServerSession session;
	
	/**
	 * The names of the files (with absolute path on local machine)
	 * to be copied from the remote machine to the local machine.
	 * This list is populated in the constructor.
	 */
	private String[] targetFileNames;
	
	/**
	 * A serialization constant (to keep the compiler happy)
	 */
	private static final long serialVersionUID = 5397447232688023756L;
}
