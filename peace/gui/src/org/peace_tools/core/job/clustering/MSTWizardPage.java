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

package org.peace_tools.core.job.clustering;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.job.ServerPanel;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as an interactive page in a JobWizard.
 * This page permits the user to provide the information about the
 * local MST file where the data is to be stored. This wizard page 
 * also permits the user to select a server on which the job must 
 * be run. In addition, the wizard page permits the user to specify
 * the number of nodes and CPUs to be used for a parallel job.
 * 
 * Note that this wizard page checks to ensure that the MST file
 * does not exist, yet.
 */
public class MSTWizardPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a text box and
	 * button to specify the MST file (this file must not yet exist
	 * as we are going to generate it from the job), a combo box to
	 * select the server to run the job (this server must have peace
	 * already installed on it), JSpinners for number of nodes and
	 * CPUs per node to run the job. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public MSTWizardPage(ClusteringJobWizard wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Server Setup", 
				"Configure MST & Server configuration");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the server info panel that helps getting
		// server and job configuration from the user.
		serverInfo = new ServerPanel();
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
			createMSTFileBox(),	Box.createVerticalStrut(5),
			serverInfo.createServerPanel(true, true, null), 
			Box.createVerticalStrut(5));
		// Set up the border to make things look good.
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * Helper method to create the MST file entry dialog box along
	 * with a "Browse" button. This method is invoked only once
	 * from the constructor. This method was introduced just to
	 * streamline the code in the constructor
	 * 
	 * @return A component containing the MST file entry text box 
	 * along with informational labels.
	 */
	private JComponent createMSTFileBox() {
		// Put the install path text file and browse button into a 
		// single horizontal panel.
		browse = Utilities.createButton(null, " Browse ", 
				"Browse", this, "Browse local file system to " +
				"choose MST file folder", true);
		mstFile = new JTextField(30);
		Utilities.adjustDimension(mstFile, 0, 4);
		JPanel horizBox = new JPanel(new BorderLayout(10, 0));
		horizBox.add(mstFile, BorderLayout.WEST);
		horizBox.add(browse, BorderLayout.EAST);
		// Create a labeled component.
		JComponent dirBox =
			Utilities.createLabeledComponents("Specify local MST file:",
					"(The *.mst file is to be created and must not exist)",
					0, false, horizBox);		
		// Return the box to the caller
		return dirBox;
	}
	
	/**
	 * Method to handle clicking of "Browse" button. This method essentially
	 * enables and disables the various inputs depending on the check box.
	 * 
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		JFileChooser jfc = new JFileChooser();
		jfc.setDialogTitle("Choose target MST File directory");
		jfc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		if (jfc.showDialog(this, "Use Directory") == JFileChooser.APPROVE_OPTION) {
			// Copy the chosen directory to the mst file entry
			String wsPath = jfc.getSelectedFile().getAbsolutePath();
			// Set the selected item in this path.
			mstFile.setText(wsPath);
		}
	}

	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the list of server entries
	 * displayed in the combo-box. It also suggests a default MST file
	 * name to the user. However, the user can edit the default suggestion
	 * and set it to an appropriate value. In addition, it updates the memory
	 * requirement for running the job. 
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Let the server info update server list and associated information
		serverInfo.updateServerList(true, false, wizard.getDataSet());
		// Suggest a default MST file to the user based on EST file name
		DataSet ds = wizard.getDataSet();
		File  estFile = new File(ds.getPath());
		String mstFileName = estFile.getName();
		
		// Remove trailing extension on the file name
		int dotPos = mstFileName.lastIndexOf('.');
		if (dotPos > 0) {
			mstFileName = mstFileName.substring(0, dotPos);
		}
		final long uid = Workspace.get().getJobList().getSeqCounter();
		mstFileName += "_" + uid;
		mstFileName += ".mst";
		// Make the path absolute with respect to workspace
		mstFileName = Workspace.get().getDirectory() + 
			File.separator + mstFileName;
		mstFile.setText(mstFileName);
	}
	
	/**
	 * This method validates the MST file set for in this page.
	 * 
	 * This method is invoked just before the user navigates to the
	 * adjacent page. This method checks to ensure that the EST file
	 * is valid before user navigates to the next page. In addition,
	 * it generates warnings about the server, if the current status
	 * of the server is not operational. 
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		// Check to ensure we have a valid server selection.
		if (serverInfo.getSelectedServer() == null) {
			return false;
		}
		// Verify mst file is valid.
		File mst = new File(mstFile.getText());
		Exception exp = null;
		// Try to create the file temporarily to check
		try {
			mst.createNewFile();
			// Delete temporary file.
			mst.delete();
		} catch (IOException e) {
			// Problem with the file. Save and report.
			exp = e;
		}

		// Ensure that the MST file is valid.
		if (exp != null) {
			JPanel msg = Utilities.collapsedMessage(INVALID_MST_MSG, 
					Utilities.toString(exp));
			JOptionPane.showMessageDialog(this, msg,
					"Invalid MST File", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		// Verify the server status is still valid?
		return true;
	}
	
	/**
	 * Obtain the path to the MST file entered by the user.
	 * 
	 * @return The path to the MST file entered by the user.
	 */
	protected String getMSTFile() {
		File file = new File(mstFile.getText());
		return file.getAbsolutePath();
	}
	
	/**
	 * Obtain the panel that contains the GUI components and data 
	 * regarding the server and its configuration.
	 * 
	 * @return The server panel associated with this wizard page.
	 */
	protected ServerPanel getServerInfoPanel() {
		return this.serverInfo;
	}

	/**
	 * Helper method to create a file entry for generated MST file.
	 * 
	 * This method is used by the JobWizard to create a complete
	 * file entry for the generated MST file.
	 * 
	 * @param id The workspace-wide unique ID to be set for this MST
	 * file entry.
	 * 
	 * @param description THe description to be set for the
	 * file entry added by this method.
	 * 
	 * @return The file entry for the cluster file to be generated
	 * by this file.
	 */
	protected FileEntry getMSTFileEntry(String id, String description) {
		FileEntry fe = new FileEntry(id, FileEntry.FileEntryType.MST, 
				DataFileType.TXT, mstFile.getText(), 
				description);
		return fe;
	}
	
	/**
	 * Obtain information about the server and job that the user has 
	 * configured.
	 *  
	 * @param indent A simple indent string for indenting the data
	 * displayed by this server.
	 * 
	 * @return A multi-line string containing various information about
	 * the job including: name of server, nodes, cpus/node, memory/node,
	 * and estimated run time.
	 */
	protected String getSummary(String indent) {
		// Obtain the ID of the server from the selected list.
		Server server   = serverInfo.getSelectedServer();
		// Compute total memory for the job and display it for convenience
		final int[] config = serverInfo.getPlatformConfiguration();
		// Compute total memory as: #Nodes * #CPUsPerNode * MemPerNode
		final int   totalMemory = config[0] * config[1] * config[2];
		return "Server Information:\n" +
			indent + "Server: " + server.getName() + "\n" +
			indent + "Nodes: " + config[0] + "\n" +
			indent + "CPUs per Node: " + config[1] + "\n" +
			indent + "Memory per Node (MB): " + config[2] +  "\n" +
			indent + "Total memory (MB)   : " + totalMemory +  "\n" + 
			indent + "Esti. run time (hours): " + config[3];
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final ClusteringJobWizard wizard;

	/**
	 * Field to read/display the install path where the MST file
	 * will be saved after it is generated on the server.
	 */
	private JTextField mstFile;

	/**
	 * The browse button to be enabled the user to choose the
	 * directory where the file is to be stored.
	 */
	private JButton browse;

	/**
	 * A generic informational message that is displayed to the user
	 * if the selected MST file is invalid.
	 */
	private static final String INVALID_MST_MSG = 
		"<html>The EST file does not meet the following criteria:" +
		"<ul>" +
		"<li>Path to MST file must be valid</li>" +
		"<li>The file must not exist</li>" +
		"<li>The file must be creatable</li>" +
		"</ul>" +
		"</html>";
	
	/**
	 * A helper panel class that provides the necessary GUI features
	 * to interactively obtain server choice and job configuration
	 * information from the user.
	 */
	private final ServerPanel serverInfo;
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1193172876744546755L;
}
