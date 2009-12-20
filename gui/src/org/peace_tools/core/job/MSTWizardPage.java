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

package org.peace_tools.core.job;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.ServerList;
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
	public MSTWizardPage(JobWizard wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Server Setup", 
				"Configure MST & Server configuration");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
			createMSTFileBox(),	Box.createVerticalStrut(5),
			createServerPanel(), Box.createVerticalStrut(5));
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
	 * Helper method to create the server information entry components
	 * in this wizard page. This includes a combo-box to select a
	 * server and spinners for nodes and cpus/node.
	 * 
	 * This method is invoked only once from the constructor. This 
	 * method was introduced just to streamline the code in the constructor
	 * 
	 * @return A component containing the input fields related to
	 * server information.
	 */
	private JPanel createServerPanel() {
		// Create and add u/v heuristic parameters.
		Box horizBox = Box.createHorizontalBox();
		horizBox.setAlignmentX(0);
		// Create the spinners for nodes and cpus/node
		nodeInfo[0] = new JSpinner(new SpinnerNumberModel(10, 1, 1024, 2));
		nodeInfo[1] = new JSpinner(new SpinnerNumberModel(2, 1, 32, 1));
		// Memory and wall clock time.
		nodeInfo[2] = new JSpinner(new SpinnerNumberModel(2048, 64, 102400, 256));
		nodeInfo[3] = new JSpinner(new SpinnerNumberModel(6, 1, 168, 1));
		
		// Place them in a 2x2 grid pattern.
		final String Labels[] = {"Compute Nodes:", "CPUs per Nodes:", 
								 "Max Memory (in MB):", "Max Run time (hours):"};
		JPanel grid = new JPanel(new GridLayout(2, 2, 10, 3));
		grid.setAlignmentX(0);
		for(int i = 0; (i < nodeInfo.length); i++) {
			Utilities.adjustDimension(nodeInfo[i], 10, 4); // Adjust size to look right
			grid.add(Utilities.createLabeledComponents(Labels[i], null, 0, false, nodeInfo[i]));
		}
		// Now put all the information into a nice titled panel.
		JPanel bigBox = Utilities.createLabeledComponents(null, null, 0, true,
			// First add server selection combo-box.
			createServerList(),
			// Second add the cpu, nodes/cpu spinners
			Box.createVerticalStrut(5),
			grid,
			// Add information at bottom.
			new JLabel(CPU_INFO_MSG, 
						Utilities.getIcon("images/16x16/Information.png"), JLabel.LEFT)
		);
		// Set border to to make things look good.
 		bigBox.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createTitledBorder("Server Info:"),
				BorderFactory.createEmptyBorder(5, 5, 5, 5)));
 		// Return the panel to the caller
		return bigBox;
	}
	
	/**
	 * Helper method to create the combo-box to select a
	 * server.
	 * 
	 * This method is invoked only once from the createServerPanel
	 * method. This method was introduced just to streamline the code 
	 * better.
	 * 
	 * @return A panel containing the server list and labels.
	 */
	private JPanel createServerList() {
		// Create and setup visual properties of the combo-box
		serverList = new JComboBox();
		serverList.setBackground(Color.white);
		// Pack the server list with a suitable label
		return 
			Utilities.createLabeledComponents("Select Server to Use for Job:",
					"(A serial or parallel job will be run on the server)", 0, 
					false, serverList);
	}

	/**
	 * Method to handle clicking of "Browse" button. This method essentially
	 * enables and disables the various inputs depending on the check box.
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
	 * displayed in the combo box.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Clear any previous entries
		serverList.removeAllItems();
		// Populate the combo-box with server entries from work space.
		ServerList list = Workspace.get().getServerList();
		ArrayList<Server> servers = list.getServers();
		// On add servers in good operational condition.
		for(Server srvr: servers) {
			if (!Server.ServerStatusType.GOOD.equals(srvr.getStatus())) {
				// Status is not good.
				continue;
			}
			// Create string with necessary information.
			String srvrInfo = srvr.getName() + " (ID: " + srvr.getID() + ")";
			serverList.addItem(srvrInfo);
		}
		// Suggest a default MST file to the user based on EST file name
		DataSet ds = wizard.getDataSet();
		File  estFile = new File(ds.getPath());
		String mstFileName = estFile.getName();
		
		// Remove trailing extension on the file name
		int dotPos = mstFileName.lastIndexOf('.');
		if (dotPos > 0) {
			mstFileName = mstFileName.substring(0, dotPos);
		}
		if (ds.getMSTList().size() > 0) {
			mstFileName += "_" + ds.getMSTList().size();
		}
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
	 * Helper method to obtain selected server ID.
	 * 
	 * This is a helper method that is used to look up the server
	 * based on its server ID string (rather than index position).
	 * The serverID string is obtained from the string representation
	 * of the server entry (assuming that the entry is in the form:
	 * serverName (ID: xxxx). 
	 * 
	 * @return The entry corresponding to the given selected server.
	 */
	protected Server getSelectedServer() {
		// Obtain entry in the form serverName (ID: xxxx)
		String entry = (String) serverList.getSelectedItem();
		// Extract the actual server ID from the entry
		String srvrID = entry.substring(entry.indexOf("ID: ") + 4);
		// Drop trailing ")"
		srvrID = srvrID.substring(0, srvrID.length() - 1);
		// Now use use the serverID to look up the entry in workspace
		ServerList list = Workspace.get().getServerList();
		return list.getServer(srvrID);
	}
	
	/**
	 * Obtain the name of the server that the user has selected.
	 *  
	 * @param indent A simple indent string for indenting the data
	 * displayed by this server.
	 * 
	 * @return This method is used to obtain the name of the server
	 * that the user has selected to run the job.
	 */
	protected String getSummary(String indent) {
		// Obtain the ID of the server from the selected list.
		Server server   = getSelectedServer();
		return "Server Information:\n" +
			indent + "Server: " + server.getName() + "\n" +
			indent + "Nodes: " + nodeInfo[0].getValue() + "\n" +
			indent + "CPUs per Node: " + nodeInfo[1].getValue() + "\n" +
			indent + "Max Memory (GB): " + nodeInfo[2].getValue() +  "\n" + 
			indent + "Max run time (hours): " + nodeInfo[3].getValue();
	}
	
	/**
	 * Helper method to create a job. This method creates a new Job
	 * entry and populates it with as much information as possible.
	 * 
	 * @param description The description for this job that was
	 * entered by the user in the JobInfoWizardPage.
	 * 
	 * @return A new job entry for the newly created job.
	 */
	protected Job createJobEntry(String description) {
		// Determine ID of the selected server.
		Server server   = getSelectedServer();
		// Get the nodes and cpus/node..
		int nodes   = ((Number) nodeInfo[0].getValue()).intValue();
		int cpus    = ((Number) nodeInfo[1].getValue()).intValue();
		int memory  = ((Number) nodeInfo[2].getValue()).intValue();
		int runTime = ((Number) nodeInfo[3].getValue()).intValue();
		// Get an unique ID for this job.
		JobList   jobList   = Workspace.get().getJobList();
		// Now we have all the info to create the new job entry. The
		// only thing critical thing that is not yet finalized is
		// the path where the job files are stored on the remote 
		// server. That will happen soon.
		Job job = new Job(jobList.reserveJobID(), description, 
				server.getID(),	null, nodes, cpus, memory, runTime);
		return job;
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final JobWizard wizard;

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
	 * The combo box that permits the user to select the server
	 * to be used for running the job.
	 */
	private JComboBox serverList;
	
	/**
	 * The array of two configuration parameter values for the CPUs
	 * and nodes per CPU on the server. The last two spinners are
	 * used for memory and run time values. 
	 */
	private JSpinner nodeInfo[] = new JSpinner[4];
		
	/**
	 * A generic informational message that is displayed to the user
	 * to provide information about the u/v heuristic.
	 */
	private static final String CPU_INFO_MSG = 
		"<html>The CPUs selected and nodes per CPU must match<br>" +
		"the configuration of the server.</html>";

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
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1193172876744546755L;
}
