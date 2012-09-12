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
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.generic.WizardPage;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Server;

/**
 * A GUI sub-component used in wizards to permit user to select an
 * appropriate server on which a job is to be run.
 * 
 * This component was primarily introduced to permit various wizards
 * to share common functionality of permitting a user to select
 * a server and specify configuration to be used for running a job.
 * 
 * This panel is not meant to serve directly as a wizard page but
 * rather as a part of a server page. This panel creates and organizes
 * its components depending on the requirements specified when the
 * panel is created. Various wizards use this component (as part of
 * their wizard pages) to display choice of servers.
 */
public class ServerPanel implements ChangeListener {
	/** The default and only constructor.
	 * 
	 * This constructor is very straightforward and sets various instance
	 * variables to their default initial value.
	 */
	public ServerPanel() {
		dataSet      = null;
		serverList   = null;
	}
	
	/**
	 * Helper method to create the server information entry components
	 * in this wizard page. This includes a ServerComboBox to select a
	 * server and spinners for nodes and cpus/node.
	 * 
	 * This method is invoked only once from the constructor. This 
	 * method was introduced just to streamline the code in the constructor
	 * 
	 * @param parallelJob If this parameter is true, then the spinners
	 * corresponding to number of nodes and CPUs/Node are enabled. Otherwise
	 * these spinners are disabled.
	 * 
	 * @param setTitle If this parameter is true, then a titled border is
	 * added around the panel created by this method.
	 * 
	 * @param subTitle If a sub-title string is specified (that is this 
	 * parameter is not null) then the sub-title is set around the cluster
	 * configuration options.
	 * 
	 * @return A component containing the input fields related to
	 * server information.
	 */
	public JPanel createServerPanel(final boolean parallelJob, final boolean setTitle,
			final String subTitle) {
		// Create spinners for job configuration.
		JPanel configPanel = new JPanel(new BorderLayout(0, 5));
		configPanel.add(createClusterConfigOptions(parallelJob), BorderLayout.NORTH);
		configPanel.add(new JLabel(parallelJob ? PARALLEL_CPU_INFO_MSG : SERIAL_CPU_INFO_MSG, 
						Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT),
						BorderLayout.SOUTH);
		// Setup a sub-title if required
		if (subTitle != null) {
			configPanel.setBorder(BorderFactory.createCompoundBorder(
					BorderFactory.createTitledBorder(subTitle),
					BorderFactory.createEmptyBorder(1, 5, 2, 5)));
		}
		// Now put all the information into a nice panel.
		JPanel bigBox = new JPanel(new BorderLayout(0, 7));
		bigBox.add(createServerList(), BorderLayout.NORTH);
		bigBox.add(configPanel, BorderLayout.SOUTH);
		/*JPanel bigBox = Utilities.createLabeledComponents(null, null, 0, true,
			// First add server selection combo-box.
			createServerList(),
			// Second add the cpu, nodes/cpu spinners
			Box.createVerticalStrut(5),
			configPanel
		);*/
		
		if (setTitle) {
			// Set border to to make things look good.
			bigBox.setBorder(BorderFactory.createCompoundBorder(
					BorderFactory.createTitledBorder("Server Info:"),
					BorderFactory.createEmptyBorder(1, 5, 5, 5)));
		}
 		// Return the panel to the caller
		return bigBox;
	}
	
	/**
	 * Helper method to create a 4x4 grid with inputs to permit the user to
	 * select job configuration for a cluster. This method creates 
	 * various spinners in the {@link #nodeInfo} array that permits the
	 * user to set the following job parameters: number of nodes, CPUs 
	 * per node, total memory, and estimated runtime. This information
	 * is used when a multi-processor job is created. If a single processor
	 * job is to be created then the compute nodes and cpus per node option
	 * is simply disabled.
	 *  
	 * @param parallelJob If this parameter is true, then the spinners
	 * corresponding to number of nodes and CPUs/Node are enabled. Otherwise
	 * these spinners are disabled.
	 * 
	 * @return This method returns a panel (organized as a 4x4 grid) that
	 * contains the four spinners that permit the user to configure 
	 * additional job parameters.
	 */
	private JPanel createClusterConfigOptions(final boolean parallelJob) {
		// Save parallel or serial configuration option
		parallelConfig = parallelJob;
		// Create the spinners for nodes and cpus/node
		nodeInfo[0] = new JSpinner(new SpinnerNumberModel(1, 1, 1024, 1));
		nodeInfo[1] = new JSpinner(new SpinnerNumberModel(1, 1, 32, 1));
		// Enable/disable these spinners depending on parallel or
		// sequential jobs.
		nodeInfo[0].setEnabled(parallelJob);
		nodeInfo[1].setEnabled(parallelJob);
		// Add listeners for the first two to update memory
		nodeInfo[0].addChangeListener(this);
		nodeInfo[1].addChangeListener(this);
		// Memory and wall clock time.
		nodeInfo[2] = new JSpinner(new SpinnerNumberModel(2048, 64, 102400, 16));
		nodeInfo[3] = new JSpinner(new SpinnerNumberModel(6, 1, 168, 1));
		
		// Place them in a 2x2 grid pattern.
		final String Labels[] = {"Compute Nodes:", "CPUs per Nodes:", 
								 "Memory per Node (in MB):", "Esti. Run time (hours):"};
		JPanel grid = new JPanel(new GridLayout(2, 2, 10, 3));
		grid.setAlignmentX(0);
		for(int i = 0; (i < nodeInfo.length); i++) {
			Utilities.adjustDimension(nodeInfo[i], 10, 4); // Adjust size to look right
			grid.add(Utilities.createLabeledComponents(Labels[i], null, 0, false, nodeInfo[i]));
		}
		return grid;
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
		serverList = new ServerComboBox();
		// A label to be displayed only in case where we don't have any
		// servers to operate with.
		emptyServerListLabel = new JLabel(EMPTY_SERVER_LIST, 
				Utilities.getIcon("images/24x24/Warning.png"), JLabel.LEFT);
		emptyServerListLabel.setVisible(false);
		emptyServerListLabel.setBackground(new Color(0xe0, 0xe0, 0xff)); // pale blue
		emptyServerListLabel.setOpaque(true);
		emptyServerListLabel.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createEtchedBorder(), BorderFactory.createEmptyBorder(1, 5, 1, 5)));
		// A label to display with the server entry when the server
		// entry to be used is locked.
		lockedServerLabel = new JLabel();
		lockedServerLabel.setBorder(BorderFactory.createEtchedBorder());
		lockedServerLabel.setVisible(false);
		// Pack the server list with a suitable label
		serverTitle = new JLabel("Select the server to use for the job:");
		return 
			Utilities.createLabeledComponents(null, null, 0, 
					false, serverTitle, serverList, emptyServerListLabel, 
					lockedServerLabel);
	}
		
	/**
	 * This method is called just before the panel is to be displayed to
	 * update the list of valid servers to be presented to the user.
	 * 
	 * This method is typically invoked from the wizard page (that logically
	 * owns this ServerPanel object) from the 
	 * {@link WizardPage#pageChanged(WizardDialog, int, int)} method.
	 *   
	 * This page essentially updates the list of server entries
	 * displayed in the combo-box to meet the requirements of the wizard page
	 * (based on the values of the parameters). In addition, it updates the memory
	 * requirement for running the job (if required).
	 * 
	 * @param needPEACE If this flag is set to true, then only server entries
	 * that have a usable PEACE install are presented to the user as valid choices.
	 * 
	 * @param needEAST If this flag is set to true, then only server entries
	 * that have a usable EAST install are presented to the user as valid choices.
	 * 
	 * @param dataSet The data set to be used by this panel to suggest memory and
	 * time needs for the job. If this parameter is null, then this panel will not
	 * provide default suggestions for these two parameters.
	 */
	public void updateServerList(final boolean needPEACE, final boolean needEAST,
			DataSet dataSet) {
		// Get list of valid servers
		ServerComboBox.getServerList(needPEACE, needEAST, serverList);
		// Update GUI elements appropriately
		checkServerListAndUpdateGUI();
		// Initialize the default memory requirement for running
		// this job via the listener method with a null event.
		stateChanged(null); // event can be null as it is not used.
	}

	/**
	 * Helper method that checks to see if the server list has at least one
	 * entry. Depending on whether the server list is empty or not it enables
	 * or disables various controls appropriately.
	 */
	private void checkServerListAndUpdateGUI() {
		// If even a single valid server entry was not found, then we change
		// overall status.
		final boolean haveServers = (serverList.getItemCount() > 0); 
		serverList.setVisible(haveServers);
		emptyServerListLabel.setVisible(!haveServers);
		// Enable/disable the config spinners as well
		nodeInfo[0].setEnabled(haveServers && parallelConfig);
		nodeInfo[1].setEnabled(haveServers && parallelConfig);
		nodeInfo[2].setEnabled(haveServers);
		nodeInfo[3].setEnabled(haveServers);
	}
	
	/**
	 * Method to set the server selection to a given server object.
	 * 
	 * This method can be used to set the current server selection
	 * (presented by this panel) and available options to a single
	 * server entry. 
	 * 
	 * @param srvr The only server to be presented to the user.
	 */
	public void lockSelectedServer(Server srvr) {
		serverList.removeAllItems();
		serverList.addItem(srvr);
		serverList.setSelectedItem(srvr);
		// Change the server title suitably.
		this.serverTitle.setText("The following pre-selected server will be used for this job:");
		// Update GUI elements appropriately
		checkServerListAndUpdateGUI();
		// Hide the job list and show the label.
		serverList.setVisible(false);
		lockedServerLabel.setVisible(true);
		// Copy the necessary information.
		final JLabel infoLabel = (JLabel) serverList.getListCellRendererComponent(new JList(), 
				srvr, 0, false, false);
		// Copy the necessary information to locked-server-label
		lockedServerLabel.setIcon(infoLabel.getIcon());
		lockedServerLabel.setText(infoLabel.getText());
		lockedServerLabel.setBackground(infoLabel.getBackground());
		lockedServerLabel.setForeground(infoLabel.getForeground());
		// Initialize the default memory requirement for running
		// this job via the listener method with a null event.
		stateChanged(null); // event can be null as it is not used.
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
	public Server getSelectedServer() {
		if (serverList.getSelectedIndex() == -1) {
			// No selections.
			return null;
		}
		// Obtain entry in the server object form.
		Server entry = (Server) serverList.getSelectedItem();
		return entry;
	}

	/**
	 * Obtain platform-specific job configuration information.
	 * 
	 * This method must be used to obtain the platform specific job configuration
	 * information.  This method converts the job information into an array and
	 * returns it. The returned array has the elements in the order: no. of. nodes,
	 * CPUs/node, total memory (in MB), and runTime (in hours).
	 * 
	 * @return An array of 4 integers that contains the job information entered
	 * by the user. The order of the elements in the array is fixed.
	 */
	public int[] getPlatformConfiguration() {
		int[] configInfo = new int[nodeInfo.length];
		for(int i = 0; (i < nodeInfo.length); i++) {
			configInfo[i] = ((Number) nodeInfo[i].getValue()).intValue();
		}
		return configInfo;
	}
	
	/**
	 * Listener for CPUs and Nodes/CPU spinner input boxes to update memory.
	 * 
	 * <p>This method is invoked whenever the user modifies the number of CPUs
	 * or nodes-per-CPU to be used for the job. This method updates the
	 * recommended amount of memory to be reserved for this job. The memory
	 * is computed using the formula:</p>
	 * 
	 * <p>Memory = (2 * FASTAFileSize) + 256MB</p>
	 * 
	 * <p>Note: Earlier the memory field represented the total aggregate memory
	 * to be allocated for the job. This value was directly passed on as part of
	 * PBS job script. However, this was confusing to users. So this field
	 * has now been revised to max memory per node.  The net memory is computed
	 * by multiplying this value with number of CPUs and number of nodes just before
	 * the job is submitted.</p> 
	 * 
	 * @param e The change event associated with this call back. Currently,
	 * this parameters is not used by this method.
	 */
	@Override
	public void stateChanged(ChangeEvent e) {
		if (dataSet == null) {
			return;
		}
		// Couple constants to make code more meaningful
		final int megaBytes   = 1024 * 1024;
		final int minMemory   = 256;  // Megabytes

		// Determine FASTA file size in megabytes
		final String fileName = dataSet.getPath();
		final File   tempFile = new File(fileName);
		// Ensure file resolves to at least 1 MB to ease computations
		final long size       = (tempFile.length() / megaBytes) + 1;
		
		// Compute minimum memory in MB.
		final long memorySize = size + size + minMemory;
		// Setup the suggested memory footprint
		nodeInfo[2].setValue(new Long(memorySize));
	}
	
	/**
	 * Convenience method to summarize information about the
	 * information entered in this server panel.
	 * 
	 * This method is typically used to display summary information
	 * about the data entered by the user in this panel.
	 * 
	 * @param sw The summary writer object to be used for writing
	 * summary information.
	 * 
	 * @param subSection If this flag is true then this method creates
	 * summary information as a sub-section. Otherwise the summary
	 * information is created as a section.
	 */
	public void summarize(SummaryWriter sw, boolean subSection) {
		if (subSection) {
			sw.addSubSection("Job Configuration", "", "");	
		} else {
			sw.addSection("Job Configuration");
		}
		// Add server information from this server panel.
		final Server srvr = getSelectedServer();
		srvr.summarize(sw, false, true);
		// Add information about the memory and other information.
		final int jobConfig[] = getPlatformConfiguration();
		summarize(sw, subSection, "Nodes", "" + jobConfig[0], 
				  "Number of nodes in cluster to reserve for the job");
		summarize(sw, subSection, "CPUs/Nodes", "" + jobConfig[1], 
				"Number of CPUs per node to reserve for the job");
		summarize(sw, subSection, "Max Memory", "" + jobConfig[3] + " MB", 
				 "Peak memory (in Megabytes) to request for the job");
		summarize(sw, subSection, "Run time", "" + jobConfig[1] + " hours", 
				  "Maximum runtime to request for the job");
	}
	
	/**
	 * Helper method to add a line of summary or sub-summary.
	 * 
	 * This is a helper method that is internally used (by
	 * the {@link #summarize(SummaryWriter, boolean)} method) to
	 * add an summary or sub-summary entry depending on a flag.
	 * 
	 * @param sw The summary writer object to be used to create
	 * a summary or sub-summary entry.
	 * 
	 * @param subSection If this flag is true then the given 
	 * information is used to create a sub-summary entry. Otherwise
	 * the given information is used to add a summary entry.
	 * 
	 * @param name The name of the parameter or property for which a summary
	 * line is being added. This parameter can be null.
	 * 
	 * @param value The value for the given parameter to be added.
	 * This parameter can be null.
	 * 
	 * @param desc An optional description to be associated with the summary
	 * line. This parameter can be null. 
	 */
	private void summarize(SummaryWriter sw, boolean subSection, 
			final String name, final String value, final String desc) {
		if (!subSection) {
			sw.addSummary(name, value, desc);
		} else {
			sw.addSubSummary(name, value, desc);
		}
	}
	
	/**
	 * The combo box that permits the user to select the server
	 * to be used for running the job.
	 */
	private ServerComboBox serverList;
	
	/**
	 * The array of two configuration parameter values for the CPUs
	 * and nodes per CPU on the server. The last two spinners are
	 * used for memory and run time values. 
	 */
	private JSpinner nodeInfo[] = new JSpinner[4];

	/**
	 * Flag to indicate if the current setup of the various spinners in the
	 * {@link #nodeInfo} array is setup for parallel (multi-node, multi-CPU)
	 * or serial (single-node, single-CPU) job. This flag is used to enable
	 * and disable spinners in the {@link #checkServerListAndUpdateGUI()}
	 * method. The values are set/initialized in the 
	 * {@link #createClusterConfigOptions(boolean)} method.
	 */
	private boolean parallelConfig;
	
	/**
	 * The data set that is going to be used in this job. This data set
	 * is used to suggest memory requirements (and possibly time later on)
	 * whenever the user changes the number of CPUs or Nodes/CPUs. 
	 */
	private DataSet dataSet;
	
	/**
	 * A generic informational message that is displayed to the user for
	 * parallel jobs. This information is displayed only for parallel
	 * and not for sequential jobs.
	 */
	private static final String PARALLEL_CPU_INFO_MSG = 
		"<html><font size=\"-2\">" + 
		"The nodes selected and CPUs-per-node must match the server<br>" +
		"configuration. Refer to user manual for setting up hostfile<br>" +
		"for jobs using openmpi/mpich.</font></html>";
	
	/**
	 * A generic informational message that is displayed to the user for
	 * sequential jobs (at the bottom of the panel).
	 */
	private static final String SERIAL_CPU_INFO_MSG = 
		"<html><font size=\"-2\">" + 
		"This is a serial job. Consequenlty setting nodes and CPUs-per-node<br/>" +
		"is disabled (they are enabled only for parallel jobs).</font></html>";

	/**
	 * This is a convenience label that is used to display a message (instead
	 * of the server selection list) when a suitable server was not found.
	 * This label displays the message in {@link #EMPTY_SERVER_LIST} string
	 * along with a nice looking icon. This label is created and setup
	 * in the {@link #createServerList()} method.
	 */
	private JLabel emptyServerListLabel;
	
	/**
	 * Label displayed on top of the server selection ({@link #serverList}) 
	 * combo-box. This label is changed if the server selection is locked
	 * to a specific entry via call to {@link #lockSelectedServer(Server)}
	 * method. The label is changed to provide the user with more
	 * meaningful message.
	 */
	private JLabel serverTitle;
	
	/**
	 * A simple label that is created to display the server information 
	 * when the server entry for this panel is locked.
	 */
	private JLabel lockedServerLabel;
	
	/**
	 * A simple informational message that is displayed (instead of the server
	 * selection list) when a suitable server (in good operating condition that
	 * has not had a connection issue recently) with the necessary software
	 * components (PEACE, EAST, or DECAGON) installed on it. This string is
	 * used in the {@link #updateServerList(boolean, boolean, DataSet)} method.
	 */
	private static final String EMPTY_SERVER_LIST =
		"<html>" +
		"Could not find a server in good operating status with the required<br/>" +
		"software components installed. Please create a suitable server entry.<br/>" +
		"</html>";	
}
