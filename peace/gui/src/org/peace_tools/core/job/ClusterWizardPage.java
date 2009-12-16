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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.JobSummary;
import org.peace_tools.workspace.MSTClusterData;

/**
 * This class serves as an interactive page in a JobWizard.
 * This page permits the user to provide clustering information. This
 * information includes the local cluster file where the clustering
 * data is to be stored. The user is permitted to enter some generic
 * description for their reference. In addition, an optional threshold
 * value for clustering can be specified by the user.  
 * 
 * Note that this wizard page checks to ensure that the cluster file
 * does not exist, yet.
 */
public class ClusterWizardPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a text box and
	 * button to specify the cluster file (this file must not yet exist
	 * as we are going to generate it from the job), a text area for
	 * job description,  and a JSpinner for the threshold value.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public ClusterWizardPage(JobWizard wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Clustering Setup", 
				"Configure cluster file information");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Pack the input fields into a box
		JPanel subPanel = new JPanel();
		subPanel.setLayout(new BoxLayout(subPanel, BoxLayout.Y_AXIS));
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		subPanel.add(Box.createVerticalGlue());
		subPanel.add(createClusterFileBox());
		subPanel.add(Box.createVerticalStrut(10));
		// Create other controls
		createOthers(subPanel);
		// Finally add a filler to take up some vertical space.
		subPanel.add(Box.createVerticalGlue());
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * Helper method to create the cluster file entry dialog box along
	 * with a "Browse" button. This method is invoked only once
	 * from the constructor. This method was introduced just to
	 * streamline the code in the constructor
	 * 
	 * @return A component containing the cluster file entry text box 
	 * along with informational labels.
	 */
	private JComponent createClusterFileBox() {
		// Put the install path text file and browse button into a 
		// single horizontal panel.
		browse = Utilities.createButton(null, " Browse ", 
				"Browse", this, "Browse local file system to " +
				"choose cluster file folder", true);
		clusterFile = new JTextField(30);
		Utilities.adjustDimension(clusterFile, 200, 4);
		Box horizBox = Box.createHorizontalBox();
		horizBox.add(clusterFile);
		horizBox.add(Box.createHorizontalStrut(10));
		horizBox.add(browse);
		// Create a labeled component.
		JComponent dirBox =
			Utilities.createLabeledComponents("Specify local cluster file:",
					"(The *.clr file is to be created and must not exist)",
					0, horizBox);		
		// Return the box to the caller
		return dirBox;
	}
	
	/**
	 * Helper method to create the description and threshold controls.
	 * 
	 * This method is invoked only once from the constructor. This 
	 * method was introduced just to streamline the code in the 
	 * constructor
	 * 
	 * @param subPanel The sub-panel to which the controls are to be
	 * added. 
	 */
	private void createOthers(JPanel subPanel) {
		// Create panels with description, install folder, and
		// polling time.
		description = new JTextArea(3, 10);
		JComponent descBox = 
			Utilities.createLabeledComponents("Description for cluster file:",
					"(This is for your reference & can be anything)", 0, 
					new JScrollPane(description));
		subPanel.add(descBox);
		subPanel.add(Box.createVerticalStrut(10));
		// Create the spinner.
		threshold = new JSpinner(new SpinnerNumberModel(130, 1, 1024, 10));
		Utilities.adjustDimension(threshold, 100, 3);
		JComponent threshBox = 
			Utilities.createLabeledComponents("Threshold for clustering:",
					"(See note below for tips)", 0, threshold);
		subPanel.add(threshBox); 
		JLabel note = new JLabel(THRESHOLD_INFO, 
				Utilities.getIcon("images/24x24/Warning.png"), JLabel.LEFT);
		Utilities.adjustFont(note, -2, 8, -1);
		note.setAlignmentX(0);
		subPanel.add(Box.createVerticalStrut(5));
		subPanel.add(note);
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * This method provides the user with a default cluster file name
	 * based on the MST file name specified by the user.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Create a default cluster file name as a suggestion to the user
		String mstFileName = wizard.getMSTFileName();
		// Remove trailing ".mst" extension if present
		if (mstFileName.endsWith(".mst")) {
			mstFileName = mstFileName.substring(0, mstFileName.length() - 4);
		}
		// Compute cluster file name by simply tagging ".cls" to mst name
		String clusterFileName = mstFileName + ".cls";
		clusterFile.setText(clusterFileName);
	}
		

	/**
	 * Method to handle clicking of "Browse" button. This method essentially
	 * enables and disables the various inputs depending on the check box.
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		JFileChooser jfc = new JFileChooser();
		jfc.setDialogTitle("Choose target cluster file");
		jfc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		if (jfc.showDialog(this, "Use Directory") == JFileChooser.APPROVE_OPTION) {
			// Copy the chosen directory to the mst file entry
			String wsPath = jfc.getSelectedFile().getAbsolutePath();
			// Set the selected item in this path.
			clusterFile.setText(wsPath);
		}
	}
	
	/**
	 * This method validates the cluster file set for in this page.
	 * 
	 * This method is invoked just before the user navigates to the
	 * adjacent page. This method checks to ensure that the cluster
	 * file is valid before user navigates to the next page.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		// Verify cluster file is valid.
		File cluster = new File(clusterFile.getText());
		Exception exp = null;
		// Try to create the file temporarily to check
		try {
			cluster.createNewFile();
			// Delete temporary file.
			cluster.delete();
		} catch (IOException e) {
			// Problem with the file. Save and report.
			exp = e;
		}

		// Ensure that the MST file is valid.
		if (exp != null) {
			JPanel msg = Utilities.collapsedMessage(INVALID_CLR_MSG, 
					Utilities.toString(exp));
			JOptionPane.showMessageDialog(this, msg,
					"Invalid Cluster File", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		// Go a head.
		return true;
	}
	
	/**
	 * Helper method to create a MSTCluster data.
	 * 
	 * This method is used by the JobWizard to create a complete
	 * MSTCluster data entry.
	 * 
	 * @param id The workspace-wide unique ID to be set for this cluster
	 * file entry.
	 * 
	 * @param mstID The ID associated with the MST data file that was
	 * used to generate this cluster.
	 * 
	 * @param summary The summary of the job that is going to be used
	 * to create this cluster file.
	 * 
	 * @return The cluster file entry to be added to the workspace.
	 */
	protected MSTClusterData getClusterEntry(String id, String mstID, 
			JobSummary summary) {
		int thresh = ((Number) threshold.getValue()).intValue();
		return new MSTClusterData(id, mstID, clusterFile.getText(), "", 
				thresh, summary);
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
	private JTextField clusterFile;

	/**
	 * The browse button to be enabled the user to choose the
	 * directory where the file is to be stored.
	 */
	private JButton browse;
	
	/**
	 * The spinner for setting the threshold. 
	 */
	private JSpinner threshold;

	/**
	 * A brief description about the cluster file. This information
	 * is meaningful only to the user.
	 */
	private JTextArea description;
	
	/**
	 * A generic informational message that is displayed to the user
	 * if the selected MST file is invalid.
	 */
	private static final String INVALID_CLR_MSG = 
		"<html>The cluster file does not meet the following criteria:" +
		"<ul>" +
		"<li>Path to cluster file must be valid</li>" +
		"<li>The file must not exist</li>" +
		"<li>The file must be creatable</li>" +
		"</ul>" +
		"</html>";
	
	/**
	 * A generic informational message that provides additional information
	 * on the significance of the threshold value.
	 */
	private static final String THRESHOLD_INFO = 
		"<html>This is an <b>important</b> value. Small values discriminate " +
		"ESTs more making more clusters. Larger thresholds make larger clusters. " +
		"Optional threshold value is important  for best results. " +
		"See help for more details." +
		"</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1193172876744546755L;
}
