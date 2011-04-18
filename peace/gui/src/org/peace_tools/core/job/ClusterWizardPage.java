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
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.FWAnalyzer;
import org.peace_tools.workspace.FileEntry;

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
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
				createClusterFileBox(), 
				Box.createVerticalStrut(10),
				createOthers());
		// Separate components from border a bit.
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
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
		clusterFile = new JTextField(10);
		Box horizBox = Box.createHorizontalBox();
		horizBox.add(clusterFile);
		horizBox.add(Box.createHorizontalStrut(10));
		horizBox.add(browse);
		// Create a labeled component.
		JComponent dirBox =
			Utilities.createLabeledComponents("Specify local cluster file:",
					"(The *.clr file is to be created and must not exist)",
					0, false, horizBox);		
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
	 * @return The sub-panel to which the controls are added. 
	 */
	private JPanel createOthers() {
		// Create panels with description, install folder, and
		// polling time.
		description = new JTextArea(3, 3);
		JScrollPane jsp = new JScrollPane(description);
		jsp.setMinimumSize(description.getPreferredSize());
		JComponent descBox = 
			Utilities.createLabeledComponents("Description for cluster file:",
					"(This is for your reference & can be anything)", 0, false,
					jsp);
		// Create the spinner.
		threshold = new JSpinner(new SpinnerNumberModel(130, -1, 1024, 10));
		JComponent threshBox = 
			Utilities.createLabeledComponents("Threshold for clustering:",
					"(See note below for tips)", 4, false, threshold);
		JPanel notePanel = createNotesPanel(); 
		notePanel.setAlignmentX(0);
		// Organize components in a vertical panel
		return Utilities.createLabeledComponents(null, null, 0, false, 
				descBox, Box.createVerticalStrut(10), 
				threshBox, Box.createVerticalStrut(5),
				notePanel);
	}
	
	/**
	 * Helper method to create a panel with multiple informational notes
	 * in them.
	 * 
	 * This method is invoked from the {@link #createOthers()} method to create
	 * a single panel with multiple, mutually exclusive informational notes
	 * (or JLabels) in them. Only one note is visible at any given time. The
	 * choice of the visible note is made based on the analyzer being currently
	 * used. The appropriate note is made visible when this wizard page is 
	 * displayed by the {@link #pageChanged(WizardDialog, int, int)} method.
	 * 
	 * @return A JPanel containing the various optional notes to be displayed
	 * to the user.
	 */
	private JPanel createNotesPanel() {
		JPanel notePanel = new JPanel(new BorderLayout(0, 0));
		// Create the default note and add it to the panel.
		infoLabels[0] = new JLabel(THRESHOLD_INFO, 
				Utilities.getIcon("images/24x24/Warning.png"), JLabel.LEFT);
		Utilities.adjustFont(infoLabels[0], -2, 8, -1);
		notePanel.add(infoLabels[0], BorderLayout.NORTH);
		// Create the two pass d2 note but make it invisible for now.
		infoLabels[1] = new JLabel(TWO_PASS_INFO, 
				Utilities.getIcon("images/24x24/Warning.png"), JLabel.LEFT);
		infoLabels[1].setVisible(false);
		Utilities.adjustFont(infoLabels[1], -2, 8, -1);
		notePanel.add(infoLabels[1], BorderLayout.CENTER);
		// Create the CLU note but make it invisible for now.
		infoLabels[2] = new JLabel(CLU_INFO, 
				Utilities.getIcon("images/24x24/Warning.png"), JLabel.LEFT);
		infoLabels[2].setVisible(false);
		Utilities.adjustFont(infoLabels[2], -2, 8, -1);
		notePanel.add(infoLabels[2], BorderLayout.SOUTH);
		// Return the panel with all three notes back to the caller.
		return notePanel;
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
		// Next setup the note to be displayed appropriately based on the
		// analyzer currently being used.
		// First hide all the notes.
		infoLabels[0].setVisible(false); infoLabels[1].setVisible(false);
		infoLabels[2].setVisible(false);
		threshold.setEnabled(false);
		// Now set the appropriate label to be visible.
		if (wizard.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2)) {
			infoLabels[1].setVisible(true);
			threshold.setValue(new Integer(1));
		} else if (wizard.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.CLU)) {
			infoLabels[2].setVisible(true);
			threshold.setValue(new Integer(-1));
		} else {
			// The default message regarding setting threshold value.
			infoLabels[0].setVisible(true);
			threshold.setValue(130);
			threshold.setEnabled(true);
		}
		// Relayout to handle changes in labels
		validate();
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
	 * Helper method to create a file entry for generated cluster file.
	 * 
	 * This method is used by the JobWizard to create a complete
	 * file entry for the generated cluster file.
	 * 
	 * @param id The workspace-wide unique ID to be set for this cluster
	 * file entry.
	 * 
	 * @return The file entry for the cluster file to be generated
	 * by this file.
	 */
	protected FileEntry getClusterFileEntry(String id) {
		FileEntry fe = new FileEntry(id, FileEntry.FileEntryType.CLS, 
				DataFileType.TXT, clusterFile.getText(), 
				description.getText());
		return fe;
	}

	/**
	 * Get the threshold value set by the user.
	 * 
	 * This is a helper method that is used by the job wizard to obtain
	 * the threshold value set by the user.
	 * 
	 * @return The threshold value set by the user.
	 */
	protected int getThreshold() {
		return ((Number) threshold.getValue()).intValue();
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
	 * This array contains three informational labels containing the 
	 * information in {@link #THRESHOLD_INFO}, {@link #TWO_PASS_INFO}, and
	 * {@link #CLU_INFO} respectively in that order. The labels are
	 * created in {@link #createNotesPanel()} method and are 
	 * suitably displayed by the {@link #pageChanged(WizardDialog, int, int)}
	 * method.
	 */
	private JLabel[] infoLabels = new JLabel[3];
	
	/**
	 * A generic informational message that provides additional information
	 * on the significance of the threshold value.
	 */
	private static final String THRESHOLD_INFO = 
		"<html>This is an <b>important</b> value. Small values result in many<br>" +
		"small clusters. Large thresholds make few, large clusters.<br>" +
		"Optimal threshold value is important  for best results.<br>" +
		"See help for more details." +
		"</html>";
	
	/**
	 * A simple text message that is displayed to the user when the user 
	 * selects the two-passed D2 configuration. The text is present to ensure
	 * that the user is clear about the operations being performed. The text is
	 * used in the {@link #pageChanged(WizardDialog, int, int)} method.
	 */
	private static final String TWO_PASS_INFO = "<html>" +
		"Two Pass D2 analyzer uses normalized distances and thresholds.<br>" +
		"This value must be set to 1 (one) for Two Pass D2 (and it <br>" +
		"cannot be changed). See help for more details." +
		"</html>";
	
	/**
	 * A simple text message that is displayed to the user when the user 
	 * selects the CLU analyzer. The text is present to ensure
	 * that the user is clear about the operations being performed. The text is
	 * used in the {@link #pageChanged(WizardDialog, int, int)} method.
	 */
	private static final String CLU_INFO = "<html>" +
		"The CLU analyzer automatically computes threshold based on mean.<br>" +
		"and variance of similarity values. Thefore the threshold is set to<br>" +
		"-1 as a sentinel value. See help for more details" +
		"</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1193172876744546755L;
}
