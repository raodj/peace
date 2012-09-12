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

package org.peace_tools.decagon.sdg;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.dataset.DataSetComboBox;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves typically as the second interactive page in
 * a EAST Job wizard. This page is used in one of the following
 * two different ways:
 * 
 * <ol>
 * 
 * <li>It is used by the EASTJobWizard to list the various data sets
 * (that contain valid MST and clustering file) that can be used for
 * assembly by EAST.</li>
 * 
 * <li>It can be used by the PEACEEASTJobWizard class to list the
 * various valid data sets that can be used for clustering followed
 * by assembly.</li>
 * 
 * </ol>
 * 
 * The choice of operation in the above two options is decided by 
 * the flag value that is passed in via the constructor.
 */
public class SourceGenesSelectionPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: combo box to
	 * select data set, a text area for job description, and a check
	 * box for clustering. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 */
	public SourceGenesSelectionPage(SyntheticDataSetGenerator wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Set Source Genes", 
				"Select the data set with source genes");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the combo-box with list of data sets (FASTA only) or MSTs
		dataSetList = new DataSetComboBox(true, true);
		dataSetList.setActionCommand("dataSetList");
		dataSetList.addActionListener(this);
		// Create panel with the combo box for the user to choose
		JComponent dataSetBox = 
			Utilities.createLabeledComponents("Select Data Set with source genes:",
					"The genes will be used to generate synthetic reads of different kinds", 0, false, 
					dataSetList);
	
		// Create the text box to enter path to target FASTA file.
		targetFilePath = new JTextField();
		browse         = Utilities.createButton(null, "Browse...", "browse", this, 
			"Select target FASTA file", true);
		JPanel horizBox = new JPanel(new BorderLayout(10, 0));
		horizBox.add(targetFilePath, BorderLayout.CENTER);
		horizBox.add(browse, BorderLayout.EAST);
		// Create a labeled component.
		JComponent dirBox =
			Utilities.createLabeledComponents("Enter target FASTA File Name:",
					"(If file exists it will be overwritten)",
					0, false, horizBox);

		// Create panel to permit user to enter description for the data set.
		description = new JTextArea(3, 4);
		JScrollPane jsp = new JScrollPane(description);
		jsp.setMinimumSize(description.getPreferredSize());
		JComponent descPanel = 
			Utilities.createLabeledComponents("Description for synthetic data set:",
					"(This is for your reference & can be anything)", 0, 
					false, jsp); 

		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
			dataSetBox, Box.createVerticalStrut(5),	dirBox, 
			Box.createVerticalStrut(5), descPanel);
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * Obtain the full path to the FASTA file that contains source genes
	 * to be used to generate synthetic reads.
	 * 
	 * @return The path to the FASTA file that contains source
	 * genes to be used to generate synthetic reads.
	 */
	protected String getSourceGeneFile() {
		Object selection = dataSetList.getSelection();
		if ((selection == null) || !(selection instanceof DataSet)) {
			// No valid selection? This is not right.
			return null;
		}
		// Extract the path for the files we care about.
		final String dataSetPath = ((DataSet) selection).getPath();
		return dataSetPath;
	}
	
	/**
	 * Obtain the path to the target FASTA file set by the user.
	 * 
	 * @return The path to the target FASTA file set by the user.
	 */
	protected String getTargetFASTAFile() {
		return targetFilePath.getText();
	}
	
	/**
	 * Obtain the description entered by the user for the 
	 * synthetic data set to be generated.
	 * 
	 * @return The description entered by the user.
	 */
	protected String getDescription() {
		return description.getText().trim();
	}
	
	/**
	 * Method to setup this page just before it is displayed.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before this page is displayed. This method updates the
	 * check box status based on the current selection in the
	 * combo box {@link #dataSetList}.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		int selIdx = Math.max(dataSetList.getSelectedIndex(), 0);
		if (dataSetList.getItemCount() > selIdx) {
			dataSetList.setSelectedIndex(selIdx);
		}
	}
	
	/**
	 * Method to ensure a valid data set is selected before proceeding to
	 * the next page.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before the user navigates to a different page. This method ensures
	 * that a valid data set has been selected for generating synthetic
	 * reads.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int prevPage) {
		if (currPage < prevPage) {
			return true;
		}
		// Before letting the user move on, ensure that the files
		// we are looking for are valid.
		Object selection = dataSetList.getSelection();
		if ((selection == null) || !(selection instanceof DataSet)) {
			// No valid selection? This is not right.
			return false;
		}
		// Extract the path for the files we care about.
		final String dataSetPath = ((DataSet) selection).getPath();
		// In all cases data set file must exist.
		File cDNAFile = new File(dataSetPath);
		if (!cDNAFile.exists() || !cDNAFile.canRead()) {
			return false;
		}
		// Check to ensure that the target file name is set
		final File targetFile = new File(targetFilePath.getText());
		if ((targetFile.getParentFile() != null) &&
			!targetFile.getParentFile().exists()) {
			JOptionPane.showMessageDialog(this, 
				"<html><i>The following directory specified for the<br/>" +
				"target FASTA file does not exist:</i><br/>" +
				Utilities.wrapStringToHTML(targetFile.getParentFile().getAbsolutePath(), 70) +
				"<br/><b>Please choose a different directory.</b></html>",
					"Invalid directory", JOptionPane.ERROR_MESSAGE);
			return false;
		} else if (targetFile.exists()) {
			final int choice = JOptionPane.showConfirmDialog(this, 
					"<html>The specified target FASTA file already exists.<br/>" +
					"The contents of this file will be overwritten.</br>" +
					"Do you want to proceed?</html>",
						"Invalid directory", JOptionPane.YES_NO_OPTION);
			return choice == JOptionPane.YES_OPTION;
		}
		// Proceed to next page.
		return true;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		final String cmd = e.getActionCommand();
		if ("browse".equals(cmd)) {
			JFileChooser jfc = new JFileChooser();
			jfc.setDialogTitle("Choose target FASTA File");
			if (jfc.showDialog(this, " OK ") == JFileChooser.APPROVE_OPTION) {
				// Copy the chosen file as the target FASTA file
				String wsPath = jfc.getSelectedFile().getAbsolutePath();
				targetFilePath.setText(wsPath);
			}
		} else if ("dataSetList".equals(cmd)) {
			// User changed the data set selected. Suggest new target
			// file name for user's convenience.
			final File srcFile = new File(getSourceGeneFile());
			String srcFileName = srcFile.getName();
			// Remove trailing extension on the file name
			int dotPos = srcFileName.lastIndexOf('.');
			if (dotPos > 0) {
				srcFileName = srcFileName.substring(0, dotPos);
			}
			// Get a unique file with in the default workspace directory.
			final String destFileName = 
				Utilities.getUniqueFileName(Workspace.get().getDirectory(), 
						srcFileName, ".fasta");
			targetFilePath.setText(destFileName);
		}
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final SyntheticDataSetGenerator wizard;
	
	/**
	 * A combo-box to select the data set to be used for this job.
	 * This field is created and populated in the constructor.
	 */
	private final DataSetComboBox dataSetList;

	/**
	 * Field to read and edit a brief description about the data set
	 * being generated. This information can be anything the user 
	 * desires and is meaningful only to the user.
	 */
	private final JTextArea description;

	/**
	 * Field to read/display the target path where the final generated
	 * FASTA file is to be placed.
	 */
	private final JTextField targetFilePath;
	
	/**
	 * The browse button that the user can use to select the target
	 * FASTA file.
	 */
	private final JButton browse;

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8586405984974688128L;
}
