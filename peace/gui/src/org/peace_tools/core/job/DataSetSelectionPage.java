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

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.dataset.DataSetComboBox;
import org.peace_tools.core.dataset.DataSetComboBox.Choice;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FileEntry;

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
 * <li>It can be used by the DECAGON Job Creator (DJC) wizard class 
 * to list the various valid data sets that can be used for 
 * verification of assemblers.</li>
 * 
 * </ol>
 * 
 * The choice of operation in the above two options is decided by 
 * the flag value that is passed in via the constructor.
 */
public class DataSetSelectionPage extends GenericWizardPage implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: combo box to
	 * select data set, a text area for job description, and a check
	 * box for clustering. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param dataSetFlag If this flag is true, then this page presents
	 * the data sets in the workspace as primary inputs for assembly (or
	 * clustering+assembly as the case may be). If this flag is false, then 
	 * MST files are presented.
	 */
	public DataSetSelectionPage(WizardDialog wizard, boolean dataSetFlag) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("DataSet to Assemble", 
				"Select the data set to be clustered/assembled");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the combo-box with list of data sets (FASTA only) or MSTs
		dataSetList = new DataSetComboBox(dataSetFlag, true);
		dataSetList.addActionListener(this);
		
		// Create panel with the combo box for the user to choose
		JComponent dataSetBox = 
			Utilities.createLabeledComponents("Select Data Set for Job:",
					"(The cDNA file in data set will be clustered/assembled)", 0, false, 
					dataSetList);
		
		// Create the quality file usage check box and associated information label.
		useQualData = new JCheckBox("Use quality data for assembly");
		Utilities.adjustFont(useQualData, 0, 8, 1);
		useQualDataInfo = new JLabel(QUAL_FILE_AVAILABLE, 
				Utilities.getIcon("images/32x32/Note.png"), 
				JLabel.LEFT);
		// Pack the check box and the associated informational label into a panel
		JPanel qualPanel = new JPanel(new BorderLayout(0, 2));
		qualPanel.add(useQualData, BorderLayout.NORTH);
		qualPanel.add(useQualDataInfo, BorderLayout.SOUTH);
		
		// Create the informational label.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/32x32/Information.png"), 
				JLabel.LEFT);
		
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
			info, Box.createVerticalStrut(10),
			dataSetBox, Box.createVerticalStrut(10),
			qualPanel);
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * Obtain the DataSet of MSTData object that has been selected 
	 * by the user.
	 * 
	 * This method is used by the EASTJobWizard to obtain the
	 * source entry selected by the user. 
	 * 
	 * @return The DataSet or MSTData object corresponding to the 
	 * selection made by the user. If the user does not have a 
	 * valid selection to make, then this method returns null.
	 */
	public Object getSelection() {
		return this.dataSetList.getSelection();
	}
	
	/**
	 * This method is called whenever the user selects a 
	 * data set or MST from the {@link #dataSetList} combo box.
	 * This method appropriately enables/disables the
	 * {@link #useQualData} check box and updates the
	 * informational label {@link #useQualDataInfo} (that is
	 * displayed below the check box).
	 * 
	 * @param arg0 The action event associated with a call.
	 * Currently this event is ignored (and can be null).
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		Choice currDS = (Choice) dataSetList.getSelectedItem();
		boolean haveQualFile = currDS.hasQualFile();
		// Enable/disable check box depending on whether we have
		// quality file or not
		useQualData.setEnabled(haveQualFile);
		useQualData.setSelected(haveQualFile);
		// Update the info-label to provide the user with additional
		// information
		useQualDataInfo.setText(haveQualFile ? QUAL_FILE_AVAILABLE : 
			QUAL_FILE_NOT_AVAILABLE);
		// Update tool tip based on selection?
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
	 * Method to setup this page just before it is displayed.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before this page is displayed. This method updates the
	 * check box status based on the current selection in the
	 * combo box {@link #dataSetList}.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int prevPage) {
		if (currPage < prevPage) {
			return true;
		}
		// Before letting the user move on, ensure that the files
		// we are looking for are valid.
		Object selection = this.getSelection();
		if (selection == null) {
			// No valid selection? This is not right.
			return false;
		}
		// Extract the path for the files we care about.
		String dataSetPath = null;
		String mstFilePath = null;
		if (selection instanceof DataSet) {
			dataSetPath = ((DataSet) selection).getPath();
		} else {
			FileEntry mst = ((FileEntry) selection);
			mstFilePath = mst.getPath();
			dataSetPath = mst.getGFL().getDataSet().getPath();
		}
		// In all cases data set file must exist.
		File cDNAFile = new File(dataSetPath);
		if (!cDNAFile.exists() || !cDNAFile.canRead()) {
			return false;
		}
		// Check MST file if set
		if (mstFilePath != null) {
			File mstFile = new File(mstFilePath);
			if (!mstFile.exists() || !mstFile.canRead()) {
				return false;
			}
		}
		return true;
	}

	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final WizardDialog wizard;
	
	/**
	 * A combo-box to select the data set to be used for this job.
	 * This field is created and populated in the constructor.
	 */
	private DataSetComboBox dataSetList;
	
	/**
	 * A check box that provides a Boolean value to indicate if the quality
	 * file associated with a given data set should be used for assembly.
	 * This check box is enabled and disabled (depending on whether the
	 * data set has a quality file associated with it) so that the user
	 * can select or deselect) the option.
	 */
	private JCheckBox useQualData;
	
	/**
	 * A simple HTML formatted string that provides some basic information
	 * to help the user about as to why quality data option usage is
	 * currently unavailable.
	 */
	private static final String QUAL_FILE_NOT_AVAILABLE =
		"<html><font size=\"-2\">" +
		"This option is currently disabled because a suitable quality<br/>" +
		"file is not available for this data set. You need to add a<br/>" +
		"suitable quality file to this data set to enable this option." +
		"</font></html>";

	/**
	 * A simple HTML formatted string that provides some basic information
	 * to inform the user  as to how the quality data is going to be 
	 * used purely for assembly purposes.
	 */
	private static final String QUAL_FILE_AVAILABLE =
		"<html><font size=\"-2\">" +
		"If the above box is checked, then the quality file associated<br/>" +
		"with this data set will be used during assembly to improve the<br/>" +
		"overall quality the contigs generated after assembly." +
		"</font></html>";

	/**
	 * A label that is used to display information to the user as to why 
	 * the {@link #useQualData} check box is enabled or disabled. The
	 * text in this label is changed between {@link #QUAL_FILE_AVAILABLE}
	 * and {@link #QUAL_FILE_NOT_AVAILABLE} depending on the currently
	 * selected Choice in the {@link #dataSetList}
	 */
	private JLabel useQualDataInfo;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user. This information is worded such that it can be
	 * used for an assembly or clustering job.
	 */
	private static final String INFO_MSG = 
		"<html>Select the data set that contains the cDNA file to be<br>" +
		"processed. Subsequent wizard pages will permit setting<br>" +
		"suitable parameters &amp; server configuration for processing.</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1437405859971218347L;
}
