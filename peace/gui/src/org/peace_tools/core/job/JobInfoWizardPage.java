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

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as the second interactive page in a JobWizard.
 * This page permits the user to provide the basic information about
 * the data set to be used, provide a brief description of the job,
 * and decide if the job must perform clustering.
 * 
 * When the "Next >" button is clicked this wizard creates a blank
 * clustering data (at the job wizard level) for use by other pages
 * if the job must perform clustering. 
 */
public class JobInfoWizardPage extends GenericWizardPage {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: combo box to
	 * select data set, a text area for job description, and a check
	 * box for clustering. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public JobInfoWizardPage(JobWizard wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Job Information", 
				"Provide core information about the job");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the combo-box with list of data sets.
		dataSetList = new JComboBox();
		dataSetList.setBackground(Color.white);
		Utilities.adjustDimension(dataSetList, 500, 0);
				
		for(DataSet ds: Workspace.get().getDataSets()) {
			dataSetList.addItem(ds.toString());
		}
		// Create panel with the combo box for the user to choose
		JComponent dataSetBox = 
			Utilities.createLabeledComponents("Select Data Set for Job:",
					"(The EST file in data set will be processed)", 0, 
					dataSetList);
		
		// Create panels with description of the job
		description = new JTextArea(3, 10);
		JComponent descBox = 
			Utilities.createLabeledComponents("Description for job:",
					"(This is for your reference & can be anything)", 0, 
					new JScrollPane(description));
		// Create the informational label.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/16x16/Information.png"), 
				JLabel.LEFT);
		
		// Pack the input fields into a box
		JPanel subPanel = new JPanel();
		subPanel.setLayout(new BoxLayout(subPanel, BoxLayout.Y_AXIS));
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Pack the default fields into the subPanel.
		subPanel.add(info);
		subPanel.add(Box.createVerticalStrut(20));
		subPanel.add(dataSetBox);
		subPanel.add(Box.createVerticalStrut(20));
		subPanel.add(descBox);
		subPanel.add(Box.createVerticalStrut(20));
		// Finally add a filler to take up some vertical space.
		subPanel.add(Box.createVerticalGlue());
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the data being displayed in
	 * the GUI fields from the data stored in the in-memory
	 * Server object.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Nothing specific to do here.
	}

	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		// Save the necessary information prior to switching to 
		// the next page in the wizard?
		return true;
	}

	/**
	 * Obtain the data set that has been selected by the user.
	 * 
	 * This method is used by the JobWizard to add MST and Cluster
	 * data entries to the data set, if the wizard operations are
	 * successfully completed by the user.
	 * 
	 * @return The data set corresponding to the selection made by
	 * the user.
	 */
	protected DataSet getDataSet() {
		return Workspace.get().getDataSets().get(dataSetList.getSelectedIndex());
	}
	
	/**
	 * Obtain the description for the job.
	 * 
	 * @return The description entered by the user for this job.
	 */
	protected String getDescription() { return description.getText(); }
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final JobWizard wizard;
	
	/**
	 * A combo-box to select the data set to be used for this job.
	 * This field is created and populated in the constructor.
	 */
	private JComboBox dataSetList;

	/**
	 * Field to read and edit a brief description about the job.
	 * This information can be anything the user desires and is
	 * meaningful only to the user.
	 */
	private JTextArea description;

	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"<html>Select the data set that contains the EST file to be<br>" +
		"processed by this job. Subsequent wizard pages will permit<br>" +
		"setting up additional information for the job such as: MST<br>" +
		"generation parameters, server to run the job, etc.</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
