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

package org.peace_tools.core.dataset;

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.MainFrame;
import org.peace_tools.data.ESTList;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as the top-level class for adding a new 
 * data set entry to the work space. This top-level class merely
 * creates the various pages and adds them to the wizard. Each
 * page performs a specific task required to create or edit
 *  a complete data set.
 */
public class DataSetWizard extends WizardDialog {
	/**
	 * The constructor for the data set wizard.
	 * 
	 * The constructor lays out the wizard and creates all the
	 * wizard pages. The wizard pages are created as blank pages
	 * and get populated with the necessary information just before
	 * the pages get displayed to the user.
	 * 
	 * @param title The title to be set for the main wizard frame.
	 * @param parent The main frame that logically owns the dialog
	 */
	public DataSetWizard(String title, MainFrame parent) {
		super(parent);
		setTitle(title);
		setResizable(false);
		// Set up the title image we want to use.
		setTitleBackground("images/peace_wizard_header.png", Color.white);
		// Set up the column image we want to use.
		setSequenceBackground("images/peace_wizard_column.png");
		// Create the data set that this wizard is going to be
		// modifying.
		dataSet = new DataSet("", null, "", DataSet.DataFileType.FASTA);
		// First setup the overview page.
		createOverview();
		// Create page to read/edit EST FASTA File
		ESTInfoWizardPage eiwp = new ESTInfoWizardPage(this, dataSet, false);
		addPage(eiwp);
		// Create the final summary page.
		SummaryWizardPage swp = new SummaryWizardPage(this, dataSet);
		addPage(swp);
	}

	/**
	 * Helper method invoked when user clicks cancel button.
	 * 
	 * This is a helper method that is overridden in this class.
	 * This method is invoked when the user clicks the cancel
	 * button in the wizard. This method is used to display a 
	 * confirmation dialog to ensure that the user really wants to
	 * exit out of the wizard.
	 * 
	 * @return This method returns true if the user wants to quit
	 * out of the wizard dialog.
	 */
	@Override
	protected boolean cancel() {
		// Check if user really want's to quit.
		int result = JOptionPane.showConfirmDialog(this,
				"Are you sure you want to exit from this wizard?",
				"Confirm", JOptionPane.YES_NO_OPTION);
		if (result == JOptionPane.NO_OPTION) {
			// The user does not want to quit.
			return false;
		}
		// Yes, the user wants to quit. Ensure no threads are left.
		return super.cancel();
	}
	
	/**
	 * Helper method to create the overview page. This method was
	 * introduced to keep the code clutter in the constructor to
	 * a bare minimum.
	 */
	private void createOverview() {
		JLabel message = new JLabel(OVERVIEW_MSG);
		Utilities.adjustFont(message, 0, 10, -1);
		GenericWizardPage overview = new GenericWizardPage();
		overview.add(message, BorderLayout.NORTH);
		overview.setTitle("Overview", "Overview of tasks in this wizard.");
		overview.setBorder(new EmptyBorder(20, 15, 10, 5));
		addPage(overview);
	}
	
	/**
	 * This method overrides the final notification method in this
	 * wizard to add the new data set (if wizard was successfully
	 * completed) to the workspace. 
	 */
	@Override
	public void done(boolean success) {
		if (!success) {
			return;
		}
		// Add the data set to the workspace. But first update ID
		dataSet.setID(Workspace.get().reserveID());
		Workspace.get().addDataSet(dataSet);
		// Save the workspace automatically.
		MainFrame mf = (MainFrame) getParent();
		mf.saveDelayedWorkspace();
	}
	
	/**
	 * This method can be used to obtain the ESTs associated 
	 * with this data set.
	 * 
	 * @return The list of ESTs associated with this data set.
	 */
	public synchronized ESTList getESTList() { return estList; }
	
	/**
	 * Set the list of ESTs associated with this data set.
	 * 
	 * @param estList The list of ESTs associated with this data set.
	 */
	public synchronized void setESTList(ESTList estList) { this.estList = estList; }
	
	/**
	 * This is the server entry that is being currently edited via 
	 * this wizard. This object is accessed by the various pages
	 * in this wizard and updated as the user navigates through
	 * the wizard pages.
	 */
	private DataSet dataSet;

	/**
	 * The list of ESTs associated with this data set. This object
	 * is set whenever one of the pages loads the EST data.
	 */
	private ESTList estList;
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user.
	 */
	private static final String OVERVIEW_MSG = "<html>" +
		"This wizard guides you through the process of adding a data<br>" +
		"set to this workspace. A data set is a collection of a related<br>" +
		"set of data data files that provide the following information: <br>"+ 
		"<ul><li>A EST file (in FASTA format) containing EST data</li>" +
		"<li>MST file(s) generated by PEACE for the above EST file</li>" +
		"<li>Clustering file(s) generated from the above MST file</li></ul>" +
		
		"Once a basic data set containing just an EST sequence<br>" +
		"file (in FASTA format) is added to the workspace, the data set<br>" +
		"can be used to perform clustering (generates MST and clusters)<br>" +
		"on a server by scheduling a suitable job via other wizards.<br><br>" +
		
		"<b>Note</b>: Existing MST and cluster files can be added to the data<br>" +
		"set or can be generated after a basic data set containing just an<br>" +
		"EST data file (in FASTA format) has been added to the workspace.<br>" +
		"</html>";

	/**
	 * The generated serialization UID (need to keep the compiler happy) 
	 */
	private static final long serialVersionUID = 804993573751301886L;
}
