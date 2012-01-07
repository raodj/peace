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

package org.peace_tools.core.job.baton;

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.MainFrame;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;


/**
 * This class serves as the top-level class for creating a new 
 * job to run BATON. This class will merely
 * create the various pages and add them to the wizard. Each
 * page performs a specific task required to create/configure
 * a complete job to be run on a server.
 *  
 * <p><b>Note:</b> In order to create a job, there must be at least one
 * data set and one server entry in the work space.</p> 
 */
public class BatonJobWizard extends WizardDialog {
	
	/**
	 * The constructor for the job wizard.
	 * 
	 * The constructor lays out the wizard and creates all the
	 * wizard pages. So far, there is only one. The wizard pages are created as templates
	 * and get populated with the necessary information just before
	 * the pages get displayed to the user.
	 * 
	 * @param title The title to be set for the main wizard frame.
	 * @param parent The main frame that logically owns this wizard.
	 */
	public BatonJobWizard(String title, MainFrame parent) {
		super(parent);
		this.mainFrame = parent;
		setTitle(title);
		setResizable(false);
		// Set up the title image we want to use.
		setTitleBackground("images/peace_wizard_header.png", Color.white);
		// Set up the column image we want to use.
		setSequenceBackground("images/peace_wizard_column.png");
		// First setup the overview page.
		createOverview();
		
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
	 * wizard. Currently this method has no specific task to 
	 * perform but is present as a place holder for future 
	 * enhancements. 
	 */
	@Override
	public void done(boolean success) {
		if (!success) {
			return;
		}
	}
	
	/**
	 * Obtain the instance of the main frame class that owns this wizard.
	 * 
	 * This method is currently used by the submit job wizard page 
	 * to create a job monitor once a job has been successfully
	 * submitted. 
	 * 
	 * @return The main frame that owns this wizard.
	 */
	protected MainFrame getMainFrame() {
		return mainFrame;
	}
	
	/**
	 * The main frame that logically owns this job. This value is
	 * used to create a job monitor if a job is successfully 
	 * submitted to run on a server.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user.
	 */
	private static final String OVERVIEW_MSG = "<html>" +
		"This wizard guides you through the process of creating a job.<br>" +
		"A job runs on a selected server (in parallel is so configured)<br>"+
		"and runs BATON.<br></html>";

	/**
	 * The generated serialization UID (need to keep the compiler happy) 
	 */
	private static final long serialVersionUID = -6700963924724398996L;
}
