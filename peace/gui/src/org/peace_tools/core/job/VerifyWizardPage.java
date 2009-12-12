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

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

/**
 * This class serves as the penultimate page in a JobWizard. This page
 * provides the user with a brief summary of the job that is going
 * to be scheduled. This is the last page that the user can back
 * track or cancel the job.
 * 
 */
public class VerifyWizardPage extends GenericWizardPage {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a copule of labels
	 * and a text area for job summary information.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public VerifyWizardPage(JobWizard wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Job Summary", 
				"Verify information and submit job");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the job summary information text area.
		summary = new JTextArea(3, 10);
		JComponent descBox = 
			Utilities.createLabeledComponents("Job information summary:",
					"(This information cannot be edited)", 0, 
					new JScrollPane(summary));
		// Create the informational labels.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/16x16/Information.png"), 
				JLabel.LEFT);
		JLabel caution = new JLabel(CAUTION_MSG, 
				Utilities.getIcon("images/16x16/Warning.png"), 
				JLabel.LEFT);
		
		// Pack the display fields into a suitable panel
		JPanel subPanel = new JPanel(new BorderLayout(5, 5));
		subPanel.add(info, BorderLayout.NORTH);
		subPanel.add(descBox, BorderLayout.CENTER);
		subPanel.add(caution, BorderLayout.SOUTH);
		// Finally add the sub panel to this panel.
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the summary information displayed
	 * in this wizard page.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Update the summary information in the summary box.
		summary.setText(wizard.getSummary());
		summary.setCaretPosition(0);
	}
		
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final JobWizard wizard;

	/**
	 * Field to read and edit a brief description about the job.
	 * This information can be anything the user desires and is
	 * meaningful only to the user.
	 */
	private JTextArea summary;

	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some contextual information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"<html>You have provided the following information regarding<br>" +
		"the job to be scheduled. Please verify the information. Next,k<br>" +
		"'the wizard will perform the initial job setup.</html>";

	/**
	 * A simple caution message that is displayed at the bottom of the
	 * wizard page.
	 */
	private static final String CAUTION_MSG = 
		"<html><b>On clicking the 'Next' button the job will be submitted!</b><br>" +
		"You will not be able to backtrack from the next page.</html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
