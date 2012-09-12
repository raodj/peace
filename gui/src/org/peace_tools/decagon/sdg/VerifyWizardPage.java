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

import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.views.GenericHTMLView;

/**
 * This class serves as the penultimate page in the SyntheticDataSetGenerator. 
 * This page provides the user with the summary of the parameters that will
 * be used to run MetaSim to generate synthetic reads. This is the last page 
 * that the user can backtrack in this wizard.
 */
public class VerifyWizardPage extends GenericWizardPage {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: couple of labels
	 * and a JEditorPane to display HTML summary information.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public VerifyWizardPage(SyntheticDataSetGenerator dcm) {
		this.sdg = dcm;
		assert(this.sdg != null);
		// Setup the title(s) for this page and border
		setTitle("Verification", 
				"Verify parameters for generating synthetic reads");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the job summary information text area.
		summary = GenericHTMLView.createHTMLPane();
		JScrollPane jsp = new JScrollPane(summary);
		// Create the informational labels.
		JLabel caution = new JLabel(CAUTION_MSG, 
				Utilities.getIcon("images/24x24/Warning.png"), 
				JLabel.LEFT);
		
		// Pack the display fields into a suitable panel
		JPanel subPanel = new JPanel(new BorderLayout(5, 5));
		subPanel.add(jsp, BorderLayout.CENTER);
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
		summary.setText(sdg.getSummary());
		summary.setCaretPosition(0);
	}
		
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to obtain summary information
	 * from various pages constituting this wizard.
	 */
	private final SyntheticDataSetGenerator sdg;

	/**
	 * This editor pane displays a summary information as HTML
	 * document. The data in this component is populated by the
	 * {@link #pageChanged(WizardDialog, int, int)} method. 
	 */
	private final JEditorPane summary;

	/**
	 * A simple caution message that is displayed at the bottom of the
	 * wizard page.
	 */
	private static final String CAUTION_MSG = 
		"<html><b>You will not be able to backtrack from the next page.</b></html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
