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

package org.peace_tools.decagon.configMgr;

import java.awt.BorderLayout;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.PEACEProperties;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

/**
 * This class serves as the last page in the DecagonConfigurationManager. 
 * This page commits the configuration information for DECAGON to the
 * global properties and saves the properties. All the operations of
 * this class are performed once the page has been displayed.
 */
public class CommitConfigPage extends GenericWizardPage 
implements Runnable, PropertyChangeListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: couple of labels
	 * and a JEditorPane to display HTML summary information.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public CommitConfigPage(DecagonConfigurationManager dcm) {
		this.dcm = dcm;
		assert(this.dcm != null);
		// Setup the title(s) for this page and border
		setTitle("Commit Configuration", 
				"Update and save DECAGON configuration");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// This page has just one label whose information is changed
		// once the operations have completed.
		infoLabel = new JLabel(START_INFO_MSG, 
				Utilities.getIcon("images/32x32/Information.png"), 
				JLabel.LEFT);
		// Pack the display fields into a suitable panel
		JPanel subPanel = new JPanel(new BorderLayout(5, 5));
		subPanel.add(infoLabel, BorderLayout.CENTER);
		// Finally add the sub panel to this panel.
		add(subPanel, BorderLayout.CENTER);
		// Add ourselves as property change listener to receive
		// notification once background task completes.
		addPropertyChangeListener(this);
	}
	

	@Override
	public void run() {
		dcm.commitProperties();
		final Exception exp = PEACEProperties.get().save(this, false);
		// Done. Notify others on Swing's event processing thread
		firePropertyChange("done", null, exp);
	}

	@Override
	public void propertyChange(PropertyChangeEvent pce) {
		final String prop = pce.getPropertyName();
		if ("done".equals(prop)) {
			// Enable the finish button.
			dcm.setButtonStatus(0, 1, 0);
			// Report any exception that we may have encountered.
			final Exception exp = (Exception) pce.getNewValue();
			if (exp != null) {
				// An error occurred. Report it to the user in the panel.
				final String htmlMessage = String.format(ERROR_INFO_MSG,
						Utilities.wrapStringToHTML(exp.getMessage(), 65));
				final JPanel info = Utilities.collapsedMessage(htmlMessage, Utilities.toString(exp));
				this.add(info, BorderLayout.CENTER);
			} else {
				// Successfully completed saving. Update the informational
				// label in this panel.
				infoLabel.setText(SUCCESS_INFO_MSG);
			}
		}
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the summary information displayed
	 * in this wizard page.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Perform the commit operations in the background to ensure that
		// the GUI continues to remain responsive. First disable all the
		// buttons while we commit the changes.
		dcm.setButtonStatus(0, 0, 0);
		// Start the operations in the background.
		SwingUtilities.invokeLater(this);
	}
		
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to get the various pages
	 * constituting this wizard to commit their configuration.
	 */
	private final DecagonConfigurationManager dcm;

	/**
	 * The label that is used to display informational messages.
	 */
	private final JLabel infoLabel;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some contextual information
	 * to the user.
	 */
	private static final String START_INFO_MSG = 
		"<html>The DECAGON configuration information is being updated<br>" +
		"and saved. Please wait..." +
		"</html>";

	/**
	 * A simple message that is displayed in the panel if the property
	 * update and saving is successfully completed.
	 */
	private static final String SUCCESS_INFO_MSG = 
		"<html>The DECAGON configuration has been successfully<br>" +
		"updated and saved. DECAGON configuration is complete.<br>" +
		"You may now use the various features that you have<br>" +
		"enabled via the DECAGON main menu options." +
		"</html>";

	/**
	 * A simple message that is displayed in the panel if the property
	 * update and saving experienced errors.
	 */
	private static final String ERROR_INFO_MSG = 
		"<html><b>The DECAGON configuration could not be saved!</b><br>" + 
		"<dl><dt>The following error occured:</dt>" +
		"<dd>%s</dd></dl><br>" +
		"<i>You will need to appropriately remedy the cause of<br>" +
		"of this error and reconfigure DECAGON again.</i>" +
		"</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
