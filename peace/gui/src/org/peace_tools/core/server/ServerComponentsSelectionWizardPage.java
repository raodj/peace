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

package org.peace_tools.core.server;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.Box;
import javax.swing.Icon;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;

/**
 * This class serves as one of the (namely 3rd page) interactive 
 * pages in the ServerWizard. This page permits the user to select
 * the software components to be installed on the remote server. 
 * The software components available as choices are:
 * 
 * <ul>
 *   <ol>Fresh install of PEACE (installed from the tar ball in 
 *   the GUI). Optionally (later on) the user will be permitted to
 *   reuse an existing install of PEACE or a SVN version of PEACE.</ol>
 *   
 *   <ol>Install the EAST assembler from the tar ball (or use an 
 *   existing copy of EAST, similar to PEACE). This available only
 *   on Linux boxes right now.</ol>
 *   
 *   <ol>An install of DECAGON (Distributed Environment for Comparative
 *    Analysis of Genomic-Assemblers) for comparing various assemblers. 
 *    This is available on Linux machines.</ol>
 *    
 * </ul>
 * 
 */
public class ServerComponentsSelectionWizardPage extends GenericWizardPage 
implements ItemListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include a combo box to
	 * choose between local and remote servers. In addition, there is
	 * a separate panel that is setup for credentials associated
	 * with the remote host. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param server The in-memory server object that contains the
	 * core information about the server.
	 */
	public ServerComponentsSelectionWizardPage(ServerWizard wizard, Server server) {
		this.wizard = wizard;
		this.server = server; 
		// Setup the title(s) for this page and border
		setTitle("Select Components", "Select software components to install");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create a panel with the type combo box for the user to
		// choose from two different options.
		String[] options = {"Install from local GUI package"};
		softCompSrc = new JComboBox(options);
		softCompSrc.setEditable(false);
		// Pack the options along with a pretty icon and label into 
		// a sub-panel.
		JPanel subPanel = new JPanel(new BorderLayout(2, 2));
		subPanel.setBorder(new EmptyBorder(10, 15, 10, 10));
		subPanel.add(new JLabel("Select source of server software:"), BorderLayout.NORTH);
		subPanel.add(softCompSrc, BorderLayout.SOUTH);
		// Add sub-panel to the main panel
		add(subPanel, BorderLayout.NORTH);
		
		// Create a sub panel with choices for the source of software
		// components to install on the server.
		boolean[] defSelect= {true, false, false};
		softCompSelections = defSelect; // Setup default selections
		guiSoftComps       = new JPanel[3];
		// Create component for PEACE selection.
		guiSoftComps[0] = createCompSelPanel(0, "NewCluster", 
				softCompSelections[0], true, false);
		// Create component for EAST selection (if server has 
		guiSoftComps[1] = createCompSelPanel(1, "NewMST", 
				softCompSelections[1], false, true);
		// Create component for DECAGON selection.
		guiSoftComps[2] = createCompSelPanel(2, "Jobs", 
				softCompSelections[2], false, false);
				
		// Pack the various components into a suitable panel.
		JPanel choicePanel = Utilities.createLabeledComponents(null, null, 0, true,
				guiSoftComps[0], Box.createVerticalStrut(10),
				guiSoftComps[1], Box.createVerticalStrut(10),
				guiSoftComps[2]);
		// Add the remote server Information panel to the main page
		add(choicePanel, BorderLayout.CENTER);
	}
	
	private JPanel createInfoPanel(final String infoMsg, final boolean enabled) {
		// Break the information message into individual lines (to determine
		// number of lines).
		final String[] lines = infoMsg.split("\n");
		// Next create the panel to contain labels in a 1 column grid layout.
		JPanel infoPanel = new JPanel(new GridLayout(lines.length, 1));
		// Add labels with suitable font to the infoPanel.
		for(String text: lines) {
			JLabel label = new JLabel(text);
			// Make font a bit smaller.
			Utilities.adjustFont(label, -2, 6, 0);
			// Enable/disable the label as indicated
			label.setEnabled(enabled);
			// Add to the panel
			infoPanel.add(label);
		}
		// Return panel with labels back to the caller.
		return infoPanel;
	}
	
	private JPanel createCompSelPanel(final int choiceIndex, final String iconName,
			final boolean checked, final boolean enabled, final boolean editable) {
		// Create a JCheckBox or a JLabel depending on editable flag.
		JComponent selector = null;
		if (editable) {
			// Make regular checkbox as we want the entry to be editable.
			JCheckBox cb = new JCheckBox(SoftCompLabels[choiceIndex * 2], checked);
			// Setup check box icons to make images look consistent
			cb.setIcon(Utilities.getIcon("images/16x16/Box.png"));
			cb.setSelectedIcon(Utilities.getIcon("images/16x16/CheckedBox.png"));
			// Set up the listener to handle checking/un-checking actions
			cb.addItemListener(this);
			selector = cb;
		} else {
			// Un-editable entry. Create a JLabel with a suitable checkbox-like icon.
			final String cbName = "images/16x16/" + (checked ? "Checked" : "" ) + "Box.png";
			selector = new JLabel(SoftCompLabels[choiceIndex * 2], 
					Utilities.getIcon(cbName), JLabel.LEFT);
		}
		// Make the selector font bold
		Utilities.adjustFont(selector, 0, 0, 1); // Make font bold & prominent.
		// Setup name for the selector to make event processing easier.
		selector.setName("" + choiceIndex);
		// Next create a panel with the help/info text
		JPanel infoPanel = createInfoPanel(SoftCompLabels[choiceIndex * 2 + 1], enabled);
		// Create a fancy icon to go with info message.
		Icon img = (iconName != null) ? Utilities.getIcon("images/24x24/" + iconName + ".png") : null;
		JLabel imgLabel = new JLabel(img);
		// Place brief descriptive messages in infoPanel along with a fancy icon 
		// (to make things look nice)
		JPanel infoIconPanel = new JPanel(new BorderLayout(5, 1));
		infoIconPanel.add(imgLabel,  BorderLayout.WEST);
		infoIconPanel.add(infoPanel, BorderLayout.CENTER);
		
		// Create a JPanel that contains the check box along with the infoIconPanel
		JPanel subPanel = new JPanel(new BorderLayout(16, 0));
		subPanel.add(selector, BorderLayout.NORTH);
		subPanel.add(new JLabel(" "), BorderLayout.WEST);
		subPanel.add(infoIconPanel, BorderLayout.CENTER);
		subPanel.setBorder(new EmptyBorder(0, 15, 0, 0));
		// Ensure entries are appropriately enabled/disabled.
		Utilities.setEnabled(subPanel, enabled);
		return subPanel;
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the data being displayed in
	 * the GUI fields from the data stored in the in-memory
	 * Server object.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Enable/disable check box for EAST depending on whether GCC has been found
		Utilities.setEnabled(guiSoftComps[1], wizard.haveGCC());
		// By default check the EAST check box. For this we need to find
		// the appropriate component.
		Component[] children = guiSoftComps[1].getComponents();
		for(Component child: children) {
			if (child instanceof JCheckBox) {
				JCheckBox eastCB = (JCheckBox) child;
				eastCB.setSelected(wizard.haveGCC());
				break;
			}
		}
	}

	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		// Copy information regard EAST install into the server
		server.setEASTInstalled(softCompSelections[1]);
		return true;
	}

	@Override
	public void itemStateChanged(ItemEvent event) {
		// Get the index of the checkbox that was changed from name of the check box.
		// The name was set in the this.createCompSelPanel() method
		String indexStr = ((JComponent) event.getSource()).getName();
		int compIndex   = Integer.parseInt(indexStr);
		// Now compIndex is 0 for PEACE, 1 for EAST, and 2 for DECAGON
		this.softCompSelections[compIndex] = (event.getStateChange() == ItemEvent.SELECTED);
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final ServerWizard wizard;
	
	/**
	 * Information about the actual server entry being edited. This
	 * reference is set when this class is instantiated and is never
	 * changed during the life time.
	 */
	private final Server server;
	
	/**
	 * The combo-box that permits the user to choose the source for
	 * software components --- between components from local install
	 * package or remote SVN/existing copy. 
	 */
	private JComboBox  softCompSrc;
	
	/**
	 * An array of check-boxes that permits the user to choose the
	 * sub-set of PEACE server software components to be installed
	 * on the server. The values in this array correspond to:
	 * 
	 * <ol>
	 * 	<li>PEACE enabled/disabled (always true to indicate enabled)</li>
	 *  <li>EAST enabled/disabled (set for non-Windows server)</li>
	 *  <li>DECAGON enabled/disabled (always false to indicate disabled)</li>
	 *  </ol>
	 *   
	 * The texts for these check boxes are listed in the static array below.
	 */
	private boolean[] softCompSelections;

	/**
	 * The list of GUI-components that contain a label and check box for the
	 * software components being installed. These GUI components are
	 * used to merely enable/disable labels to provide the user with
	 * better visual feedback.
	 */
	private JPanel[] guiSoftComps;
	
	/**
	 * An array of static labels that are used to provide the users
	 * with information about the software components that are to be 
	 * installed on the server. 
	 */
	private static final String SoftCompLabels[] = {
		"PEACE Clustering Engine",
		"This is a required software component that must be\n" +
		"be installed as part of the server software components.\n",
		"EAST Gene Assembler",
		"Assembles ESTs into Genes using the clustering results.\n"  +
		"This component is available only on Linux servers.",
		"DECAGON assembler evaluation environment",
		"ADistributed Environment for Comparative Analysis\n"  +
		"of Genomic-Assemblers (DECAGON) used to compare/contrast PEACE\n" +
		"with other software using various statistical metrics."
	};
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
