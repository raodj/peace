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

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;
import javax.swing.text.JTextComponent;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataFileStats;
import org.peace_tools.workspace.DataSet;

/**
 * This class serves as the last page in a DataSetWizard. This page
 * gives a summary information to the user about the EST file and 
 * the data set to be created.  This page provides the user the
 * last chance to cancel out of the data set addition process. If the
 * user proceeds, then the data set wizard adds the new data set
 * entry to the work space.
 */
public class SummaryWizardPage extends GenericWizardPage { 
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include a set of labels
	 * that display information to the user. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * @param dataSet The in-memory data set object that contains the
	 * core information.
	 */
	public SummaryWizardPage(DataSetWizard wizard, DataSet dataSet) {
		this.wizard  = wizard;
		this.dataSet = dataSet;
		// Setup the title(s) for this page and border
		setTitle("Summary", "Check & verify new data set information");
		setBorder(new EmptyBorder(30, 10, 5, 10));
		// Create and set properties for the text fields.
		infoFields[0] = new JTextField();
		infoFields[1] = new JTextArea(5, 5);
		for(int i = 0; (i < infoFields.length); i++) {
			infoFields[i].setEditable(false);
			infoFields[i].setBackground(Color.white);
			infoFields[i].setBorder(BorderFactory.createEtchedBorder());
		}
		// Pack all the labels into a panel with vertical box layout.
		JPanel container = 
		Utilities.createLabeledComponents(STANDARD_MSG, null, 4, true,
			new JLabel(" "), // Blank line
			new JLabel("EST file name:"), infoFields[0],
			new JLabel(" "), // Blank line
			new JLabel("EST data summary:"), infoFields[1],
			new JLabel(" "), // Blank line
			new JLabel(INFO_MSG, 
					Utilities.getIcon("images/16x16/Information.png"),
					JLabel.LEFT)
		);
		// Add container to this page
		add(container, BorderLayout.CENTER);
	}
			
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the data being displayed in
	 * the GUI fields from the data stored in the in-memory
	 * Server object.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// First update the necessary information.
		infoFields[0].setText(dataSet.getPath());
		// Populate the text area with data on the EST file.
		DataFileStats stats = wizard.getESTList().computeStatistics();
		// Convert the stats to a string.
		String statistics = String.format(STATS_STR, stats.getCount(), 
				stats.getAvgLength(), stats.getLengthSD(), 
				stats.getMinLength(), stats.getMaxLength());
		// Set it up for display.
		infoFields[1].setText(statistics);
		// Setup stats in the data set
		dataSet.setStats(stats);
	}

	/**
	 * A simple standard message to be displayed at the top of this page.
	 */
	private final String STANDARD_MSG =  "<html>A new data set will be created " +
		"with using the<br>following EST data file:</html>";
	
	/**
	 * A simple standard message to be displayed at the bottom of this
	 * page.
	 */
	private final String INFO_MSG =  "Click the Finish button to " +
		"add the new data set.";
	
	/**
	 * A simple standard message to be displayed at the bottom of this
	 * page.
	 */
	private final String STATS_STR =  
		" Number of ESTs     :   %1$d\n"    +
		" Average  EST length:   %2$f bp (SD: %3$f)\n" +
		" Shortest EST length:   %4$d bp\n" +
		" Longest  EST length:   %5$d bp";

	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final DataSetWizard wizard;
	
	/**
	 * Information about the actual data set entry being edited. This
	 * reference is set when this class is instantiated and is never
	 * changed during the life time.
	 */
	private final DataSet dataSet;
			
	/**
	 * Field to be displayed to the user in this wizard pane. The
	 * first field contains the name of the EST data file and the
	 * second field is a multi-line text area that contains some
	 * details about the EST file.
	 */
	private JTextComponent infoFields[] = new JTextComponent[2];
		
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
