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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Filter;
import org.peace_tools.workspace.Param;

/**
 * This class serves as an interactive page in a JobWizard.
 * This page permits the user to provide the information about the
 * filters to be used in the job. Filters are used to weed out ESTs
 * with undesirable characteristics (biologically speaking) that 
 * negatively impact with the quality of clusters generated
 * by PEACE.  Currently the following two filters are available in
 * the clustering engine (and are exposed by the GUI):
 * 
 * <ul>
 *
 * <li>Length Filter: This filter is used to weed out short ESTs 
 * below a given length. Currently, the minimum length is set to 50.</li>
 * 
 * <li>LC Filter: Low Complexity Filter is used to remove ESTs that
 * have low complexity regions in them. The low complexity regions 
 * often result in unnecessary linkages between otherwise unrelated
 * clusters.</li>
 * 
 * </ul>
 * 
 */
public class FiltersWizardPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: check boxes to 
	 * enable/disable either of the two filters. In addition there are
	 * suitable controls that can be used to enter data for the filter
	 * if the filter is enabled.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public FiltersWizardPage(JobWizard wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Setup Filters", 
				"Configure the filters to use");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the informational label.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/16x16/Information.png"), 
				JLabel.LEFT);
		info.setAlignmentX(0);
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
				info, Box.createVerticalStrut(5),
				createLenFilterPanels(),
				Box.createVerticalStrut(5),
				createLCFilterPanels());
		// Set border to make layout look good.
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * Helper method to correctly format, enhance, and layout a given
	 * check box. 
	 * 
	 * @param cb The check box to be formatted.
	 * @param cmd The action command to be set.
	 * @return A JPanel containing the check box with additional
	 * decorations.
	 */
	private JPanel adjustCheckBox(JCheckBox cb, String cmd) {
		cb.setBackground(Color.white);
		cb.setSelected(true);
		cb.setActionCommand(cmd);
		cb.setAlignmentX(0);
		cb.addActionListener(this);
		JPanel bag = new JPanel(new BorderLayout(0, 0));
		bag.setBackground(Color.white);
		bag.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createEtchedBorder(), 
				BorderFactory.createEmptyBorder(1, 2, 1, 2)));
		bag.add(cb, BorderLayout.NORTH);
		bag.setAlignmentX(0);
		return bag;
	}
	/**
	 * This is a helper method that is used to create the necessary
	 * information regarding the length filter. This method is 
	 * invoked from the constructor to layout this page in the wizard.
	 * This is a helper method that was introduced to streamline the
	 * code in the constructor.
	 */
	private JPanel createLenFilterPanels() {
		// First create enable/disable checkbox with nice border.
		enableLenFilter = new JCheckBox("Enable length filter for this job");
		JPanel bag = adjustCheckBox(enableLenFilter, "lenFilter");
		// Create and add the parameter.
		Box horizBox = Box.createHorizontalBox();
		horizBox.setAlignmentX(0);
		// Create the input spinner.
		SpinnerNumberModel model = new SpinnerNumberModel(50, 50, 10000, 5);
		minLen = new JSpinner(model); 
		Utilities.adjustDimension(minLen, 50, 4); // Adjust size to look right
		horizBox.add(Utilities.createLabeledComponents("Set minimum length",
				"(Fragments shorter than this value are filtered out)",
					0, false, minLen));

		// Now put all the information into a nice titled panel.
		JPanel bigBox = Utilities.createLabeledComponents(null, null, 0, false,
				bag, Box.createVerticalStrut(5), 
				horizBox);
		// Set border to make things look good.
		bigBox.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createTitledBorder("Length Filter"),
				BorderFactory.createEmptyBorder(5, 5, 5, 5)));
		return bigBox;
	}
	
	/**
	 * This is a helper method that is used to create the necessary
	 * information regarding the Low Complexity (LC) filter.
	 * 
	 * @return A panel to which the text controls have been added. 
	 */
	private JPanel createLCFilterPanels() {
		// First create enable/disable checkbox with nice border.
		enableLCFilter = new JCheckBox("Enable Low Complexity (LC) Filter");
		JPanel bag = adjustCheckBox(enableLCFilter, "lcFilter");
		// Create the filter list parameter.
		filterList = new JTextField(10);
		filterList.setText("A,C");
		JComponent listBox  =
            Utilities.createLabeledComponents("Enter filter strings:",
                            "(comma separated list of neucleotide values)",
                            0, false, filterList);

		// Now put all the information into a nice titled panel.
		JPanel bigBox = Utilities.createLabeledComponents(null, null, 0, false,
				bag, Box.createVerticalStrut(5), listBox);
		// Set up border 
		bigBox.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createTitledBorder("LC Filter"),
				BorderFactory.createEmptyBorder(5, 5, 5, 5)));
		// Return panel back to caller
		return bigBox;
	}

	/**
	 * Method to handle clicking of check boxes. This method essentially
	 * enables and disables the various inputs depending on the check box.
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		// Check and enable/disable the minimum length input.
		minLen.setEnabled(enableLenFilter.isSelected());
		// Check and enable/disable LC filter string inputs
		filterList.setEnabled(enableLCFilter.isSelected());
	}
	
	/**
	 * This method is a convenience method that can be used to obtain
	 * a FilterList object containing information from this page. This
	 * method is used by the JobWizard to compose a complete set of
	 * entries once the wizard has successfully completed.
	 * 
	 * @return A FilterList object containing the necessary information
	 * about the filters configured by the user.
	 */
	protected ArrayList<Filter> getFilters() {
		ArrayList<Filter> filtList = new ArrayList<Filter>(2);
		// Populate list with necessary entries.
		if (enableLenFilter.isSelected()) {
			Filter filter = new Filter(Filter.FilterType.LengthFilter);
			filter.addParameter(new Param("minESTLen", minLen.getValue().toString()));
			// Add length filter entry to the list.
			filtList.add(filter);
		}
		if (enableLCFilter.isSelected()) {
			Filter filter = new Filter(Filter.FilterType.LCFilter);
			filter.addParameter(new Param("lcPatterns", filterList.getText()));
			// Add LC filter entry to the list.
			filtList.add(filter);
		}
		// Return list of heuristics to use.
		return filtList;
	}
	
	/**
	 * This method is a convenience method that can be used to obtain
	 * a String summarizing the filters setup for this job. 
	 * This method is used by the JobWizard to compose a complete set of
	 * summary information just before the job is submitted.
	 * 
	 * @return A String containing a summary information of the various
	 * filters configured for this run by the user.
	 */
	protected String getSummary() {
		StringBuilder sb = new StringBuilder();
		// Obtain heuristic list and get summary information from there.
		ArrayList<Filter> filtList = getFilters();
		// Dump summary information.
		for(Filter filter: filtList) {
			sb.append(filter.getSummary());
		}
		// Return summary string.
		return sb.toString();
	}
	
	/**
	 * Method to veto a page change if filter string is invalid.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before the user changes from this page. This method checks to
	 * ensure that the filter string specified by the user is valid.
	 * If the string is valid then this method returns true. 
	 * 
	 * @return This method returns true if the filter text is valid.
	 * Otherwise it returns false.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if ((nextPage < currPage) || (!enableLCFilter.isSelected())) {
			return true;
		}
		// Check if the filter list is valid. The following logic ensures
		// that a comma occurs only when expected.
		boolean commaOK  = false;
		boolean filterOK = true;
		String filterStr = filterList.getText();
		for(int i = 0; (i < filterStr.length()); i++) {
			if (filterStr.charAt(i) == ',')  {
				if (commaOK) {
					commaOK = !commaOK;
				} else {
					filterOK = false;
					break;
				}
			} else {
				if ("ATCG".indexOf(filterStr.charAt(i)) == -1) {
					// Invalid character.
					filterOK = false;
					break;
				} else {
					// After a character we can have a comma
					commaOK = true;
				}
			}
		}
		// At end of list, commaOK should be true. Otherwise we have trouble.
		if ((!filterOK) || (!commaOK)) {
			JOptionPane.showConfirmDialog(wizard, FILTER_LIST_MSG, 
					"Invalid Filter String", JOptionPane.ERROR_MESSAGE);
			filterOK = false;
		}
		// Return filter status
		return filterOK;
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final JobWizard wizard;

	/**
	 * Check box to indicate if the length filter must be enabled 
	 * for this job.
	 */
	private JCheckBox enableLenFilter;

	/**
	 * The minimum length of the EST associated with the length filter.
	 * ESTs below this length are filtered out by the length filter 
	 * (if it is enabled).
	 */
	private JSpinner minLen;
	
	/**
	 * Check box to indicate if the Low Complexity (LC) filter must
	 * be enabled for this job.
	 */
	private JCheckBox enableLCFilter;

	/**
	 * A list of comma separated strings (over ATCG) that are used to 
	 * generate dummy low complexity sequences for filtering out ESTs
	 * that have low complexity regions.
	 */
	private JTextField filterList;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"It is best to use both filters for high quality clustering.";

	/**
	 * A generic informational message that is displayed to the user
	 * if the string list in filterList is invalid.
	 */
	private static final String FILTER_LIST_MSG = 
		"<html>The filter list you provided is not valid.<br><br>" +
		"This filter requires a comma separated list of substrings<br>" +
		"that are repeated (to required length) to detect low complexity<br>" +
		"regions in fragments. Example: <i>A,C</i> or <i>AG,A,C,TC</i>.</html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 6691921082180723745L;
}
