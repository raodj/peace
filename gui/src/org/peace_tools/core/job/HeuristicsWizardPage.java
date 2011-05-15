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
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.FWAnalyzer;
import org.peace_tools.workspace.Heuristic;
import org.peace_tools.workspace.Param;

/**
 * This class serves as an interactive page in a JobWizard.
 * This page permits the user to provide the information about the
 * frame/word analyzer to be used to computing distance/similarity
 * metric between two given ESTs. The distance/similarity metric
 * is used to decide if two ESTs must be clustered together. 
 * 
 */
public class HeuristicsWizardPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: combo box to
	 * select the frame/word analyzer to be used, spinners for window
	 * and word size, the type of cache and size of cache to be used.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * @param awp    The analyzer wizard page from where the frame size
	 * for t/v heuristic is obtained.
	 */
	public HeuristicsWizardPage(JobWizard wizard, AnalyzerWizardPage awp) {
		this.wizard = wizard;
		this.awp    = awp;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Setup Heuristics", 
				"Configure the heuristics to use");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the informational label.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/16x16/Information.png"), 
				JLabel.LEFT);
		info.setAlignmentX(0);
		// Pack the input fields into a box
		heuristicPanel = Utilities.createLabeledComponents(null, null, 0, true,
				info, Box.createVerticalStrut(5),
				createUvPanels(),
				Box.createVerticalStrut(5),
				createTvPanels());
		// Set border to make layout look good.
		heuristicPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(heuristicPanel, BorderLayout.CENTER);
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
	 * information regarding the u/v heuristic.
	 */
	private JPanel createUvPanels() {
		// First create enable/disable check box with nice border.
		enableUV = new JCheckBox("Enable u/v sample heuristic for this job");
		JPanel bag = adjustCheckBox(enableUV, "uv");
		// Create and add u/v heuristic parameters.
		Box horizBox = Box.createHorizontalBox();
		horizBox.setAlignmentX(0);
		// Create the spinners for u, v, and fame size and add them
		// to the horiz box with a label
		final int DefaultValues[] = {8, 4, 16};
		final String Titles[] = {"V (bp):", "U (bp):", "Word Shift (bp):"};
		for(int i = 0; (i < DefaultValues.length); i++) {
			SpinnerNumberModel model = 
				new SpinnerNumberModel(DefaultValues[i], 2, 32, 1);
			// Create and configure the spinner.
			uvParams[i] = new JSpinner(model); 
			Utilities.adjustDimension(uvParams[i], 50, 4); // Adjust size to look right
			if (i > 0) {
				horizBox.add(Box.createHorizontalStrut(20));
			}
			horizBox.add(Utilities.createLabeledComponents(Titles[i], null,
					0, false, uvParams[i]));
		}
		// Now put all the information into a nice titled panel.
		JPanel bigBox = Utilities.createLabeledComponents(null, null, 0, false,
				bag, Box.createVerticalStrut(5), 
				new JLabel(UV_INFO_MSG), 
				Box.createVerticalStrut(5),
				horizBox);
		// Set border to make things look good.
		bigBox.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createTitledBorder("u/v Heuristic"),
				BorderFactory.createEmptyBorder(5, 5, 5, 5)));
		return bigBox;
	}
	
	/**
	 * This is a helper method that is used to create the necessary
	 * information regarding the t/v heuristic.
	 * 
	 * @return A panel to which the text controls have been added. 
	 */
	private JPanel createTvPanels() {
		// First create enable/disable checkbox with nice border.
		enableTV = new JCheckBox("Enable t/v sample heuristic (needs u/v)");
		JPanel bag = adjustCheckBox(enableTV, "tv");
		// Create and add t/v heuristic parameters.
		Box horizBox = Box.createHorizontalBox();
		horizBox.setAlignmentX(0);
		final int DefaultValues[] = {100, 65};
		final String Titles[] = {"v (frame size in bp):", "T (#v-word matches):"};
		for(int i = 0; (i < DefaultValues.length); i++) {
			SpinnerNumberModel model = 
				new SpinnerNumberModel(DefaultValues[i], 2, 1000, 1);
			// Create and configure the spinner.
			tvParams[i] = new JSpinner(model); 
			Utilities.adjustDimension(tvParams[i], 100, 4); // Adjust size to look right
			if (i > 0) {
				horizBox.add(Box.createHorizontalStrut(30));
			}
			horizBox.add(Utilities.createLabeledComponents(Titles[i], null,
					0, false, tvParams[i]));
		}
		// The first parameter is locked based on analyze frame size
		tvParams[0].setEnabled(false);

		// Now put all the information into a nice titled panel.
		JPanel bigBox = Utilities.createLabeledComponents(null, null, 0, false,
				bag, Box.createVerticalStrut(5), horizBox);
		// Set up border 
		bigBox.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createTitledBorder("t/v Heuristic"),
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
		// Check and set controls based on uv heuristic.
		if (!enableUV.isSelected()) {
			// Disable everything on this page but for the checkbox.
			Utilities.setEnabled(this, false);
			enableUV.setEnabled(true);
		} else {
			// Enable every control.
			Utilities.setEnabled(this, true);
			// Disable tv's first input.
			tvParams[0].setEnabled(false);
		}
		// Check and set control based on t/v heuristic
		tvParams[1].setEnabled(enableTV.isSelected());
	}

	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the frame size for the t/v
	 * heuristic from the data provided by the analyzer wizard page. 
	 * 
	 * In addition, the heuristic sub-panel is hidden (or displayed)
	 * if the heuristic is TwoPassD2 or CLU indicating that these
	 * analyzers configure their own heuristics.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Remove component at the center in preparation for new one below.
		BorderLayout bl = (BorderLayout) getLayout();
		Component currComp = bl.getLayoutComponent(this, BorderLayout.CENTER);
		this.remove(currComp);
		// Display suitable message or the heuristic configuration
		// information depending on the type of heuristic chosen earlier.
		if (awp.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2)) {
			JLabel msg = new JLabel(TWO_PASS_INFO, Utilities.getIcon("images/32x32/Information.png"),
					JLabel.LEFT);
			msg.setIconTextGap(10);
			msg.setVerticalAlignment(JLabel.CENTER);
			add(msg, BorderLayout.CENTER);
		} else if (awp.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.CLU)) {
			JLabel msg = new JLabel(CLU_INFO, Utilities.getIcon("images/32x32/Information.png"),
					JLabel.LEFT);
			msg.setIconTextGap(10);
			msg.setVerticalAlignment(JLabel.CENTER);
			add(msg, BorderLayout.CENTER);
		} else {
			// Update the fixed window size
			tvParams[0].setValue(awp.getFrameSize());
			// Disable check box to enable/disable heuristics in twopass d2
			// non-adaptive as as it uses uv and tv heuristics no-matter what.
			if (awp.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2_DONT_ADAPT)) {
				enableUV.setSelected(true);
				enableUV.setEnabled(false);
				enableTV.setSelected(true);
				enableTV.setEnabled(false);
			} else {
				// Let user decide if heuristics should be used.
				enableUV.setEnabled(true);
				enableTV.setEnabled(true);
			}
			// Display the panel for configuring heuristics
			add(heuristicPanel, BorderLayout.CENTER);
		}
		// Ensure any changes are picked up and displayed
		validate();
	}
	
	/**
	 * This method is a convenience method that can be used to obtain
	 * a HeuristicList object containing information from this page. This
	 * method is used by the JobWizard to compose a complete set of
	 * entries once the wizard has successfully completed.
	 * 
	 * @return A HeuristicList object containing the necessary information
	 * about the heuristics configured by the user. If the user has chosen
	 * two-pass-d2 then this method returns null indicating default heuristics
	 * are to be used.
	 */
	protected ArrayList<Heuristic> getHeuristics() {
		if (awp.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2)) {
			// For two-pass-d2 we use default heuristics by returning null here.
			// Note that there is a big difference between null return (indicating
			// default heuristics) versus empty heuristic list (indicating no
			// heuristics at all)
			return null;
		}
		ArrayList<Heuristic> heurList = new ArrayList<Heuristic>(2);
		final boolean useHeuristics = 
			awp.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2_DONT_ADAPT) ||
			awp.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.D2) ||
			awp.getAnalyzerType().equals(FWAnalyzer.FWAnalyzerType.D2ZIM);
		// Populate list with necessary entries.
		if (enableUV.isSelected() && (useHeuristics)) {
			Heuristic uv = new Heuristic("uv");
			uv.addParameter(new Param("uv_v", uvParams[0].getValue().toString()));
			uv.addParameter(new Param("uv_u", uvParams[1].getValue().toString()));
			uv.addParameter(new Param("uv_wordShift", uvParams[2].getValue().toString()));
			// Add uv heuristic entry to the list.
			heurList.add(uv);
		}
		if (enableUV.isSelected() && enableTV.isSelected() && (useHeuristics)) {
			Heuristic tv = new Heuristic("tv");
			tv.addParameter(new Param("tv_win", tvParams[0].getValue().toString()));
			tv.addParameter(new Param("tv_t", tvParams[1].getValue().toString()));
			// Add uv heuristic entry to the list.
			heurList.add(tv);
		}
		// Return list of heuristics to use.
		return heurList;
	}
	
	/**
	 * This method is a convenience method that can be used to obtain
	 * a String summarizing the heuristic setup for this job. 
	 * This method is used by the JobWizard to compose a complete set of
	 * summary information just before the job is submitted.
	 * 
	 * @return A String containing a summary information of the various
	 * heuristics configured for this run by the user.
	 */
	protected String getHeuristicsSummary() {
		StringBuilder sb = new StringBuilder();
		// Obtain heuristic list and get summary information from there.
		ArrayList<Heuristic> heuristicList = getHeuristics();
		// Dump summary information.
		for(Heuristic heuristic: heuristicList) {
			sb.append(heuristic.getSummary());
		}
		// Return summary string.
		return sb.toString();
	}
	
	/**
	 * This panel actually contains all the controls for setting up the
	 * necessary heuristics and their parameters required by certain 
	 * analyzers. Some of the analyzers currently do not permit setting or
	 * modifications of heuristics. In this case this panel is not displayed
	 * to the user. Instead a simple message indicating that there are no
	 * heuristics to be configured is displayed to the user. 
	 */
	private final JPanel heuristicPanel;
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final JobWizard wizard;

	/**
	 * A reference to the analyzer wizard page from where the the
	 * frame size set for the analyzer is obtained and used to
	 * seed the frame size for the t/v heuristic.
	 */
	private final AnalyzerWizardPage awp;

	/**
	 * Check box to indicate if the u/v sample heuristic is to be
	 * used for this job.
	 */
	private JCheckBox enableUV;

	/**
	 * Check box to indicate if the t/v heuristic is to be
	 * used for this job. Note that t/v heuristic needs u/v
	 * heuristic to be enabled.
	 */
	private JCheckBox enableTV;

	/**
	 * The array of three configuration values for the u/v 
	 * heuristic. The first entry is v (the length of common words
	 * to search for), the second entry is u (number of v-word matches)
	 * and word shift (number of words to shift by)
	 */
	private JSpinner uvParams[] = new JSpinner[3];

	/**
	 * The array of two configuration parameter values for the t/v 
	 * heuristic. The first entry is frame size (it is fixed based on 
	 * frame size of analyzer), the second entry is t (number of 
	 * v-word matches).
	 */
	private JSpinner tvParams[] = new JSpinner[2];

	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"It is best to use both heuristics for maximum performance.";

	/**
	 * A generic informational message that is displayed to the user
	 * to provide information about the u/v heuristic.
	 */
	private static final String UV_INFO_MSG = 
		"<html>This heursitc requires: <i>v</i> (#common words), <i>u</i> (# v-word<br>" +
		"matches) and <i>shift</i> (# of bases to shift/skip)</html>";

	/**
	 * A simple text message that is displayed to the user when the user 
	 * selects the two-passed D2 configuration. The text is present to ensure
	 * that the user is clear about the operations being performed. The text is
	 * used in the {@link #pageChanged(WizardDialog, int, int)} method.
	 */
	private static final String TWO_PASS_INFO = "<html>" +
		"Two Pass D2 analyzer is adaptive and uses several different ranges<br>" +
		"of window (aka frame) sizes that are most suitable to analyze<br>" +
		"a given pair of fragments. In addition, the analyzer suitably<br>" +
		"sets the values and thresholds for the <i>u/v</i> and the <i>t/v</i><br>" +
		"heuristics. Consequently, there is no additional heuristic setup or<br>" +
		"configuration needed for this analyzer. PEACE will automatically<br>" +
		"default to the best, adaptive choices." +
		"</html>";
	
	/**
	 * A simple text message that is displayed to the user when the user 
	 * selects the CLU analyzer. The text is present to ensure
	 * that the user is clear about the operations being performed. The text is
	 * used in the {@link #pageChanged(WizardDialog, int, int)} method.
	 */
	private static final String CLU_INFO = "<html>" +
		"This analyzer currently does not use any heuristics.<br>" +
		"Therefore there are no heuristics to configure on this page." +
		"</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
