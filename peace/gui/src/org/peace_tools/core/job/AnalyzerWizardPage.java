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
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.FWAnalyzer;

/**
 * This class serves as an interactive page in a JobWizard.
 * This page permits the user to provide the information about the
 * frame/word analyzer to be used to computing distance/similarity
 * metric between two given ESTs. The distance/similarity metric
 * is used to decide if two ESTs must be clustered together. 
 * 
 * When the "Next >" button is clicked this wizard populates the
 * user supplied information in FWAnalyzer data structure.
 */
public class AnalyzerWizardPage extends GenericWizardPage implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: combo box to
	 * select the frame/word analyzer to be used, spinners for window
	 * and word size, the type of cache and size of cache to be used.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public AnalyzerWizardPage(JobWizard wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Analyzer Setup", 
				"Configure the EST analyzer to use");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the combo-box with list of analyzers.
		analyzerList = new JComboBox(ANALYZERS);
		analyzerList.setBackground(Color.white);
 		analyzerList.addActionListener(this);	
		// Create panel with the combo box for the user to choose
		JComponent analyzerBox = 
			Utilities.createLabeledComponents("Select analyzer for EST comparison:",
					"(Analyzers provide distance or similarity metrics for clustering)", 0, 
					false, analyzerList);
		// Create the informational label.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/16x16/Information.png"), 
				JLabel.LEFT);
		
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
				info, Box.createVerticalStrut(10), 
				analyzerBox, Box.createVerticalStrut(10),
				createFrameWordPanels(), 
				Box.createVerticalStrut(10),
				createCachePanels());
		// Set border to layout things nicely
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * This is a helper method that is used to create the frame and
	 * word entry boxes. This method adds the two controls with 
	 * suitable labels to the supplied sub-panel.
	 * 
	 * @return A panel to which the text controls have been added.
	 */
	private JPanel createFrameWordPanels() {
		// Create the frame size spinner.
		frameSize = new JSpinner(new SpinnerNumberModel(100, 10, 1000, 5));
		frameSize.setVisible(false);
		Utilities.adjustDimension(frameSize, 0, 6); // Adjust size to look right
		// Create the dummy "Auto" label that is used instead of the frameSize 
		// spinner when the analyzer is two pass d2 based analyzers
		frameSizeAuto = new JLabel("Auto/Adaptive", Utilities.getIcon("images/16x16/Job.png"), JLabel.LEFT);
		frameSizeAuto.setBackground(frameSize.getBackground());
		frameSizeAuto.setOpaque(true);
		frameSizeAuto.setBorder(BorderFactory.createCompoundBorder(frameSize.getBorder(), 
				BorderFactory.createEmptyBorder(2, 5, 2, 5)));
		frameSizeAuto.setToolTipText("Window size is automatically detected in an adaptive manner");
		// Wrap the above two mutually-exclusive controls into a panel so that
		// one of them can be shown easily.
		JPanel wrapper = new JPanel(new BorderLayout(0, 0));
		wrapper.add(frameSize, BorderLayout.NORTH);
		wrapper.add(frameSizeAuto, BorderLayout.SOUTH);
		
		// Create the word size selection combo box 
		wordSize = new JComboBox(new String[]{"6 bp"});
		wordSize.setBackground(Color.white);
		// Put the above two components into a horizontal box.
		JPanel horizBox = new JPanel(new GridLayout(1, 2, 20, 0));
		horizBox.add(Utilities.createLabeledComponents("Window size (bp):", 
				null, 0, false, wrapper));
		horizBox.add(Utilities.createLabeledComponents("Word size (bp):", 
				null, 0, false, wordSize));
		horizBox.setBorder(new EmptyBorder(3, 0, 0, 0));
		// Put the horizontal box into a top-level labeled box.
		return Utilities.createLabeledComponents("Set window & frame size for analyzer:",
				"(influences speed, memory, and sensitivity)", 0, false, horizBox);
	}
	
	/**
	 * This is a helper method that is used to create the cache type
	 * and cache size entry boxes. This method adds the two controls
	 * with suitable labels to the supplied sub-panel. This method
	 * was introduced to keep the code clutter to a minimum in the
	 * constructor.
	 */
	private JPanel createCachePanels() {
		// Create the frame size spinner.
		cacheList = new JComboBox(new String[]{"Heap (best bet)", 
				"Multi List (alternative)"
		});
		cacheList.setBackground(Color.white);
		// Create the cache size selection combo box 
		cacheSize = new JSpinner(new SpinnerNumberModel(128, 32, 1024, 16)); 
		Utilities.adjustDimension(cacheSize, 0, 6); // Adjust size to look right
		// Put the above two components into a horizontal box.
		JPanel horizBox = new JPanel(new GridLayout(1, 2, 20, 0));
		horizBox.add(Utilities.createLabeledComponents("Cache Type:", null, 
				0, false, cacheList));
		horizBox.add(Utilities.createLabeledComponents("Cache size:", null, 
				0, false, cacheSize));
		horizBox.setBorder(new EmptyBorder(3, 0, 0, 0));
		// Put the horizontal box into a top-level labeled box.
		return Utilities.createLabeledComponents("Set cache data structure and maximum metrics to cache:",
				"(influences speed and memory of MST construction)", 0, false, horizBox);
	}
	
	/**
	 * Action listener to enable/disable frame depending on analyzer selected.
	 * 
	 * This method is invoked whenever the user changes the currently selected 
	 * analyzer in the {@link #analyzerList} combo box. This method enables or
	 * disables the spinner for setting the window/frame size.
	 * 
	 * @param e The action event associated with this method. This event is
	 * currently unused.
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		final boolean isTwoPassD2 = (analyzerList.getSelectedIndex() == 0);
		// The adaptive TwoPassD2 analyzer has been selected. 
		// The user cannot modify the frame size here.
		frameSize.setVisible(!isTwoPassD2);
		frameSizeAuto.setVisible(isTwoPassD2);
	}
	
	/**
	 * Obtain the currently set frame size for the analyzer.
	 * This method is currently used by the HeuristicsWizardPage class.
	 * 
	 * @return The currently set frame size value for the analyzer.
	 */
	protected Integer getFrameSize() {
		Number value = (Number) frameSize.getValue();
		return new Integer(value.intValue());
	}
	
	/**
	 * Obtain the currently selected analyzer.
	 * 
	 * This method is a convenience method that can be used by other
	 * wizard pages to determine the type of the wizard that the user
	 * has currently selected. Currently this method is used by the
	 * HeuristicsWizardPage to enable or disable heuristic configuration.
	 * 
	 * @return A predefined enumeration defining the type of analyzer that
	 * the user has currently chosen.
	 */
	protected FWAnalyzer.FWAnalyzerType getAnalyzerType() {
		return FWAnalyzer.FWAnalyzerType.values()[analyzerList.getSelectedIndex()];
	}
	
	/**
	 * This method is a convenience method that can be used to obtain
	 * a FWAnalyzer object containing information from this page. This
	 * method is used by the JobWizard to compose a complete set of
	 * entries once the wizard has successfully completed.
	 * 
	 * @return A FWAnalyzer object containing the necessary information
	 * about this analyzer.
	 */
	protected FWAnalyzer getAnalyzer() {
		// Array to translate cache type
		final String CacheTypes[] = {"heap", "mlist"};
		// Translate analyzer and cache selections to a suitable enums
		FWAnalyzer.FWAnalyzerType type = 
			FWAnalyzer.FWAnalyzerType.values()[analyzerList.getSelectedIndex()];
		// Set if it is a distance or similarity based approach
		boolean distance = !type.equals(FWAnalyzer.FWAnalyzerType.CLU);
		// Extract other analyzer related information.
		int windowSize   = ((Number) frameSize.getValue()).intValue();
		int wordSize     = 6;
		String cacheType = CacheTypes[cacheList.getSelectedIndex()];
		int size         = ((Number) cacheSize.getValue()).intValue();
		// Now create and return an analyzer object.
		FWAnalyzer info = new FWAnalyzer(type, distance, windowSize,
				wordSize, cacheType, size);
		return info;
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final JobWizard wizard;
	
	/**
	 * A combo-box to select the type of frame/word algorithm to
	 * be used for computing similarity metric.
	 */
	private JComboBox analyzerList;

	/**
	 * Formatted field to read the frame size to be used for
	 * analysis.
	 */
	private JSpinner frameSize;

	/**
	 * A simple JLabel that is used instead of the {@link #frameSize}
	 * spinner whenever the user selected TwoPassD2 to indicate that
	 * the frame size is automatically determined.
	 */
	private JLabel frameSizeAuto;
	
	/**
	 * This combo box is provided to permit the user to select a
	 * suitable word size to be used for analysis. Currently the
	 * word size is locked to 6.
	 */
	private JComboBox wordSize;

	/**
	 * A combo-box to select the type of caching that must be used
	 * to cache intermediate values to improve performance of
	 * the analyzer.
	 */
	private JComboBox cacheList;

	/**
	 * The size of the cache (number of entries maintained in the
	 * cache) to be used for analysis.
	 */
	private JSpinner cacheSize;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"<html>Configure the analyzer to be used to generate metrics<br>" +
		"(distance or similarity) for clustering.</html>";
	
	/**
	 * The list of analyzers that the user can choose from to compare
	 * two ESTs to determine similarity or distance. 
	 */
	private static final String ANALYZERS[] = {
		"Adaptive Two Pass D2 (safe & best choice)", 
		"Non-adaptive Two Pass D2 (for Sanger sequences)",
		"D2 Optimized (similar to wcd)", 
		"D2 Standard",
		"CLU (similarity metric)"
	};
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
