package org.peace_tools.core.job.clustering;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Heuristic;
import org.peace_tools.workspace.Param;

/**
 * This class serves as an interactive page in a JobWizard. This page permits
 * the user to choose whether or not to use primes heuristic. Primes heuristic
 * is a method to speed up the process of creating MST or AST by skipping D2
 * comparisons. If chose to enable, primes heuristic requires 4 options to be
 * set:
 * <ul>
 * <li>pri-heur-features: Number of features for primes heuristic</li>
 * <li>pri-heur-at: Prime value for A/T in primes heuristic</li>
 * <li>pri-heur-cg: Prime value for C/G in primes heuristic</li>
 * <li>pri-heur-topN: Restrict to top-n reads below threshold</li>
 * </ul>
 */
public class PrimesHeuristicWizardPage extends GenericWizardPage implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components on this
	 * wizard page.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public PrimesHeuristicWizardPage() {
		// Setup the title(s) for this page and border
		setTitle("Primes Heuristic Setup", "Configure the primes heuristic to use for clustering");
		setBorder(new EmptyBorder(5, 5, 5, 5));

		// Create the informational label.
		JLabel info = new JLabel("Primes heuristic is used to boost clustering performance.",
				Utilities.getIcon("images/16x16/Information.png"), SwingConstants.LEFT);

		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true, info, Box.createVerticalStrut(10),
				createOptionsPanel());
		// Set border to layout things nicely
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}

	/**
	 * A helper method to create the options panel. The user can choose values for 4
	 * different options that are required with primes heuristic.
	 */
	private JPanel createOptionsPanel() {
		// First create enable/disable checkbox with nice border.
		enablePrimesHeuristic = new JCheckBox("Enable primes heuristic for this job");
		JPanel bag = adjustCheckBox(enablePrimesHeuristic, "primesHeuristic");
		// Create and add the parameter.
		Box box = Box.createVerticalBox();

		// Create the input spinner.
		features = new JSpinner(new SpinnerNumberModel(4, 1, 100, 1));
		Utilities.adjustDimension(features, 50, 4);
		box.add(Utilities.createLabeledComponents("Number of features for primes heuristic", "", 0, false, features));
		features.setEnabled(false);

		at = new JSpinner(new SpinnerNumberModel(71, 2, 1000, 1));
		Utilities.adjustDimension(at, 50, 4);
		box.add(Utilities.createLabeledComponents("Prime value for A/T in primes heuristic", "", 0, false, at));
		at.setEnabled(false);

		at.addChangeListener(new CheckPrimeNumber(at, 71));

		cg = new JSpinner(new SpinnerNumberModel(113, 2, 1000, 1));
		Utilities.adjustDimension(cg, 50, 4);
		box.add(Utilities.createLabeledComponents("Prime value for C/G in primes heuristic", "", 0, false, cg));
		cg.setEnabled(false);

		cg.addChangeListener(new CheckPrimeNumber(cg, 113));

		topN = new JSpinner(new SpinnerNumberModel(50, 1, 100, 1));
		Utilities.adjustDimension(topN, 50, 4);
		box.add(Utilities.createLabeledComponents("Restrict to top-n reads below threshold", "", 0, false, topN));
		topN.setEnabled(false);

		// Now put all the information into a nice titled panel.
		JPanel bigBox = Utilities.createLabeledComponents(null, null, 0, false, bag, Box.createVerticalStrut(5), box);
		// Set border to make things look good.
		bigBox.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createTitledBorder("Primes Heuristic"),
				BorderFactory.createEmptyBorder(5, 5, 5, 5)));
		return bigBox;
	}

	/**
	 * Getter method to get the heuristic object constructed by this wizard page. If
	 * the user does not enable primes heuristic, null is returned.
	 * 
	 * @return if primes heuristic is enable, the constructed primes heuristic
	 *         object is returned, otherwise null.
	 */
	public Heuristic getPrimesHeuristic() {
		if (!enablePrimesHeuristic.isSelected()) {
			return null;
		} else if (primesHeuristic != null) {
			return primesHeuristic;
		}
		primesHeuristic = new Heuristic("primes");
		primesHeuristic.addParameter(new Param("pri-heur-features", features.getValue().toString()));
		primesHeuristic.addParameter(new Param("pri-heur-at", at.getValue().toString()));
		primesHeuristic.addParameter(new Param("pri-heur-cg", cg.getValue().toString()));
		primesHeuristic.addParameter(new Param("pri-heur-topN", topN.getValue().toString()));
		return primesHeuristic;
	}

	/**
	 * Helper method to correctly format, enhance, and layout a given check box.
	 * 
	 * @param cb  The check box to be formatted.
	 * @param cmd The action command to be set.
	 * @return A JPanel containing the check box with additional decorations.
	 */
	private JPanel adjustCheckBox(JCheckBox cb, String cmd) {
		cb.setBackground(Color.white);
		cb.setSelected(false);
		cb.setActionCommand(cmd);
		cb.setAlignmentX(0);
		cb.addActionListener(this);
		JPanel bag = new JPanel(new BorderLayout(0, 0));
		bag.setBackground(Color.white);
		bag.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEtchedBorder(),
				BorderFactory.createEmptyBorder(1, 2, 1, 2)));
		bag.add(cb, BorderLayout.NORTH);
		bag.setAlignmentX(0);
		return bag;
	}

	/**
	 * Enable or disable the options inputs as the user chooses to use primes
	 * heuristic or not.
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		features.setEnabled(enablePrimesHeuristic.isSelected());
		at.setEnabled(enablePrimesHeuristic.isSelected());
		cg.setEnabled(enablePrimesHeuristic.isSelected());
		topN.setEnabled(enablePrimesHeuristic.isSelected());
	}

	/**
	 * A class that implements ChangeListener to be added to {@link #at} and
	 * {@link cg}. A/T and C/G values are required to be prime numbers, therefore
	 * this listener validates the user's input for primes only, otherwise change it
	 * to the default value.
	 */
	class CheckPrimeNumber implements ChangeListener {
		private JSpinner input;
		private int defaultValue;

		public CheckPrimeNumber(JSpinner input, int defaultValue) {
			this.input = input;
			this.defaultValue = defaultValue;
		}

		@Override
		public void stateChanged(ChangeEvent arg0) {
			SpinnerNumberModel model = (SpinnerNumberModel) input.getModel();
			if (!isPrime((int) model.getValue())) {
				model.setValue(defaultValue);
			}
		}

		private boolean isPrime(int n) {
			if (n == 2) {
				return true;
			}
			if (n % 2 == 0) {
				return false;
			}
			for (int i = 3; i <= Math.sqrt(n); i += 2) {
				if (n % i == 0) {
					return false;
				}
			}
			return true;
		}
	}

	/**
	 * A checkbox to allow user to opt in for primes heuristic.
	 */
	private JCheckBox enablePrimesHeuristic;

	/**
	 * An input for the number of features for primes heuristic
	 */
	private JSpinner features;

	/**
	 * An input for prime value for A/T in primes heuristic
	 */
	private JSpinner at;

	/**
	 * An input for prime value for C/G in primes heuristic
	 */
	private JSpinner cg;

	/**
	 * An input for top-n reads restriction below threshold
	 */
	private JSpinner topN;

	/**
	 * The Heuristic object that is constructed in this class.
	 */
	private Heuristic primesHeuristic;

	/**
	 * Generated serial version UID
	 */
	private static final long serialVersionUID = 8028090850092763244L;
}
