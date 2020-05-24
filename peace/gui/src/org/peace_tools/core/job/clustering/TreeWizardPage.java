package org.peace_tools.core.job.clustering;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;

/**
 * This class serves as an interactive page in a JobWizard. This page permits
 * the user to choose between MST and AST as the type of tree to construct
 * clusters. Tree type will influence the speed and of quality the clusters
 * made. The default tree type is MST. If the user select AST, a threshold needs
 * to be specified to determine how quickly nodes are added to the tree. If the
 * threshold is 0, AST essentially reduces to MST.
 */
public class TreeWizardPage extends GenericWizardPage implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components on this
	 * wizard page.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public TreeWizardPage() {
		// Setup the title(s) for this page and border
		setTitle("Tree Setup", "Configure the tree type to use for clustering");
		setBorder(new EmptyBorder(5, 5, 5, 5));

		// Create the informational label.
		JLabel info = new JLabel("<html>Configure the tree type to generate clusters.</html>",
				Utilities.getIcon("images/16x16/Information.png"), SwingConstants.LEFT);

		// Create the combo-box with list of analyzers.
		treeList = new JComboBox<>(new String[] { MST, AST });
		treeList.setBackground(Color.white);
		treeList.addActionListener(this);

		// Create panel with the combo box for the user to choose
		JComponent treeListBox = Utilities.createLabeledComponents("Select tree type for clustering:",
				"(Tree determines the way clusters are formed and can affect clusters quality)", 0, false, treeList);

		astThreshold = new JSpinner(new SpinnerNumberModel(0.5, 0.01, 1.0, 0.01));
		Utilities.adjustDimension(astThreshold, 50, 4);
		astThresholdPanel = Utilities.createLabeledComponents("Set AST threshold",
				"<html>(Set a threshold to aggressively use metrics to create AST. If the value is zero,<br>"
						+ "then AST essentially reduces to an MST.)</html>",
				0, false, astThreshold);
		astThresholdPanel.setVisible(false);

		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true, info, Box.createVerticalStrut(10),
				treeListBox, Box.createVerticalStrut(10), astThresholdPanel);
		// Set border to layout things nicely
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
	}

	/**
	 * Getter method to retrieve the selected type of tree.
	 * 
	 * @return the name of the tree, whether {@value #MST} or {@value #AST}
	 */
	public String getTreeSelection() {
		return String.valueOf(treeList.getSelectedItem());
	}

	/**
	 * Getter method to retrieve the threshold value if AST tree is selected.
	 * 
	 * @return string value of AST threshold. It returns string because later the
	 *         value is used to set a cmd parameter, which requires string value.
	 */
	public String getASTThreshold() {
		return String.valueOf(astThreshold.getValue());
	}

	/**
	 * Show or hide the threshold panel according to the user's selection. Only
	 * enable if the user selects AST.
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		boolean isAST = (treeList.getSelectedIndex() == 1);
		astThresholdPanel.setVisible(isAST);
	}

	/**
	 * A combo-box to select the type of tree.
	 */
	private final JComboBox<String> treeList;

	/**
	 * The threshold number input.
	 */
	private final JSpinner astThreshold;

	/**
	 * The panel that contains the threshold input. It will be shown or hidden
	 * according to the user's tree selection. Only shown when user selects AST.
	 */
	private final JPanel astThresholdPanel;

	/**
	 * A static value for Minimum Spanning Tree
	 */
	public static final String MST = "Minimum Spanning Tree";

	/**
	 * A static value for Approximate Spanning Tree
	 */
	public static final String AST = "Approximate Spanning Tree";

	/**
	 * Generated serial version UID
	 */
	private static final long serialVersionUID = 3741516715139116875L;

}
