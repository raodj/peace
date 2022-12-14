package org.peace_tools.core.job.east;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.Param;
import org.peace_tools.workspace.Server;

/**
 * This wizard page is relatively straightforward wizard page 
 * that permits the user to set some of the basic parameters 
 * for an EAST job. The parameters include:
 * 
 * <ul>
 * 
 * <li><b>NUMOFLEVELS</b>: The number of levels used to verify
 *  a left end from a Minimum Spanning Tree (MST). Zero (0) is the 
 *  default, meaning doing until leaves of the MST. For two non-zero
 *  integers <i>a</i> and <i>b</i> (<i>a</i> &lt; <i>b</i>), setting 
 *  this parameter to <i>a<i> will lead to longer runtime (than setting
 *  it to <i>b</i>) while setting it to <i>b</i> will lead to better 
 *  quality in finding the "real" left end.</i>
 *  
 *  <li><b>USE_BOUNDED_NW</b>: A boolean flag indicating whether to use
 *  bounded version of Needleman-Wunsch (NW) algorithm to get NW score 
 *  when computing overlap distance. Use bounded version may reduce the
 *  accuracy, but save time.</li>
 *  
 *  <li><b>OUTPUT_ACE</b>: A boolean flag to indicate if the output
 *  of contigs (generated by EAST) are to be written in ACE file
 *  format.</li>
 *  
 * </ul>
 *
 */
public class EASTBasicParamsPage extends GenericWizardPage {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include appropriate GUI
	 * components to interactively obtain various basic parameters
	 * for EAST. Each GUI component is associated with a brief text
	 * message describing the purpose of each parameter.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param eastSWP The server wizard page that the user uses to select
	 * the server on which the EAST job is going to run. Information about
	 * the server selected by the user is used to determine the features
	 * available on the server (and the set of output file formats
	 * presented by this wizard).
	 */
	public EASTBasicParamsPage(EASTJobWizard wizard, ServerWizardPage eastSWP) {
		this.wizard  = wizard;
		this.eastSWP = eastSWP;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("EAST Parameters", 
				"Configure basic parameters for EAST");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the MST depth spinner.
		SpinnerNumberModel model = new SpinnerNumberModel(0, 0, 1000, 5);
		// Create and configure the spinner.
		numOfLevels = new JSpinner(model); 
		Utilities.adjustDimension(numOfLevels, -10, 4);
		// Create the check boxes for bounded NW and ACE output
		useBoundedNW = new JCheckBox();
		useBoundedNW.setSelected(true);
		// Create the list of fixed valid output file formats. 
		// SAM is added later on if server supports it.
		DataFileType[] FileFormats = {DataFileType.ACE, DataFileType.FASTA};
		contigOutputFormat = new JComboBox(FileFormats);
		
		// Create label indicating why SAM is not an valid output format
		noAMOSTools = new JLabel("<htmL>Server does not have AMOS tools and <br/>" +
				"therefore cannot create SAM output.</html>", 
				Utilities.getIcon("images/16x16/Warning.png"), JLabel.LEFT);
		// Add the various components along with suitable descriptive
		// labels to a suitable container panel
		JPanel subPanel = new JPanel(new GridBagLayout());
		add(subPanel, numOfLevels, NUMLEVELS_INFO_MSG, null);
		add(subPanel, useBoundedNW, BOUNDED_NW_INFO_MSG, null);
		add(subPanel, contigOutputFormat, CONTIG_OUTPUT_INFO_MSG, null);
		add(subPanel, Box.createHorizontalBox(), null, noAMOSTools);
		// Create the wizard-page-level informational message label.
		JLabel info = new JLabel(PAGE_INFO_MSG,
				Utilities.getIcon("images/32x32/Information.png"),
				JLabel.LEFT);
		// Finally create the top-level panel to contain an
		// informational label and the subPanel.
		JPanel panel = new JPanel(new BorderLayout(0, 0));
		panel.add(info, BorderLayout.NORTH);
		panel.add(subPanel, BorderLayout.CENTER);
		// Add the contents to this page
		add(panel, BorderLayout.CENTER);
	}
	
	/**
	 * A helper method to create the spinner ({@link #numOfLevels})
	 * that permits the user to set the number of levels to 
	 * explore in MST. This method has been introduced to 
	 * streamline the code in the constructor. This method
	 * is called only once from the constructor.
	 * 
	 * @param bag The panel that uses a GridBagLayout to which
	 * the left and right components are to be added.
	 * 
	 * @param left The component to be added to the left column
	 * of a grid bag layout. This component occupies as little
	 * space as possible.
	 * 
	 * @param labelText The label text to be added to the right 
	 * column of the grid bag layout. This component will
	 * occupy rest of the row. Note either labelText or label
	 * must be used (can't use both) 
	 * 
	 * @param label The label to be added to the right column of
	 * the grid bag layout. This component will occupy rest of
	 * the row. Note that either label or labelText must be used
	 * (both can't be used).
	 * 
	 */
	private void add(JPanel bag, JComponent left, String labelText, JLabel label) {
		// Create label to be placed to right of component if it is null.
		final JLabel right = (label != null) ? label : new JLabel(labelText);
		// Determine row based on current number of components
		final int row = bag.getComponentCount() / 2;
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.gridy              = row;
		gbc.gridx              = 0;
		gbc.insets             = new Insets(0, 5, 5, 0);
		gbc.anchor             = GridBagConstraints.NORTHEAST;
		// Add the left component with above constraints 
		bag.add(left, gbc);
		// Setup constraints for right component
		gbc.gridx              = 1;
		gbc.gridwidth          = GridBagConstraints.REMAINDER;
		gbc.fill               = GridBagConstraints.HORIZONTAL;
		// Add the right component with above constraints 
		bag.add(right, gbc);
	}
	
	/**
	 * Obtain the output file type currently chosen by the user.
	 * 
	 * This method is used by the OutputFilesWizardPage to provide
	 * appropriate default extension to the default files suggested
	 * by the wizard.
	 * 
	 * @return The output file type currently chosen by this method.
	 */
	protected DataFileType getContigOutputFormat() {
		return (DataFileType) contigOutputFormat.getSelectedItem();
	}
	
	/**
	 * Convenience method used by wizard to obtain the set of
	 * parameters configured on this page.
	 * 
	 * @return The set of EAST command line parameters configured
	 * by the user.
	 * 
	 */
	protected ArrayList<Param> getParamList() {
		ArrayList<Param> paramList = new ArrayList<Param>();
		final int numLevels = ((Number) numOfLevels.getValue()).intValue();
		paramList.add(new Param("-NUMOFLEVELS", "" + numLevels));
		if (useBoundedNW.isSelected()) {
			paramList.add(new Param("-USE_BOUNDED_NW", null));
		}
		if (!contigOutputFormat.getSelectedItem().equals(DataFileType.FASTA)) {
			paramList.add(new Param("-OUTPUT_ACE", null));
		}
		if (contigOutputFormat.getSelectedItem().equals(DataFileType.SAM)) {
			paramList.add(new Param("-CONVERT_TO_SAM", null));
		}
		// Add command line option to generate progress information
		paramList.add(new Param("--progress", "progress.dat"));
		return paramList;
	}
	
	/**
	 * Method to fill-in default file names and path when this page is
	 * displayed.
	 * 
	 * This method is called just before this page is to be displayed.
	 * This page merely adds/removes SAM as an optional output file
	 * format if the server has AMOS tools installed on it.
	 * 
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Convenient reference to the server on which EAST is going run.
		final Server eastServer = eastSWP.getServerInfoPanel().getSelectedServer();
		// Get the first output file format that is currently in the list.
		DataFileType firstFF = (DataFileType) contigOutputFormat.getItemAt(0);
		// Add SAM file format if the server supports it and our current
		// first entry is not SAM.
		if ((eastServer.hasAMOSTools()) && !DataFileType.SAM.equals(firstFF)) {
			contigOutputFormat.insertItemAt(DataFileType.SAM, 0);
		}
		// Remove SAM file format if the server DOES NOT support it and our current
		// first entry is not SAM.
		if ((!eastServer.hasAMOSTools()) && DataFileType.SAM.equals(firstFF)) {
			contigOutputFormat.remove(0);	
		}
		// Select first entry by default.
		contigOutputFormat.setSelectedIndex(0);
		// Show/hide warning label about missing AMOS tools
		noAMOSTools.setVisible(!eastServer.hasAMOSTools());
	}
	
	/**
	 * A boolean flag indicating whether to use  bounded 
	 * version of Needleman-Wunsch (NW) algorithm to get NW 
	 * score when computing overlap distance. Use bounded 
	 * version may reduce the accuracy, but save time.</li>
	 */
	private JCheckBox useBoundedNW;
	
	/**
	 * The set of output file formats in which the contigs from
	 * EAST are to be written.
	 */
	private JComboBox contigOutputFormat;
	
	/**
	 * The number of levels used to explore a Minimum Spanning
	 * Tree (MST) to determine a "real" left end of a contig. 
	 * Zero (0) is the default, meaning doing until leaves of
	 * the MST. For two non-zero integers <i>a</i> and <i>b</i>
	 * (<i>a</i> &lt; <i>b</i>), setting this parameter to <i>a<i>
	 * will lead to longer runtime (than setting it to <i>b</i>) 
	 * while setting it to <i>b</i> will lead to better quality 
	 * in finding the "real" left end.
	 */
	private JSpinner numOfLevels;

	/**
	 * A short message to indicate why SAM is not an available 
	 * output format as AMOS tools are not installed on the server.
	 * This label is shown/hidden by the {@link #pageChanged(WizardDialog, int, int)}
	 * method depending on the features of the selected server.
	 * The label is created by the constructor.
	 */
	private final JLabel noAMOSTools;
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final EASTJobWizard wizard;

	/**
	 * A reference to the server wizard page that contains the server
	 * selected by the user for running the EAST job. Information about
	 * the server selected by the user is used to determine the features
	 * available on the server (and the set of output file formats
	 * presented by this wizard).
	 */
	private final ServerWizardPage eastSWP;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page providing a brief overview of the
	 * purpose of this page.
	 */
	private static final String PAGE_INFO_MSG = 
		"<html>" +
		"Configure the parameters that play an important role in<br/>" +
		"the overall operations of the EAST assembler." +
		"</html>";
	
	/**
	 * A generic informational message that is displayed along with
	 * the MST depth numeric input. This message is meant to provide
	 * a brief description of the purpose of the input.
	 */
	private static final String NUMLEVELS_INFO_MSG = 
		"<html><b>Set MST-neighborhood exploration threshold</b><br/>" +
		"<font size=\"-2\">Larger values increases exploration space to provide<br/>"+
		"better contigs but increase runtime. Zero implies expore full MST." +
		"</font></html>";
	
	/**
	 * A generic informational message that is displayed along with
	 * the choice of using (or not using) bounded needleman-wunsch.
	 */
	private static final String BOUNDED_NW_INFO_MSG = 
		"<html><b>Use bounded Needleman-Wunsch (NW)</b><br/>" +
		"<font size=\"-2\">When checked, notably decreases assembly time at the<br/>"+
		"cost of some minor reduction in overall quality of contigs." +
		"</font></html>";
	
	/**
	 * A generic informational message that is displayed along with
	 * the choice of output file format for contigs generated by
	 * EAST.
	 */
	private static final String CONTIG_OUTPUT_INFO_MSG = 
		"<html><b>Select output file format for contigs</b><br/>" +
		"<font size=\"-2\">Select the file format in which the contigs<br/>" +
		"generated by the assembler are to be stored. Output file path<br/>" +
		"will be set in subsequent steps in this wizard." +
		"</font></html>";

	/**
	 * A generated serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1560174917966066396L;
}
