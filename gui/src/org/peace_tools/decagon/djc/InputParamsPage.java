package org.peace_tools.decagon.djc;

import java.awt.BorderLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;

import org.peace_tools.decagon.helpers.DADXHelper;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

/**
 * A wizard page to obtain inputs for parameters from the user.
 * 
 * This class serves as a wizard page in the DECAGON Job Creator.
 * This wizard page obtains inputs for various parameters
 * in the DADX file from the user. This class also populates
 * various DECAGON variables that can be used for conditional
 * processing of variables in the DADX file.
 * 
 */
public class InputParamsPage extends GenericWizardPage {
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to obtain information about
	 * the DADX file and the assembler pipeline.
	 */
	private final DecagonJobCreator wizard;
	
	/**
	 * The constructor.
	 * 
	 * The constructor initializes the wizard page and creates a couple
	 * of static GUI components providing static information to the
	 * user. The actual GUI components used to obtain inputs from the
	 * user are created later on (just before the page is actually
	 * displayed in the wizard) by the {@link #pageChanged(WizardDialog, int, int)}
	 * method.
	 * 
	 * @param wizard The DECAGON job creation wizard that logically owns
	 * this page.
	 */
	InputParamsPage(DecagonJobCreator wizard) {
		// Save DADX helper object for later use
		this.wizard = wizard;
		// Set the title for this page.
		setTitle("Parameter Inputs", "Set values for various pipeline parameters");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Set an informational label for now. The actual contents
		// for obtaining input will be created later
		JLabel infoLabel = new JLabel(INPUT_PARAMETERS_INFO, 
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		infoLabel.setBorder(new EmptyBorder(5, 5, 5, 5));
		add(infoLabel, BorderLayout.NORTH);
	}
	
	/**
	 * This method is called just before this page is to be displayed
	 * in the  wizard.
	 * 
	 * This method creates the actual GUI components to obtain
	 * inputs for various parameters from the user. First this method
	 * updates values for various DECAGON variables. Next this method
	 * calls the {@link DADXHelper#createGUI(JPanel)} method to create 
	 * the GUI inputs.
	 * 
	 * @param dialog The wizard dialog that logically owns this page.
	 * Currently this parameter is not used.
	 * 
	 * @param currPage The zero-based logical sequence index of this
	 * wizard page in the overall wizard. Currently this parameter
	 * is not used.
	 *  
	 * @param prevPage The zero-based logical sequence index of the 
	 * previous page from where the user was directed to this page.
	 * Currently this parameter is not used by this method.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Update pertinent DECAGON variables.
		wizard.updateDecagonVariables();
		// Create a sub-panel to obtain inputs (the following method
		// uses updated DECAGON variables for processing).
		JPanel inputSubPanel = wizard.getDADXHelper().createGUI(null);
		inputSubPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		JScrollPane jsp = new JScrollPane(inputSubPanel);
		// Add the various parameter-inputs to this panel.
		add(jsp, BorderLayout.CENTER);
	}
	
	/**
	 * A simple string description that is associated with clustering job
	 * entries created by this wizard. This string merely provides the
	 * user with a future reference as to why/how a output file was created. 
	 */
	private final String INPUT_PARAMETERS_INFO = 
		"<html>Set values for various parameters to be used for running<br/>" +
		"the various jobs (and processes) for the selected assembler<br/>" +
		"pipeline. Subsequent steps will obtain additional job-parameters." +
		"</html> ";
	
	/**
	 * A generated serialization GUID to keep the compiler happy. 
	 */
	private static final long serialVersionUID = -4262659158369134810L;

}