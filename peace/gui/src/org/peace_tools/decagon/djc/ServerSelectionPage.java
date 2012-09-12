package org.peace_tools.decagon.djc;

import java.awt.BorderLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.job.ServerComboBox;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;

/**
 * A wizard page that permits the user to select the server
 * on which all jobs are to be run. 
 *
 * In DECAGON the software pipeline associated with a given
 * assembler is completely run on a single server. That is
 * different stages of the pipeline cannot be split across 
 * multiple servers. This wizard page permits the user to
 * select the server on which jobs are to be run. 
 * Subsequent wizard pages use the selected server to 
 * obtain additional information about the configurations to
 * be used for running each step of the pipeline.
 * 
 */
public class ServerSelectionPage extends GenericWizardPage {
	/**
	 * Constructor to configure this page to select a server for 
	 * all jobs for a given pipeline. 
	 * 
	 * This constructor is primarily used by the DecagonJobCreator 
	 * to setup a wizard page for permitting the user to select a server
	 * for running the pipeline. The constructor sets up the 
	 * various components on this wizard page.
	 * 
	 * @param wizard The DECAGON job creation wizard that logically 
	 * owns this page. This parameter cannot be null.
	 */
	public ServerSelectionPage(DecagonJobCreator wizard) {
		this.wizard       = wizard;
		// Setup the title(s) for this page and border
		setTitle("Select Server", "Set server to run assembler pipeline");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		
		// The list that displays servers that can be
		// used to run the assembler pipeline
		serverList   = new ServerComboBox();
		JPanel listPanel = Utilities.createLabeledComponents("Select server", 
				"This server will be used for all job(s) for this pipeline", 
				0, false, serverList);
		// Add the sub panel with all the GUI components to the wizard page.
		add(listPanel, BorderLayout.NORTH);
		
		// Create a simple label giving user some additional information.
		JLabel infoMsg = new JLabel(GENERAL_INFO_MSG, 
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		add(infoMsg, BorderLayout.CENTER);
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the list of server entries
	 * displayed in the combo-box. The user may select the server entry
	 * to be used for running jobs.
	 * 
	 * @param dialog The wizard dialog that logically owns this page.
	 * Currently this parameter is not used.
	 * 
	 *  @param currPage The zero-based logical sequence index of this
	 *  wizard page in the overall wizard.
	 *  
	 *  @param prevPage The zero-based logical sequence index of the 
	 *  previous page from where the user was directed to this page.
	 *  
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Let the server info update server list and associated information
		ServerComboBox.getServerList(false, false, serverList);
		// Ensure we have at least one.
		if (serverList.getItemCount() < 1) {
			// No valid server is available. The user cannot proceed further.
			wizard.setButtonStatus(-1, 0, -1);
		}
	}
	
	/**
	 * Obtain the server entry selected by the user (if any).
	 * 
	 * @return The server entry selected by the user. If there
	 * isn't a valid server selected this method will return null.
	 */
	protected Server getSelectedServer() {
		if (serverList.getSelectedIndex() == -1) {
			// No selections.
			return null;
		}
		// Obtain entry in the server object form.
		Server entry = (Server) serverList.getSelectedItem();
		return entry;
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final DecagonJobCreator wizard;
	
	/**
	 * The combo box that permits the user to select the server
	 * to be used for running the job.
	 */
	private final ServerComboBox serverList;

	/**
	 * A simple informational message that is displayed to the user giving
	 * the user some information about subsequent steps in this 
	 * wizard.
	 */
	private static final String GENERAL_INFO_MSG =
		"<html>" +
		"In the next page the operational status of various executables<br/>" +
		"(required for the pipeline) will be verified. Subsequent steps<br/>" +
		"will permit configuration additional properties for each job<br/>" +
		"such as: CPUs, memory, and run time." +
		"</html>";

	/**
	 * A generated GUID to keep the compiler happy
	 */
	private static final long serialVersionUID = 2804769325422036361L;
}
