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

package org.peace_tools.decagon.sdg;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.util.ArrayList;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.MainFrame;
import org.peace_tools.core.PEACEProperties;
import org.peace_tools.decagon.DecagonMenuHelper;
import org.peace_tools.decagon.helpers.ParameterSetGUIHelper;
import org.peace_tools.decagon.helpers.ParameterSetLoader;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.SeqTechType;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves as the top-level class for managing various
 * configuration options associated with DECAGON.  Note that each
 * of the configuration options is independent of the other and 
 * they are all optional. The user can mix-and-match these options
 * to suit their needs.  Currently, this wizard permits the user 
 * to configure the following DECAGON-specific features:
 * 
 * <ul>
 * 
 * <li>Setup the location where MetaSim is installed. The MetaSim
 * path is used by DECAGON to generate synthetic reads.</li>
 *   
 * <li>Setup the location of a shared storage where shared data sets
 * and source-transcripts are stored. This storage space is typically
 * used to share the information between various team members.</li>
 * 
 * <li>Setup access to MySQL data base that is used to store summary
 * information about tests. MySQL data base provides a convenient
 * mechanism to catalog, process, and report results.</li>
 * 
 * <li>Setup SSH-tunneling option for MySQL data base. This feature
 * is used when MySQL requires login to a remote server prior to
 * accessing the MySQL data base.</li>
 * 
 * </ul>
 * 
 * The aforementioned configuration options are handled by independent
 * wizard-pages that are implemented as separate wizard pages in
 * this package.
 * 
 */
public class SyntheticDataSetGenerator extends WizardDialog {
	/**
	 * Helper method to validate necessary preconditions prior to creating
	 * this wizard.
	 *
	 * The source genes to be used to generate a synthetic data set
	 * must be present in the workspace as a data set. Consequently, 
	 * This is a helper method that is used to ensure that there is 
	 * at least one data set associated with a FASTA file in the workspace.
	 * 
	 * @param parent The parent window for the wizard dialog box. This parameter
	 * cannot be null.
	 * 
	 * @return If the workspace has the necessary/valid data sets 
	 * then this method instantiates a valid wizard and returns
	 * the newly created wizard for further use. Otherwise this method
	 * displays a suitable error message and returns null.
	 */
	public static SyntheticDataSetGenerator create(MainFrame parent) {
		// Check to ensure we have a valid data set that can be used
		// as the source genes in the workspace.
		String[] msg = SyntheticDataSetGenerator.NO_DATASET_MSG; // Set invalid condition
		final ArrayList<DataSet> dataSets = Workspace.get().getDataSets();
		// Check to ensure we have at least one good data set entry
		for(DataSet ds: dataSets) {
			if (ds.isGood() && (ds.getFileType().equals(DataFileType.FASTA))) {
				msg = null; // We got a good data set
				break;      // no further checks needed.
			}
		}
		// Display error message if set
		if (msg != null) {
			WizardDialog.showMessage(parent, msg);
			return null;
		}
		// The following operations may take some time. So we change
		// cursor to let user know about potential delay
		SyntheticDataSetGenerator sdg = null;
		final Cursor prevCursor       = parent.getCursor();
		try {
			parent.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			// Load the set of parameters used to provide parameters to MetaSim
			// for generating sanger, illumina, and 454 sequences using a helper method.
			final String[] XMLFileList = {"decagon/MetaSimSanger.xml", 
					"decagon/MetaSim454.xml", "decagon/MetaSimIllumina.xml"};
			ParameterSetLoader psl = new ParameterSetLoader(XMLFileList);
			ArrayList<ParameterSetGUIHelper> paramList = psl.loadParameterSets(parent);
			if (paramList == null) {
				// An error occurred when loading parameters
				return null;
			}
			// Return new wizard only if msg is null indicating no errors.
			sdg = new SyntheticDataSetGenerator(parent, paramList);
		} finally {
			// Reset cursor back to the previous cursor
			parent.setCursor(prevCursor);
		}
		return sdg;
	}
	
	/**
	 * The constructor for the configuration manager.
	 * 
	 * <p>The constructor lays out the wizard and creates all the
	 * wizard pages. The wizard pages are created as templates
	 * and get populated with the necessary information just before
	 * the pages get displayed to the user.</p>
	 * 
	 * <p><b>Note: The SDG wizard must be created via a call to
	 * the {@link SyntheticDataSetGenerator#create(MainFrame)} method.</p>
	 * 
	 * @param parent The main frame that logically owns this wizard.
	 * 
	 * @param paramList The list of parameters to be used for generating
	 * sanger, 454, and illumina sequences respectively.
	 * 
	 * @see DecagonMenuHelper
	 * 
	 */
	private SyntheticDataSetGenerator(MainFrame parent, 
			ArrayList<ParameterSetGUIHelper> paramList) {
		super(parent);
		this.mainFrame = parent;
		this.paramsList = paramList;
		setTitle("Synthetic DataSet Generator");
		setResizable(false);
		setPreferredSize(new Dimension(600, 400));
		// Set up the title image we want to use.
		setTitleBackground("images/peace_wizard_header.png", Color.white);
		// Set up the column image we want to use.
		setSequenceBackground("images/peace_wizard_column.png");
		// First setup the overview page.
		createOverview();
		// Next create wizard page to select source gene data set and to
		// set a description for the synthetic data set to be generated.
		sgsp = new SourceGenesSelectionPage(this);
		addPage(sgsp);
		// Create the wizard page for setting up sanger parameters for Metasim
		sangerParams = new MetaSimParametersPage(this, paramsList.get(0),
				SeqTechType.SANGER, "Configure Sanger parameters");
		addPage(sangerParams);
		// Create the wizard page for setting up 454 parameters for Metasim
		roche454Params = new MetaSimParametersPage(this, paramsList.get(1),
				SeqTechType.R454, "Configure 454 parameters");
		addPage(roche454Params);
		// Create the wizard page for setting up sanger parameters for Metasim
		illuminaParams = new MetaSimParametersPage(this, paramsList.get(2),
				SeqTechType.ILLUMINA, "Configure Illumina parameters");
		addPage(illuminaParams);
		// Create the verification page.
		VerifyWizardPage vwp = new VerifyWizardPage(this);
		addPage(vwp);
		// Create the actual synthetic data set generation page.
		dsgp = new DataSetGenerationPage(this);
		addPage(dsgp);
	}
		
	/**
	 * Helper method invoked when user clicks cancel button.
	 * 
	 * This is a helper method that is overridden in this class.
	 * This method is invoked when the user clicks the cancel
	 * button in the wizard. This method is used to display a 
	 * confirmation dialog to ensure that the user really wants to
	 * exit out of the wizard.
	 * 
	 * @return This method returns true if the user wants to quit
	 * out of the wizard dialog.
	 */
	@Override
	protected boolean cancel() {
		// Check if user really want's to quit.
		int result = JOptionPane.showConfirmDialog(this,
				"<html><b>Are you sure you want to exit?</b><br>" +
				"A synthetic data set for evaluating assembers will<br>" +
				"<i>not</i> be generated.</html>",
				"Exit Wizard?", JOptionPane.YES_NO_OPTION);
		if (result == JOptionPane.NO_OPTION) {
			// The user does not want to quit.
			return false;
		}
		// Yes, the user wants to quit. Ensure no threads are left.
		return super.cancel();
	}
	
	/**
	 * Helper method to create the overview page. This method was
	 * introduced to keep the code clutter in the constructor to
	 * a bare minimum.
	 */
	private void createOverview() {
		JLabel message = new JLabel(OVERVIEW_MSG);
		Utilities.adjustFont(message, 0, 10, -1);
		GenericWizardPage overview = new GenericWizardPage();
		overview.add(message, BorderLayout.NORTH);
		overview.setTitle("Overview", "Overview of configuration options.");
		overview.setBorder(new EmptyBorder(20, 15, 10, 5));
		addPage(overview);
	}
	
	/**
	 * This method overrides the final notification method in this
	 * wizard. Currently this method saves the workspace if
	 * creation and addition of new synthetic data set was
	 * successfully completed.
	 * 
	 * @param success If this flag is true then the wizard completed
	 * its tasks successfully.
	 */
	@Override
	public void done(boolean success) {
		if (success) {
			// Save the workspace automatically.
			mainFrame.saveDelayedWorkspace();
		}
	}
	
	/**
	 * Obtain the instance of the main frame class that owns this wizard.
	 * 
	 * This method is currently used by the submit job wizard page 
	 * to create a job monitor once a job has been successfully
	 * submitted. 
	 * 
	 * @return The main frame that owns this wizard.
	 */
	protected MainFrame getMainFrame() {
		return mainFrame;
	}
	
	/**
	 * Obtain the full path to the FASTA file that contains source genes
	 * to be used to generate synthetic reads.
	 * 
	 * @return The path to the FASTA file that contains source
	 * genes to be used to generate synthetic reads.
	 */
	protected String getSourceGeneFile() {
		return sgsp.getSourceGeneFile();
	}

	/**
	 * Obtain the full path to the target FASTA file to contain the
	 * generated synthetic reads.
	 * 
	 * @return The path to the FASTA file that contains the generated
	 * synthetic reads.
	 */
	protected String getTargetGeneratedFile() {
		return sgsp.getTargetFASTAFile();
	}
	
	/**
	 * Obtain the description entered by the user for the 
	 * synthetic data set to be generated.
	 * 
	 * @return The description entered by the user.
	 */
	protected String getDescription() {
		return sgsp.getDescription();
	}
	
	/**
	 * Obtain estimate of number of synthetic reads generated for a given
	 * technology type.
	 * 
	 * @param techType The technology type for which the estimated
	 * number of synthetic reads is to be returned.
	 * 
	 * @return The number of reads generated for the given sequencing
	 * technology. If the technology is not enabled or if it is invalid,
	 * then this method returns -1.
	 */
	protected long getCount(SeqTechType techType) {
		if (SeqTechType.SANGER.equals(techType)) {
			return sangerParams.getCountSummary(true, null);
		} else if (SeqTechType.R454.equals(techType)) {
			return roche454Params.getCountSummary(false, null);
		} else if (SeqTechType.ILLUMINA.equals(techType)) {
			return illuminaParams.getCountSummary(false, null);
		}
		return -1;
	}
	
	/**
	 * Returns summary information as a HTML document 
	 * 
	 * This method is invoked by the VerifyWizardPage to obtain
	 * summary information to be displayed to the user for
	 * verification purposes. This method essentially collates
	 * information regarding sanger, 454, and illumina sequences
	 * from
	 * 
	 * @return A HTML-sub-fragment providing information about the
	 * configuration to be committed by this class.
	 */
	protected String getSummary() {
		StringBuilder sb = new StringBuilder(1024);
		// Startup the HTML summary.
		sb.append("<html>" + VERIFICATION_MSG);
		// First create summary information about number of sequences that
		// will be generated as a table to show columnar data.
		sb.append("<b>Sequences summary information</b><br/>");
		sb.append("<table><tr><td><i>Technology</i></td>" +
				"<td><i>Count</i></td><td><i>Avg. Read Len</i></td></tr>");
		// Get read-count summary information for each of the technologies.
		long totalReadCount = sangerParams.getCountSummary  (true, sb);
		totalReadCount     += roche454Params.getCountSummary(false, sb);
		totalReadCount     += illuminaParams.getCountSummary(false, sb);
		// Create row for total read count.
		sb.append("<tr><td><i>Total</i></td><td>~" + totalReadCount + 
				"</td><td></td></tr>");
		sb.append("</table>");
		// Obtain more detailed information about Metasim configurations for
		// each technology area.
		sb.append(sangerParams.getSummary());
		sb.append("<br/>");
		// Obtain summary information about shared storage configuration.
		sb.append(roche454Params.getSummary());
		sb.append("<br/>");
		// Obtain summary information about SSH-tunnel for MySQL
		sb.append(illuminaParams.getSummary());
		sb.append("<br/>");
		// Wrap up the HTML document.
		sb.append("<br/></html>");
		return sb.toString();
	}
	
	/**
	 * Getter method used by the various wizard pages.
	 * 
	 * @param techType The type of sequencing-technology for which this class 
	 * is managing parameters. This value must be one of <code>Sanger</code>,
	 * <code>454</code>, or <code>Illumina</code>. This value is also reflected
	 * as the prefixes for various parameter-names in the corresponding parameter
	 * configuration XML file.
	 * 
	 * @return The MetaSimParametersPage associated with the given sequencing
	 * technology.
	 */
	protected MetaSimParametersPage getParametersPage(final SeqTechType techType) {
		if (SeqTechType.SANGER.equals(techType)) {
			return sangerParams;
		} else if (SeqTechType.R454.equals(techType)) {
			return roche454Params;
		} else if (SeqTechType.ILLUMINA.equals(techType)) {
			return illuminaParams;
		}
		return null;
	}
	
	/**
	 * The main frame that logically owns this job. This value is
	 * used to create a job monitor if a job is successfully 
	 * submitted to run on a server.
	 */
	private final MainFrame mainFrame;

	/**
	 * The second page in this wizard that permits user to select source
	 * gene data set and provide description for the synthetic data set
	 * to be generated.
	 */
	private final SourceGenesSelectionPage sgsp;
	
	/**
	 * The wizard page that contains parameters associated with synthetic
	 * sanger-type reads in the generated data set.
	 */
	private final MetaSimParametersPage sangerParams;
	
	/**
	 * The wizard page that contains parameters associated with synthetic
	 * 454-type reads in the generated data set.
	 */
	private final MetaSimParametersPage roche454Params;
	
	/**
	 * The wizard page that contains parameters associated with synthetic
	 * illumina-type reads in the generated data set.
	 */
	private final MetaSimParametersPage illuminaParams;
	
	/**
	 * The parameters for Sanger, 454, and illumina sequences
	 * respectively. This array is set when the SDG is instantiated
	 * from the {@link #create(MainFrame)} method.
	 */
	private final ArrayList<ParameterSetGUIHelper> paramsList;
	
	/**
	 * The wizard page that performs the core operations required
	 * to generate synthetic data sets using the information
	 * supplied by the user.
	 */
	private final DataSetGenerationPage dsgp;
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user.
	 */
	private static final String OVERVIEW_MSG = "<html>" +
		"This wizard guides you through the process of generating<br>" +
		"a synthetic data set from a given set of genes: <ol>"+
		"<li>Select data set in workspace containing genes</li>" +
		"<li>Setup parameters for synthetic sanger data (if needed)</li>" +
		"<li>Setup parameters for synthetic 454 sequences (if needed)</li>" +
		"<li>Setup parameters for synthetic illumina data (if needed)</li>" +
		"<li>Verify the data set to be generated</li>" +
		"<li>Create the data set and add it to workspace</li></ol>" +
		"<b>Note:</b> This wizard does not expose all MetaSim options.<br>" +
		"You may edit the final command-line to add options as needed." +
		"</html>";
	
	/**
	 * A simple static message that is displayed to the user if
	 * a valid data set (that can be used as the source gene set for
	 * creating synthetic data set) is not found in the workspace.
	 */
	private static final String[] NO_DATASET_MSG = {
		"A dataset referring to a FASTA file must be present in<br/>" +
		"the workspace to be used as the source gene set for <br/>" + 
		"generating a synthetic data set.",
		"<b>Currently, the workspace does not have a dataset with<br/>" +
		"a FASTA file (MetaSim can work only with FASTA files).</b>",
		"Add the FASTA file with source genes to be used<br/>" +
		"to the workspace prior to launching this wizard.<br/>" +
		"You may add a new data set using the data set wizard."
	};
	
	/**
	 * A simple HTML-style message fragment that is added to the top
	 * of the summary information generated by the {@link #getSummary()}
	 * method in this class.
	 */
	private static final String VERIFICATION_MSG =
		"<font size=+1>Synthetic Dataset Summary</font><br/>" +
		"<i>The information cannot be edited.</i><br/><br/>" +
		"<table><tr><td  bgColor='#d0d0d0'>" +
		"You have provided the following information regarding " +
		"various parameters. Please verify the information. Next " +
		"the wizard will use the parameters to run MetaSim to generate " +
		"synthetic reads and assemble them into a data set in the workspace" +
		"</td</tr></table>" +
		"<br/>";
	
	/**
	 * The generated serialization UID (need to keep the compiler happy) 
	 */
	private static final long serialVersionUID = 804993573751301886L;
	
	public static void main(String args[]) {
		// Turn off metal's use of bold fonts
		UIManager.put("swing.boldMetal", Boolean.FALSE);
		PEACEProperties.get().load(null, true);
		SyntheticDataSetGenerator dcm = SyntheticDataSetGenerator.create(null);
		dcm.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		dcm.showWizard(null);
	}
}
