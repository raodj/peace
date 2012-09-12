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

package org.peace_tools.decagon.configMgr;

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.AbstractMenuHelper.HelperType;
import org.peace_tools.core.MainFrame;
import org.peace_tools.core.PEACEProperties;
import org.peace_tools.decagon.DecagonMenuHelper;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

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
public class DecagonConfigurationManager extends WizardDialog {
	/**
	 * The constructor for the configuration manager.
	 * 
	 * The constructor lays out the wizard and creates all the
	 * wizard pages. The wizard pages are created as templates
	 * and get populated with the necessary information just before
	 * the pages get displayed to the user.
	 * 
	 * @param parent The main frame that logically owns this wizard.
	 */
	public DecagonConfigurationManager(MainFrame parent) {
		super(parent);
		this.mainFrame = parent;
		setTitle("Decagon Configuration");
		setResizable(false);
		// Load inline help information from DecagonConfigInfo.html
		loadHelpInformation();
		// Set up the title image we want to use.
		setTitleBackground("images/peace_wizard_header.png", Color.white);
		// Set up the column image we want to use.
		setSequenceBackground("images/peace_wizard_column.png");
		// First setup the overview page.
		createOverview();
		// Create and add the Metasim path configuration page.
		mcp = new MetasimConfigPage(this, createHelpComponent(helpInformation[0]));
		addPage(mcp);
		// Create and add the shared storage configuration page.
		sscp = new SharedStorageConfigPage(this, createHelpComponent(helpInformation[1]));
		addPage(sscp);
		// Create and add the MySQL tunnel configuration page
		tcp = new MySQLTunnelConfigPage(this, createHelpComponent(helpInformation[2]));
		addPage(tcp);
		// Create and add the actual MySQL configuration page
		mscp = new MySQLConfigPage(this, tcp, createHelpComponent(helpInformation[3]));
		addPage(mscp);
		// Create the wizard page to configure the database
		msdbcp = new MySQLDBConfigPage(this, mscp, createHelpComponent(helpInformation[4]));
		addPage(msdbcp);
		// Create the verification page.
		VerifyWizardPage vwp = new VerifyWizardPage(this);
		addPage(vwp);
		// Finally create the page that commits the configuration.
		CommitConfigPage ccp = new CommitConfigPage(this);
		addPage(ccp);
	}

	/**
	 * Helper method to wrap in-line HTML information into a GUI component
	 * for use by web-pages.
	 * 
	 * This is a generic helper method that is called from the constructor
	 * for each wizard-page associated with this configuration manager.
	 * The help information to be presented to the user is passed-in as
	 * the parameter. This method wraps the help text in a GenericHTMLView
	 * component (that handles HTML hyper links, scrolling etc.) and returns
	 * the GUI component for use by various wizard pages.
	 * 
	 * @param helpInfo The help information to be wrapped in a suitable
	 * GUI component for presentation to the user.
	 * 
	 * @return The GUI component that can be readily used by various 
	 * wizard-pages for displaying in-line help to the user.
	 */
	private JPanel createHelpComponent(final String helpInfo) {
		JComponent helpView = mainFrame.createHTMLComponent(helpInfo, mainFrame.getBackground());
		// Create a panel to hold a title and the help information.
		JPanel helpPanel = new JPanel(new BorderLayout(0, 5));
		helpPanel.add(new JLabel("Helpful Information:"), BorderLayout.NORTH);
		helpPanel.add(helpView, BorderLayout.CENTER);
		return helpPanel;
	}
	/**
	 * Helper method to load informational messages from an multi-HTML
	 * file for use by various pages into {@link #helpInformation}.
	 * 
	 * This method is a helper method that is invoked from the constructor
	 * to load data from the <code>DecagonConfigInfo.html</code> file for
	 * use by various wizard pages. The <code>DecagonConfigInfo.html</code>
	 * file contains multiple segments of HTML, one for each of the wizard
	 * pages. The fragments of HTML are delimited by <code>&lt;html&gt;</code>
	 * and <code>&lt;/html&gt;</code>.  This method uses these delimiters
	 * to extract the HTML fragments into the {@link #helpInformation} array.
	 * If any error occurs when loading the help information this method 
	 * displays the error to the user.
	 */
	private void loadHelpInformation() {
		try {
			// Load the full file.
			final String htmlInfo = Utilities.readSmallTextFile("decagon/DecagonConfigInfo.html");
			// Extract each set of <html>..</html> segments separately.
			int startPos = 0, entryID = 0;
			while (startPos < htmlInfo.length()) {
				final int endPos = htmlInfo.indexOf("</html>", startPos);
				if (endPos != -1) {
					helpInformation[entryID] = htmlInfo.substring(startPos, endPos + 7);
					entryID++;
				}
				startPos += (endPos + 7);
			}
		} catch (Exception e) {
			// Log the exception.
			ProgrammerLog.log(e);
			// Display error message to the user.
			JPanel fullInfo = Utilities.collapsedMessage(INFO_LOAD_ERROR_MSG, Utilities.toString(e));
			JOptionPane.showMessageDialog(mainFrame, fullInfo,
					"Unable to load inline help", JOptionPane.ERROR_MESSAGE);
		} 
	}
	
	/**
	 * Returns summary information as a HTML document 
	 * 
	 * This method is invoked by the VerifyWizardPage
	 * 
	 * @return A HTML-sub-fragment providing information about the
	 * configuration to be committed by this class.
	 */
	protected String getSummary() {
		StringBuilder sb = new StringBuilder(1024);
		// Startup the HTML summary.
		sb.append("<html><b>DECAGON Configuration</b><br>");
		sb.append("The information cannot be edited.<br><br>");
		sb.append("<font size=\"-1\">");
		// Obtain summary information about Metasim configuration.
		sb.append(mcp.getSummary());
		sb.append("<br>");
		// Obtain summary information about shared storage configuration.
		sb.append(sscp.getSummary());
		sb.append("<br>");
		// Obtain summary information about SSH-tunnel for MySQL
		sb.append(tcp.getSummary());
		sb.append("<br>");
		// Obtain summary information about MySQL connection & database
		sb.append(mscp.getSummary());
		sb.append(msdbcp.getSummary());
		// Wrap up the HTML document.
		sb.append("<br></font></html>");
		return sb.toString();
	}
	
	/**
	 * Commit configuration to global properties.
	 * 
	 * This is a convenience method that can be used to have the
	 * various wizard pages commit their configuration to the
	 * global configuration properties.
	 * 
	 */
	protected void commitProperties() {
		mcp.commitProperties();
		sscp.commitProperties();
		tcp.commitProperties();
		mscp.commitProperties();
		msdbcp.commitProperties();
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
				"<html>Are you sure you want to exit?<br>" +
				"<b>Note: Your changes will not be saved</b></html>",
				"Discard Configuration Changes?", JOptionPane.YES_NO_OPTION);
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
	 * wizard. Currently this method has no specific task to 
	 * perform but is present as a place holder for future 
	 * enhancements. 
	 */
	@Override
	public void done(boolean success) {
		if (!success) {
			return;
		}
		// Enable/disable menu items in DECAGON menu now that the
		// configuration has been successfully updated.
		final DecagonMenuHelper dcm = (DecagonMenuHelper) mainFrame.getMenuHelper(HelperType.DECAGON_MENU);
		dcm.updateMenuItems();
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
	 * The main frame that logically owns this job. This value is
	 * used to create a job monitor if a job is successfully 
	 * submitted to run on a server.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user.
	 */
	private static final String OVERVIEW_MSG = "<html>" +
		"This wizard guides you through the process of configuring the following<br>" +
		"independent features of DECAGON: <ol>"+
		"<li>Path to MetaSim for generating synthetic reads</li>" +
		"<li>Location of shared storage for sharing data sets &amp; results</li>" +
		"<li>SSH-tunnel for MySQL Database connectivity</li>" +
		"<li>MySQL Database credentials</li></ol>" +
		"<b>Note</b>: You may configure just the features you need." +
		"</html>";

	/**
	 * The wizard configuration page that is used to configure the use
	 * of MetaSim.
	 */
	private final MetasimConfigPage mcp;
	
	/**
	 * The wizard page that is used to configure the use of 
	 * shared storage space on a given server accessible over 
	 * the network.
	 */
	private final SharedStorageConfigPage sscp;
	
	/**
	 * Wizard page that permits configuration of SSH-tunnel to access
	 * MySQL database.
	 */
	private final MySQLTunnelConfigPage tcp;
	
	/**
	 * Wizard page that permits configuration of MySQL access.
	 */
	private final MySQLConfigPage mscp;
	
	/**
	 * Wizard page that permits configuration of MySQL database to be
	 * used by DECAGON.
	 */
	private final MySQLDBConfigPage msdbcp;
	
	/**
	 * An array of HTML fragments that contain information that is displayed
	 * to the user in various wizard pages.
	 * 
	 * This array is initialized by the {@link #loadHelpInformation()} method
	 * from the <code>DecagonConfigInfo.html</code> file. The informational
	 * messages in the HTML fragments are passed to various wizard pages
	 * for their display and use.
	 */
	private final String helpInformation[] = new String[5];

	/**
	 * A static string that is displayed to the user if an error occurs
	 * when loading the DECAGON information from the <code>DecagonConfigInfo.html</code>
	 * file.
	 */
	private static final String INFO_LOAD_ERROR_MSG = 
		"<html>" +
		"Error loading help information for DECAGON configuration pages.<br>" +
		"Inline help information may be missing for some or all entries.<br>" +
		"Refer to online manual for details. <b>You may still proceed with<br>" + 
		"configuration settings despite this error.</b>" +
		"</html>";
	
	/**
	 * The generated serialization UID (need to keep the compiler happy) 
	 */
	private static final long serialVersionUID = 804993573751301886L;
	
	public static void main(String args[]) {
		// Turn off metal's use of bold fonts
		UIManager.put("swing.boldMetal", Boolean.FALSE);
		PEACEProperties.get().load(null, true);
		DecagonConfigurationManager dcm = new DecagonConfigurationManager(null);
		dcm.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		dcm.showWizard(null);
	}
}
