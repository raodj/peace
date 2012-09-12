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
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.PEACEProperties;
import org.peace_tools.decagon.MySQLConnector;
import org.peace_tools.decagon.PropertyKeys;
import org.peace_tools.generic.BorderWithComponent;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;

/**
 * Class to configure access to a MySQL data base.
 * 
 * This page permits the user to configure the necessary information
 * for accessing the MySQL database used by DECAGON for collating
 * results from various runs of genomic-assemblers.
 */
public class MySQLConfigPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include a set of input
	 * GUI components that permit the user to enter MySQL server 
	 * configuration.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param helpInfo A piece of HTML fragment loaded from the
	 * DecagonConfigInfo.html file. This information is extracted
	 * from the data file by DecagonConfigurationManager, wrapped
	 * in a suitable GUI component, and passed-in to this method.
	 */
	public MySQLConfigPage(DecagonConfigurationManager dcm,
			MySQLTunnelConfigPage tunnelConfig, JComponent helpInfo) {
		assert(dcm != null);
		this.tunnelConfig = tunnelConfig;
		// Setup the title(s) for this page and border
		setTitle("MySQL Configuration", 
				"Configure access to MySQL database");
		setBorder(new EmptyBorder(5, 5, 5, 5));

		// Create the custom component to supply the MySQL server
		// information.
		dbInfoPanel = createCustomServerInfoPanel();
		// Create an informational label to provide summary information
		// to the user.
		JLabel info = new JLabel(MYSQL_INFO_MSG, 
				Utilities.getIcon("images/16x16/Blank.png"), JLabel.LEFT);
		// Create label providing information about MySQL-tunnel
		tunnelLabel = new JLabel("filled in later", 
				Utilities.getIcon("images/16x16/Information.png"), JLabel.LEFT);
		tunnelLabel.setFont(tunnelLabel.getFont().deriveFont(Font.BOLD));
		// Put the labels within a panel.
		JPanel labelPanel = new JPanel(new BorderLayout(5, 1));
		labelPanel.add(info, BorderLayout.NORTH);
		labelPanel.add(tunnelLabel, BorderLayout.SOUTH);
		// Create top-level panel to contain the various fields
		JPanel mainPanel = new JPanel(new BorderLayout(0, 5));
		mainPanel.add(labelPanel, BorderLayout.NORTH);
		mainPanel.add(dbInfoPanel, BorderLayout.CENTER);
		mainPanel.add(helpInfo, BorderLayout.SOUTH);
		
		// Create the top-level enable/disable check box
		enableMySQL = new JCheckBox("Use MySQL Database", true);
		enableMySQL.setFont(enableMySQL.getFont().deriveFont(Font.BOLD));
		enableMySQL.setActionCommand("enableMySQL");
		enableMySQL.addActionListener(this);
		// Setup a custom border with the check-box part of the border
		BorderWithComponent border = new BorderWithComponent(enableMySQL, mainPanel, 
				BorderFactory.createCompoundBorder(BorderFactory.createEtchedBorder(),
						BorderFactory.createEmptyBorder(5, 10, 5, 10)));
		// Set up the border to make things look good.
		mainPanel.setBorder(border);
		// Add the contents to this page
		add(mainPanel, BorderLayout.CENTER);
		// Setup initial values
		setupInitialValues();
	}

	/**
	 * Helper method to create the GUI components to enter information 
	 * about MySQL database access.
	 * 
	 * This method is invoked only once from the constructor to create the
	 * GUI components used to enter information required to access MySQL. This
	 * method creates the host name, port, userID, and database name entries. 
	 * It wraps the GUI elements in a suitable panel and returns it back to 
	 * the caller.
	 *  
	 * @return A panel containing custom server information entry boxes.
	 */
	private JPanel createCustomServerInfoPanel() {
		// Create the various input fields that we are going to use.
		dbHost   = new JTextField(15);
		dbPort   = new JSpinner(new SpinnerNumberModel(3306, 1, 65535, 1));
		dbUserID = new JTextField(9);
		// Wrap the above components into a panel.
		JPanel dbInfoPanel = new JPanel();
		dbInfoPanel.setLayout(new BoxLayout(dbInfoPanel, BoxLayout.X_AXIS));
		dbInfoPanel.setAlignmentX(0);
		// Add components to the custom server information panel
		dbInfoPanel.add(Utilities.createLabeledComponents("MySQL Host:", null, 5, false, dbHost));
		dbInfoPanel.add(Box.createHorizontalStrut(5));
		dbInfoPanel.add(Utilities.createLabeledComponents("MySQL Port:", null, 5, false, dbPort));
		dbInfoPanel.add(Box.createHorizontalStrut(5));
		dbInfoPanel.add(Utilities.createLabeledComponents("MySQL User ID:", null, 5, false, dbUserID));
		// Return the customized panel back to the caller.
		return dbInfoPanel;
	}

	/**
	 * This method is invoked from the constructor to setup the default
	 * values for this page.
	 */
	private void setupInitialValues() {
		// Short cut to the current global properties
		final PEACEProperties props = PEACEProperties.get();
		// Setup current values from properties into GUI components
		dbHost.setText(props.getProperty(PropertyKeys.MYSQL_SERVER_NAME, "localhost"));
		final String portStr = props.getProperty(PropertyKeys.MYSQL_SERVER_PORT, "3306");
		dbPort.setValue(new Integer(portStr));
		dbUserID.setText(props.getProperty(PropertyKeys.MYSQL_SERVER_USERID, 
				System.getProperty("user.name")));
		// Obtain the enable/disable flag values
		final boolean enableFlag = props.getBooleanProperty(PropertyKeys.MYSQL_ENABLED, true);
		enableMySQL.setSelected(enableFlag);
		Utilities.setEnabled(dbInfoPanel, enableFlag);
	}
		
	/**
	 * Method to handle clicking of "Browse" button. This method essentially
	 * enables and disables the various inputs depending on the check box.
	 * 
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		// Enable/disable inputs.
		final boolean enabled = enableMySQL.isSelected();
		Utilities.setEnabled(dbInfoPanel, enabled);
	}

	/**
	 * Method to check and set SSH-tunnel usage message.
	 * 
	 * This method is invoked just before this wizard page is about
	 * to be displayed. This method checks and updates the {@link #tunnelLabel}
	 * information.
	 * 
	 * @param dialog The wizard dialog that logically owns this page.
	 * 
	 * @param currPage The zero-based index of the current page.
	 * 
	 * @param prevPage The zero-based index of the previous page from where
	 * the user navigated to this page.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		if (tunnelConfig.getTunnelServer() != null) {
			tunnelLabel.setText("A SSH-tunnel will be used");
		} else {
			tunnelLabel.setText("Using direct connection to the database");
		}
		super.pageChanged(dialog, currPage, prevPage);
	}
	
	/**
	 * Convenience method to get the server entry for MySQL database.
	 * 
	 * This is a convenience method that can be used to obtain the 
	 * a server entry that encapsulates the necessary information 
	 * to access the MySQL database.
	 * 
	 * @return A server entry that encapsulates the information
	 * entered by the user on this wizard page.
	 */
	protected Server getServer() {
		final int portNum = (Integer) dbPort.getValue();
		return new Server("MySQLServer", dbHost.getText(), 
				"The MySQL database server", dbUserID.getText(), 
				"", Server.getDefaultPollDuration(), true, portNum); 
				
	}
	
	/**
	 * Convenience method to obtain a MySQLConnector to connect to a 
	 * given MySQL database.
	 * 
	 * This method is a convenience method that is used by this class
	 * and the MySQLDBCheckPage to obtain a connection to MySQL server.
	 * 
	 * @param forceCreate If this parameter is true then the a new
	 * connector is always created. Otherwise if an existing connector
	 * is available it is returned.
	 * 
	 * @return A MySQLConnector object (that may not yet be connected to 
	 * the MySQL database) that can be used to connect to a MySQL database.
	 */
	protected MySQLConnector getMySQLConnector(final boolean forceCreate) {
		if ((sqlConnector != null) && (!forceCreate)) {
			return sqlConnector;
		}
		sqlConnector = new MySQLConnector(tunnelConfig.getTunnelServer(), getServer(), "");
		return sqlConnector;
	}
	
	/**
	 * Determine if MySQL has been enabled.
	 * 
	 * This is a convenience method that is used by 
	 * {@link MySQLDBConfigPage#pageChanged(WizardDialog, int, int)} to
	 * determine if the user has enabled the use of MySQL. 
	 * 
	 * @return This method returns true if MySQL has been enabled.
	 */
	protected boolean isMySQLEnabled() {
		return enableMySQL.isSelected();
	}
	
	/**
	 * This method validates the shared directory is valid.
	 * 
	 * This method is invoked just before the user navigates to the
	 * adjacent page. This method checks to ensure that the shared
	 * storage path is valid before user navigates to the next page.
	 * 
	 * @param dialog The wizard dialog that owns this page. This parameter
	 * is currently unused.
	 * 
	 * @param currPage The zero-based logical index position of this page
	 * in the set of pages in the wizard. This parameter is not used 
	 * by this method.
	 * 
	 * @param nextPage The zero-based logical index of the subsequent page
	 * to which the user is attempting to switch. If this page is less than
	 * the currPage (the user clicked on "Prev" button) then this method
	 * does not perform any checks and freely permits the user to go back.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		if (!enableMySQL.isSelected()) {
			// The user has disabled MySQL. No need to check connection.
			// On to the next page.
			return true;
		}
		try {
			// Ensure that the setup for accessing MySQL is operational.
			// For this we create temporary MySQLConnector object.
			final MySQLConnector connector = getMySQLConnector(true);
			if (!connector.connect(this, "Testing connectivity to MySQL")) {
				// The connection to MySQL was unsuccessful. Do not move to
				// next page.
				return false;
			}
			// Show successful connection message.
			JOptionPane.showMessageDialog(this, MYSQL_SUCCESS_MSG, 
					"MySQL connection successful", JOptionPane.INFORMATION_MESSAGE);
		} catch (Exception e) {
			ProgrammerLog.log(e);
			JPanel msg = Utilities.collapsedMessage(MYSQL_ERROR_MSG, Utilities.toString(e));
			JOptionPane.showMessageDialog(this, msg, 
					"Invalid MySQL Configuration", JOptionPane.ERROR_MESSAGE);			
			return false;
		}
		// MySQL connectivity is good.
		return true;
	}

	/**
	 * Returns summary information as a HTML fragment. 
	 * 
	 * This method is invoked by the top-level wizard for generating 
	 * summary information about the changes to be committed by this
	 * wizard page. 
	 * 
	 * @return A HTML-sub-fragment providing information about the
	 * configuration to be committed by this class.
	 */
	protected String getSummary() {
		StringBuilder summary = new StringBuilder(512);
		summary.append("<b>Use of MySQL : <i>" + 
				(enableMySQL.isSelected() ? "Enabled" : "Disabled") + "</i></b><br>");
		if (enableMySQL.isSelected()) {
			final String dbSrvrInfo = String.format(DB_SERVER_SUMMARY_INFO, dbHost.getText(), 
					dbPort.getValue().toString(), dbUserID.getText());
			summary.append(dbSrvrInfo);
		} else {
			summary.append(Utilities.HTML_4SPACES);
			summary.append("Additional details are not applicable");
		}
		return summary.toString();
	}
	
	/**
	 * Convenience method to commit user-entered properties to the global
	 * properties. This method is invoked by the top-level Wizard once
	 * the user has verified the configuration.
	 */
	protected void commitProperties() {
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		// Setup the enable/disable flag values
		final boolean enabled = enableMySQL.isSelected();
		props.setProperty(PropertyKeys.MYSQL_ENABLED, enabled);
		if (enabled) {
			props.setProperty(PropertyKeys.MYSQL_SERVER_NAME, dbHost.getText());
			props.setProperty(PropertyKeys.MYSQL_SERVER_PORT, dbPort.getValue().toString());
			props.setProperty(PropertyKeys.MYSQL_SERVER_USERID, dbUserID.getText());
		} else {
			props.removeProperty(PropertyKeys.MYSQL_SERVER_NAME);
			props.removeProperty(PropertyKeys.MYSQL_SERVER_PORT);
			props.removeProperty(PropertyKeys.MYSQL_SERVER_USERID);			
		}
	}
	
	/**
	 * The top-level check-box that enables or disables the items in
	 * this wizard-page. If the check-box is checked, then other GUI
	 * items are enabled and permit the user to configure MySQL 
	 * database information. 
	 */
	private JCheckBox enableMySQL;

	/**
	 * A simple label that is used to display whether a SSH-tunnel 
	 * will be used to establish MySQL connection.
	 */
	private final JLabel tunnelLabel;

	/**
	 * The configuration page that contains MySQL tunnel configuration.
	 * This information is used to test the connection to the 
	 * database.
	 */
	private final MySQLTunnelConfigPage tunnelConfig;
	
	/**
	 * Field to read/display the IP or name of the machine on which
	 * MySQL is running.
	 */
	private JTextField dbHost;

	/**
	 * Field to read/display the port number for MySQL. The default
	 * MySQL port is 3306. However, this port can be changed through
	 * MySQL configuration. Consequently, it is provided as a
	 * configurable value.
	 */
	private JSpinner dbPort;

	/**
	 * Field to read/display the user/login ID for MySQL. Note that
	 * the DB user ID is different than the userID to login to a
	 * given machine.
	 */
	private JTextField dbUserID;

	/**
	 * A top-level JPanel that wraps various database (DB) configuration
	 * GUI components. This panel wraps the {@link #dbHost}, {@link #dbPort},
	 * {@link #dbUserID}, and {@link #dbName} components. This information
	 * is used to conveniently enable/disable the four GUI elements. 
	 */
	private JPanel dbInfoPanel;
	
	/**
	 * The current connector being used by this class. This connector
	 * can be obtained by calling the {@link #getMySQLConnector()}
	 * method in this class.
	 */
	private MySQLConnector sqlConnector;
	
	/**
	 * A simple message that is displayed to the user at the top of this
	 * wizard page to provide abbreviated information. This message
	 * is formatted and then displayed.
	 */
	private static final String MYSQL_INFO_MSG = 
		"<html><font size=\"-2\">" +
		"Provide information to be used to access a MySQL<br>" +
		"database. The user ID is the MySQL DB user ID." +
		"</font></html>";
	
	/**
	 * A simple message that is displayed to the user if the SSH-Tunnel
	 * configuration failed to establish connection with the specified
	 * server.
	 */
	private static final String MYSQL_ERROR_MSG = 
		"<html>" +
		"The specified MySQL server configuration is invalid.<br>" +
		"Please verify that the following requirements are met:<ul>" +
		"<li>Ensure that server information (hostname/IP, port) etc. are valid</li>"+
		"<li>Ensure you have typed in the correct database UserID &amp; password.</li>" +
		"</ul>Please appropriately edit the necessary information."+
		"</html>";

	/**
	 * A simple message that is displayed to the user if the MySQL
	 * configuration worked as expected.
	 */
	private static final String MYSQL_SUCCESS_MSG = 
		"<html>" +
		"The specified MySQL server configuration is good.<br>" +
		"In the next page the actual MySQL database <br>" +
		"will be validated. If database does not exist you <br>"+ 
		"will need to create it in the next step." +
		"</html>";
	
	/**
	 * A HTML fragment that is formatted to fill-in the necessary information
	 * to generate summary in {@link #getSummary()} method.
	 */
	private static final String DB_SERVER_SUMMARY_INFO = 
		Utilities.HTML_4SPACES + "MySQL Host name: %s<br>" +
		Utilities.HTML_4SPACES + "MySQL Port: %s<br>" +
		Utilities.HTML_4SPACES + "MySQL User ID: %s<br>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8629449833492205394L;
}
