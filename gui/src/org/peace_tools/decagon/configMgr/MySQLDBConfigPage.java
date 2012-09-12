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
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.Statement;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.PEACEProperties;
import org.peace_tools.decagon.MySQLConnector;
import org.peace_tools.decagon.PropertyKeys;
import org.peace_tools.generic.BackgroundTask;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Server;

/**
 * Class to create (if needed) and verify MySQL database tables. 
 * 
 * This page occurs after the MySQL configuration page and basic
 * connection to MySQL has been verified. This page permits the
 * user to select a database to use or create a new database
 * for DECAGON to use. Once a database has been selected this 
 * class verifies that the database is valid and has the
 * necessary tables.
 */
public class MySQLDBConfigPage extends GenericWizardPage implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include subset of
	 * ServerInfoPanel object that permits the user to enter
	 * MySQL server configuration.
	 * 
	 * @param dcm The wizard that logically owns this page.
	 * 
	 * @param mySqlConfigPage The MySQL configuration page from where
	 * the necessary connection and details are to be obtained for
	 * verifying the database.
	 * 
	 * @param helpInfo A piece of HTML fragment loaded from the
	 * DecagonConfigInfo.html file. This information is extracted
	 * from the data file by DecagonConfigurationManager, wrapped
	 * in a suitable GUI component, and passed-in to this method.
	 */
	public MySQLDBConfigPage(DecagonConfigurationManager dcm,
			MySQLConfigPage mySqlConfigPage, JComponent helpInfo) {
		assert(dcm != null);
		this.dcm = dcm;
		this.mySqlConfigPage = mySqlConfigPage;
		// Setup the title(s) for this page and border
		setTitle("Database Setup", 
				"Select or Create MySQL database & tables");
		setBorder(new EmptyBorder(5, 5, 5, 5));

		// Create a combo-box that will contain list of databases that
		// the user can choose from.
		dbList = new JComboBox(new String[]{"Database List (to be filled)"});
		dbList.setEditable(true);
		// Create the "Refresh" and "Create DB" buttons (just icons)
		JButton refresh = Utilities.createButton("images/16x16/Refresh.png", 
				"", "Refresh", this, "Refresh list of databases displayed in list above", true);
		JButton createDB = Utilities.createButton("images/16x16/DatabaseAdd.png", 
				"", "Create", this, "Create a new database (you will prompted for name of DB)", true);
		// Create a button panel to hold buttons.
		JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 2));
		buttonPanel.add(refresh);
		buttonPanel.add(createDB);
		// Create wrapper panel to contain the various controls
		wrapperPanel = new JPanel(new BorderLayout(5, 5));
		wrapperPanel.add(dbList, BorderLayout.CENTER);
		wrapperPanel.add(buttonPanel, BorderLayout.EAST);
		
		// Create an informational label to provide summary information
		// to the user.
		infoLabel = new JLabel("Filled-in later");
		// Create a sub-panel to hold the info-label and wrapper panel.
		JPanel subPanel = new JPanel(new BorderLayout(0, 5));
		subPanel.add(infoLabel, BorderLayout.NORTH);
		subPanel.add(wrapperPanel, BorderLayout.SOUTH);
		// Create top panel to contain label (at top), wrapper panel
		// (in the middle) and the help information in the bottom.
		JPanel mainPanel = new JPanel(new BorderLayout(0, 5));
		mainPanel.add(subPanel, BorderLayout.NORTH);
		mainPanel.add(helpInfo, BorderLayout.CENTER);
		// Add the contents to this page
		add(mainPanel, BorderLayout.CENTER);
	}

	/**
	 * Helper method to enable/disable GUI components on this page.
	 * 
	 * This method is invoked from the {@link #pageChanged(WizardDialog, int, int)}
	 * method. This method enables or disables the various input GUI
	 * components on this page. In addition, it also updates the 
	 * {@link #infoLabel} with appropriate message. 
	 * 
	 *  @param enabled Flag to indicate if use of MySQL has been
	 *  enabled. If this flag is true then MySQL has been enabled and
	 *  various components must also be enabled to permit the user to
	 *  configure the database.
	 */
	private void updateGUIComponents(boolean enabled) {
		if (enabled) {
			// MySQL use is enabled. Show the user a detailed message.
			final Server mySqlSrvr = mySqlConfigPage.getServer();
			final String msg = String.format(MYSQL_DB_CHECK_INFO_MSG,
					mySqlSrvr.getName(), mySqlSrvr.getPort(), mySqlSrvr.getUserID());
			// wrap the msg using HTML tags and convert message to HTML.
			final String htmlMsg = "<html><font size=\"-1\">" +
				Utilities.wrapStringToHTML(msg, 65) + "</font></html>";
			infoLabel.setText(htmlMsg);
		} else {
			// MySQL use is disabled. Show the user a suitable message.
			infoLabel.setText(MYSQL_DB_NOT_USED_MSG);
		}
		// Enable/disable various components
		Utilities.setEnabled(wrapperPanel, enabled);
	}
	
	/**
	 * Helper method to load/refresh the list of databases displayed to
	 * the user.
	 * 
	 * This is a helper method used refresh the list of databases that
	 * are currently visible to the user. This method essentially runs
	 * a "show databases;" SQL command and displays the list of tables.
	 */
	private void refreshDBList() throws Exception {
		// Clear out any existing database entries.
		dbList.removeAllItems();
		// Ensure we have usable connection to database.
		if ((connector == null) || (connector.getConnection() == null)) {
			// Connection has not yet been setup!
			throw new IOException("Connection to MySQL database has not been setup successfully");
		}
		// Create a statement to obtain list of databases.
        Statement st = connector.getConnection().createStatement();
        ResultSet rs = st.executeQuery("show databases;");
        // Populate combo-box with list of databases
        while (rs.next()) {
        	final String dbName = rs.getString(1);
        	dbList.addItem(dbName);
        }
        // Close the statement
        st.close();
	}
	
	/**
	 * Helper method to create a table.
	 * 
	 * <p>This is a helper method that is invoked from the 
	 * {@link #createDB(String)} method to create a given table. This
	 * method builds a SQL statement using the supplied information
	 * and creates a table.<p>
	 * 
	 * <p><b>NOTE: The first column is assumed to be the primary key.</b></p>
	 * 
	 * @param statement The statement to be used to create the table.
	 * 
	 * @param tableName The name of the table to be created.
	 * 
	 * @param columnInfo A 2-D array (see {@link #DataSetTableStructure} 
	 * for example) that provides information about each
	 * column in the table. Each row must have the following 4 non-null
	 * columns:
	 * 
	 * <ol>
	 * 
	 * <li> The first entry is the column name. The first entry in the
	 * first row is assumed to be the primary key.</li>
	 * 
	 * <li> the second entry is MySQL data type (such as: MEDIUMINT, 
	 * MEDIUMTEXT, etc.)</li>
	 * 
	 * <li> The third entry is additional SQL attributes that are to be
	 * used when creating a table which the given column entry.</li>
	 * 
	 * <li>The fourth entry is a a brief comment about the contents of 
	 * the column and its purpose.</li>
	 * 
	 * </ol>
	 * 
	 * @param description The overall description of the contents of the
	 * table.
	 * 
	 * @param bt The background task object to which logs are to be generated.
	 * 
	 * @exception Exception This method does not catch any exceptions
	 * but simply exposes them.
	 */
	private void createTable(final Statement statement, final String tableName,
			final String[][] columnInfo, final String description, BackgroundTask bt) throws Exception {
		StringBuilder createStmt = new StringBuilder(1024);
		createStmt.append("CREATE TABLE ");
		createStmt.append(tableName);
		createStmt.append("(");
		// Append information about the columns
		for(int col = 0; (col < columnInfo.length); col++) {
			 // Add "," as needed after each column
			createStmt.append((col > 0) ? "',\n" : "");
			// Create statement fragment for each column
			createStmt.append(columnInfo[col][0]); createStmt.append(' ');
			createStmt.append(columnInfo[col][1]); createStmt.append(' ');
			createStmt.append(columnInfo[col][2]); createStmt.append(' ');
			if (col == 0) {
				// We assume first row is primary key
				createStmt.append("PRIMARY KEY ");
			}
			// Add brief description for column
			createStmt.append("COMMENT '");
			createStmt.append(columnInfo[col][3]);
		}
		// Add overall comment/description for the table.
		createStmt.append("')\n COMMENT='");
		createStmt.append(description);
		createStmt.append("';");
		// Now we have a complete SQL statement to create the table.
		final String createTableSQL = createStmt.toString();
		// Cut logs indicating progress
		bt.log("Attempting to create SQL Table " + tableName + 
				"via SQL statement:\n" + createTableSQL + "\n");
		statement.executeUpdate(createTableSQL);
		bt.log("SQL Table " + tableName + " created successfully.\n");
	}
	
	/**
	 * Method to create a DECAGON database.
	 * 
	 * This method is called from the {@link #createDatabase()} method 
	 * in a background thread. This method creates the database and then
	 * creates tables in the database.
	 * 
	 * @param dbName The name of the database to be created.
	 * 
	 * @param bt The background task object to which progress information
	 * and logs are to be generated.
	 * 
	 * @throws Exception This method exposes any exceptions that may occur
	 * during the various operations.
	 */
	private void createDB(final String dbName, BackgroundTask bt) throws Exception {
		bt.log("Checking connection to MySQL server...\n");
		if ((connector == null) || (connector.getConnection() == null)) {
			// Connection has not yet been setup!
			throw new IOException("Connection to MySQL server has not been setup successfully");
		}
		bt.log("Connection to MySQL server is available.\n");
		// Create the database.
		final String CreateDBStmt = "CREATE DATABASE " + dbName + ";";
		bt.log("Attempting to create database " + dbName + 
				" via SQL statement:\n\t" + CreateDBStmt + "\n");
		Statement st = connector.getConnection().createStatement();
		st.executeUpdate(CreateDBStmt);
		bt.log("Database created successfully.\n");
		bt.updateProgress();
		// Connect to the database to create the tables.
		bt.log("Creating new SQL connection to new database " + dbName + "...\n");
		MySQLConnector dbConnector = new MySQLConnector(connector, dbName);
		dbConnector.connect(bt.getDialog(), "Creating tables in new database: " + dbName);
		Statement dbSt = dbConnector.getConnection().createStatement();
		bt.log("Successfully obtained connection to the new database.\n");
		// Create the DataSet table.
		createTable(dbSt, "DataSet", DataSetTableStructure, 
				"Contains pertinent details about data sets", bt);
		bt.updateProgress();
		// Create the Job table.
		createTable(dbSt, "Job", JobTableStructure, 
				"Contains details about jobs/tests on genome-assemblers", bt);
		bt.updateProgress();
		// The entries were created successfully. Add the newly created
		// database as an option in the list of usable databases
		dbList.addItem(dbName);
		dbList.setSelectedItem(dbName);
	}
	
	/**
	 * Helper method to validate columns in a table.
	 * 
	 * <p>This is a helper method that is invoked from the 
	 * {@link #validateDB(String)} method to validate a given table. 
	 * This method runs a simple query and uses the column names 
	 * and SQL data types to validate column names and types.<p>
	 * 
	 * <p><b>NOTE: The first column is assumed to be the primary key.</b></p>
	 * 
	 * @param statement The statement to be used to run a simple query
	 * to validate columns in the table.
	 * 
	 * @param tableName The name of the table to be created.
	 * 
	 * @param columnInfo A 2-D array (see {@link #DataSetTableStructure} 
	 * for example) that provides information about each
	 * column in the table. Each row must have the following 4 non-null
	 * columns:
	 * 
	 * <ol>
	 * 
	 * <li> The first entry is the column name. The first entry in the
	 * first row is assumed to be the primary key.</li>
	 * 
	 * <li> the second entry is MySQL data type (such as: MEDIUMINT, 
	 * MEDIUMTEXT, etc.)</li>
	 * 
	 * <li> The third entry is additional SQL attributes that are to be
	 * used when creating a table which the given column entry.</li>
	 * 
	 * <li>The fourth entry is a a brief comment about the contents of 
	 * the column and its purpose.</li>
	 * 
	 * </ol>
	 * 
	 * @param bt The background task object to which logs are to be generated.
	 * 
	 * @exception Exception This method does not catch any exceptions
	 * but simply exposes them.
	 */
	private void validateTable(final Statement statement, final String tableName,
			final String[][] columnInfo, BackgroundTask bt) throws Exception {
		final String sqlStmt = "SELECT * FROM " + tableName + " LIMIT 1";
		bt.log("Verifying SQL table " + tableName + 
				" using statement:\n" + sqlStmt + "\n");
		final ResultSet rs   = statement.executeQuery(sqlStmt);
		final ResultSetMetaData rsmd = rs.getMetaData();
		for  (int col = 1; (col<= rsmd.getColumnCount()); col++){
			final String colName = rsmd.getColumnName(col);
			final String colType = rsmd.getColumnTypeName(col);
			final int    colSize = rsmd.getColumnDisplaySize(col);
			final String colFullType = ("VARCHAR".equals(colType) ? "VARCHAR(" + colSize + ")" : colType);
			bt.log("    Validating Column #" + col + ": Name='" + colName + 
				   "' and DataType='" + colFullType + "' ... ");
			if (!columnInfo[col - 1][0].equals(colName) || 
				(!columnInfo[col - 1][1].equals(colFullType) && 
				 !("BOOLEAN".equals(columnInfo[col - 1][1]) && "TINYINT".equals(colType)))) {
				throw new IOException("Column #" + col + " is inconsistent with expectation. " +
						"[expected: Name='" + columnInfo[col - 1][0] + "' with DataType='" +
						columnInfo[col - 1][1] + "' but got Name='" + colName + 
						"' with DataType='" + colFullType + "']"); 
			
			}
			bt.log("Good.\n");
		}
		bt.log("Table " + tableName + " validated successfully.\n");
	}

	/**
	 * Method to validate a DECAGON database.
	 * 
	 * This method is called from the {@link #validateDatabase()} method 
	 * in a background thread. This method validates the given database
	 * if the database is valid it switches the wizard to the next
	 * page.
	 * 
	 * @param dbName The name of the database to be validated.
	 * 
	 * @param bt The background task object to which progress information
	 * and logs are to be generated.
	 * 
	 * @throws Exception This method exposes any exceptions that may occur
	 * during the various operations.
	 */
	private void validateDB(final String dbName, BackgroundTask bt) throws Exception {
		// Connect to the database to create the tables.
		bt.log("Creating new SQL connection to new database " + dbName + "...\n");
		MySQLConnector dbConnector = new MySQLConnector(connector, dbName);
		dbConnector.connect(bt.getDialog(), "Validating tables in database: " + dbName);
		Statement dbSt = dbConnector.getConnection().createStatement();
		bt.log("Successfully obtained connection to the new database.\n");
		bt.updateProgress();
		// Validate the DataSet table.
		validateTable(dbSt, "DataSet", DataSetTableStructure, bt);
		bt.updateProgress();
		// Validate the Job table.
		validateTable(dbSt, "Job", JobTableStructure, bt); 
		bt.updateProgress();
		// Cut overall logs.
		bt.log("Database " + dbName + " is a valid DECAGON database.\n");
	}

	/**
	 * Method to check and enable various background 
	 * operations to verify database setup.
	 * 
	 * This method is invoked just before this wizard page is about
	 * to be displayed. This method checks to see if MySQL use has been
	 * enabled (by calling) and triggers background operations that
	 * verify the MySQL database tables. If the database does not exist
	 * then the user is prompted to create the database. 
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
		// Update informational label and enable/disable GUI components.
		updateGUIComponents(mySqlConfigPage.isMySQLEnabled());
		
		if (mySqlConfigPage.isMySQLEnabled()) {
			Exception exp = null; // Exception if any is generated.
			try {
				// Update our connection to the database.
				connector = mySqlConfigPage.getMySQLConnector(false);
				// Connect to MySQL as needed. This can throw an exception.
				if (!connector.connect(this, "Connecting to verify MySQL database for DECAGON")) {
					// Not connected to the database! The user has been notified about the error.
					connector = null;
				} else {
					// Connected successfully. Refresh list of databases
					refreshDBList();
					// Select the current database name (if any)
					final String currDbName = PEACEProperties.get().getProperty(PropertyKeys.MYSQL_DB_NAME);
					if (currDbName != null) {
						dbList.setSelectedItem(currDbName);
					}
				}
			} catch (Exception e) {
				// Error occurred. Can't proceed further.
				ProgrammerLog.log(e);
				exp = e;
				connector = null;
			} finally {
				if ((connector == null) || (exp != null)) {
					// Could not establish connection to the database. The user
					// has to go back to previous page.
					Object msg = (exp != null) ? Utilities.collapsedMessage(MYSQL_CONNECT_ERROR_MSG, 
							Utilities.toString(exp)) : MYSQL_CONNECT_ERROR_MSG;
					JOptionPane.showMessageDialog(this, msg, 
							"Can't connect to MySQL", JOptionPane.ERROR_MESSAGE);
					// Disable the next button.
					dcm.setButtonStatus(-1, 0, -1);
				}
			}
		}
		super.pageChanged(dialog, currPage, prevPage);
	}

	/**
	 * This method validates the database schema is as expected.
	 * 
	 * This method is invoked just before the user navigates to the
	 * adjacent page. This method checks to ensure that the database
	 * schema is as expected. The check is performed in a background
	 * thread. So this method always returns false (indicating not
	 * to change the page) and eventually if the schema is valid it
	 * automatically switches to the next page.
	 * 
	 * @param dialog The wizard dialog that owns this page. This parameter
	 * is currently unused.
	 * 
	 * @param currPage The zero-based logical index position of this page
	 * in the set of pages in the wizard. 
	 * 
	 * @param nextPage The zero-based logical index of the subsequent page
	 * to which the user is attempting to switch. If this page is less than
	 * the currPage (the user clicked on "Prev" button) then this method
	 * does not perform any checks and freely permits the user to go back.
	 * 
	 * @return This method always returns false preventing page change.
	 * However, upon successful database validation the page is changed
	 * programatically.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		if (!mySqlConfigPage.isMySQLEnabled()) {
			// The user has disabled MySQL. No need to check database.
			// On to the next page.
			return true;
		}
		// This page validates the MySQL database in the background and
		// eventually switches to the next page.
		final MySQLDBConfigPage wizPage = this;
		final String dbName             = (String) dbList.getSelectedItem();
		// Create a runnable background task as an anonymous class.
		BackgroundTask.UserTask validateDBTask = new BackgroundTask.UserTask() {	
			@Override
			public void run(BackgroundTask bTask) throws Exception {
				wizPage.validateDB((String) bTask.getParameter(), bTask);
			}
			@Override
			public void done(BackgroundTask bTask) {
				// Display the overall status to the user.
				bTask.showMessage("The database is valid DECAGON database",
						"The database is not a valid DECAGON database.<br>" +
						"Please select a different database or create a new database.");
			}
			@Override
			public void dialogClosed(BackgroundTask bTask) {
				// Update the button status appropriately
				dcm.setButtonStatus(1, 1, -1);
				// Switch the wizard page if successful.
				if (bTask.getResult() == null) {
					final int currPage = wizPage.dcm.getCurrentPage();
					wizPage.dcm.changePage(currPage, currPage + 1); 
				}
			}
		};
		
		// Create the background task object to perform database creation in 
		// the background.
		final JLabel infoLabel = new JLabel("<htmL>Please wait while the database and required " +
				"tables are validated.</html>", Utilities.getIcon("images/32x32/Information.png"),
						JLabel.LEFT);
		final String[] StepList = {"Connect to database " + dbName, "Validate DataSet Table",
				"Validate Job Table"};
		BackgroundTask bt = new BackgroundTask(validateDBTask, infoLabel, true, 
				StepList, true, false, false, 1, true, true, false, dbName);
		// Kick start the background operations
		bt.start(true, this, "Validating DECAGON Database: " + dbName,
				Utilities.getIcon("images/16x16/Database.png"));
		
		// Always return false for now here. The wizard page will be changed
		// programatically if the database validation is successful.
		// First disable the buttons
		dcm.setButtonStatus(0, 0, -1);
		return false;
	}

	/**
	 * Helper method to prompt and create a DECAGON database.
	 * 
	 * This is a helper method that is invoked from the {@link #actionPerformed(ActionEvent)}
	 * method to create a new DECAGON database.
	 */
	private void createDatabase() {
		final String dbName = JOptionPane.showInputDialog(this, DB_CREATE_MSG, 
				"Create new DECAGON Database", JOptionPane.QUESTION_MESSAGE);
		if (dbName == null) {
			// The user canceled the operation. 
			return;
		}
		// Create the database in a background thread to ensure GUI is responsive.
		final MySQLDBConfigPage wizPage = this;
		BackgroundTask.UserTask createDBTask = new BackgroundTask.UserTask() {	
			@Override
			public void run(BackgroundTask bTask) throws Exception {
				wizPage.createDB((String) bTask.getParameter(), bTask);
			}
			@Override
			public void done(BackgroundTask bTask) {
				// Display the overall status to the user.
				bTask.showMessage("The database was successfully created",
						"The database could not be created.");
			}
			@Override
			public void dialogClosed(BackgroundTask bTask) {}
		};
		
		// Create the background task object to perform database creation in 
		// the background.
		final JLabel infoLabel = new JLabel("<htmL>Please wait while the database and required " +
				"tables are created.</html>", Utilities.getIcon("images/32x32/Information.png"),
						JLabel.LEFT);
		final String[] StepList = {"Create database " + dbName, "Create DataSet Table",
				"Create Job Table"};
		BackgroundTask bt = new BackgroundTask(createDBTask, infoLabel, true, 
				StepList, true, false, false, 1, true, true, false, dbName);
		// Kick start the background operations
		bt.start(true, this, "Creating DECAGON Database: " + dbName, 
				Utilities.getIcon("images/16x16/Database.png"));
	}

	@Override
	public void actionPerformed(ActionEvent ae) {
		final String cmd = ae.getActionCommand();
		
		if ("Refresh".equals(cmd)) {
			// Refresh list of databases shown in the dbList combo-box
			// and report any errors we encounter.
			try {
				refreshDBList();
			} catch (Exception exp) {
				ProgrammerLog.log(exp);
				final String errMsg = "refresh database entries:</dt>" +
					"<dd>" + Utilities.wrapStringToHTML(exp.getMessage(), 60) +
					"</dd></dl></html>";
				final Object msg = Utilities.collapsedMessage(errMsg, Utilities.toString(exp));
				JOptionPane.showMessageDialog(this, msg, 
						"Unable to refresh database", JOptionPane.ERROR_MESSAGE);
			}
		} else if ("Create".equals(cmd)) {
			// Get helper method to create database.
			createDatabase();
		}
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
		if (mySqlConfigPage.isMySQLEnabled()) {
			return Utilities.HTML_4SPACES + "MySQL database: " +
			dbList.getSelectedItem() + "<br>";
		} else {
			return "";
		}
	}
	
	/**
	 * Convenience method to commit user-entered properties to the global
	 * properties. This method is invoked by the top-level Wizard once
	 * the user has verified the configuration.
	 */
	protected void commitProperties() {
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		if (mySqlConfigPage.isMySQLEnabled()) {
			props.setProperty(PropertyKeys.MYSQL_DB_NAME, (String) dbList.getSelectedItem());
		} else {
			props.removeProperty(PropertyKeys.MYSQL_DB_NAME);
		}
	}
	
	/**
	 * A informational label that is displayed at the top
	 * of this page to inform the user about the operation
	 * being performed. This label is setup in the
	 * {@link #updateControls(boolean)} method in this page.
	 */
	private final JLabel infoLabel;
	
	/**
	 * The configuration page that contains MySQL configuration.
	 * This information is used to connect to the MySQL database.
	 */
	private final MySQLConfigPage mySqlConfigPage;

	/**
	 * The combo-box that permits the user to either select
	 * a database or simply type-in the name of the database to be used
	 * by DECAGON. Entries are added to this GUI component by the  
	 * {@link #pageChanged(WizardDialog, int, int)} method.
	 */
	private final JComboBox dbList;
	
	/**
	 * The panel that is used to wrap the {@link #dbList},
	 * "Refresh" and "Create DB" buttons at the bottom of the
	 * list of databases. This panel is used to enable/disable 
	 * the GUI components in one shot. 
	 */
	private final JPanel wrapperPanel;
	
	/**
	 * The connector to be used to interact with the MySQL database.
	 * This connector is created each time this page is displayed
	 * in the {@link #pageChanged(WizardDialog, int, int)} method.
	 */
	private MySQLConnector connector;
	
	/**
	 * The configuration wizard that logically owns this wizard page.
	 * This reference is used to enable/disable next button 
	 * depending on the results on this page.
	 */
	private final DecagonConfigurationManager dcm;
	
	/**
	 * This is a simple message that is formatted and displayed at the 
	 * top of this page when MySQL use has been enabled and the database
	 * is to be validated. Note that this message is converted to an HTML
	 * message and wrapped appropriately by the {@link #updateGUIComponents(boolean)}
	 * method.
	 */
	private static final String MYSQL_DB_CHECK_INFO_MSG = 
		"The list of databases in MySQL instance on host " +
		"<i>%s (port: %d)</i> for user ID <i>%s</i> " +
		"are listed below. Select existing DECAGON database from " +
		"list or create a new DECAGON database (you will be prompted " +
		"for the name of the database)";
	
	/**
	 * This is a simple message that is displayed as is at the 
	 * top of this page when MySQL use has been disabled and no
	 * further operations are needed. This value is used in 
	 * the {@link #updateGUIComponents(boolean)} method.
	 */
	private static final String MYSQL_DB_NOT_USED_MSG = "<html>" +
		"Use of MySQL has been disabled (in the previous page).<br>" +
		"No further operations are needed on this page.<br>" +
		"You may now proceed further with this wizard." +
		"</html>";

	/**
	 * This is a simple message that is displayed when prompting the
	 * user to enter database name. This variable is used in the
	 * {@link #createDatabase()} method.
	 */
	private static final String DB_CREATE_MSG = "<html>" +
		"Enter name for the new DECAGON database:<br>" +
		"<font size=\"-2\">" + 
		"If database with same name already exists it will not be created" +
		"</font></html>";

	/**
	 * A simple message that is displayed to the user if connection
	 * to MySQL could not be established.
	 */
	private static final String MYSQL_CONNECT_ERROR_MSG = 
		"<html>" +
		"Unable to reconnect to MySQL database.<br>" +
		"Please navigate to the previous page and reconfigure<br>" +
		"connectivity to MySQL database." +
		"</html>";
	
	/**
	 * Definition for the columns in the DataSet table.
	 * Each column in the DataSet table is defined using a set of
	 * following 4 entries: 
	 * 
	 * <ol>
	 * 
	 * <li> The first entry is the column name. The first entry in the
	 * first row is assumed to be the primary key.</li>
	 * 
	 * <li> the second entry is MySQL data type (such as: MEDIUMINT, 
	 * MEDIUMTEXT, etc.)</li>
	 * 
	 * <li> The third entry is additional SQL attributes that are to be
	 * used when creating a table which the given column entry.</li>
	 * 
	 * <li>The fourth entry is a a brief comment about the contents of 
	 * the column and its purpose.</li>
	 */
	private static final String[][] DataSetTableStructure = {
		{"DataSetID", "INT", "NOT NULL AUTO_INCREMENT", 
		 "A unique identifier used to cross reference in Job table"},
		{"Host", "VARCHAR(256)", "NOT NULL", 
		 "Name/IP & port of machine with shared storage that contains the dataset"},
		{"Path", "VARCHAR(1204)", "NOT NULL", 
		 "The shared storage path on Host where dataset is located"},
		{"FileName", "VARCHAR(1024)", "NOT NULL", 
		 "The file name for the data set"},
		{"FileType", "VARCHAR(16)", "NOT NULL", 
		 "The storage file format like FASTA, SFF, SAM, BAM etc. (use NA if format is not known)"},
		{"Description", "VARCHAR(8192)", "", 
		 "User supplied description for this data set"},
		{"Source", "BOOLEAN", "NOT NULL", 
		 "Flag to indicate if this file contains source genes used to generate other datasets"},
		{"Generated", "BOOLEAN", "NOT NULL", 
		 "Flag to indicate if the file is a synthetic/generated file"},
		// Columns for Sanger sequences
		{"SangerCount", "INT", "NOT NULL", 
		 "Number of Sanger fragments in data set (zero if none)"},
		{"SangerAvgLen", "FLOAT", "NOT NULL", 
		 "Average length of sanger fragments"},
		{"SangerSD", "FLOAT", "NOT NULL", 
		 "Standard deviation in length of sanger sequences"},
		// Columns for Illumina sequences
		{"IlluminaCount", "INT", "NOT NULL", 
		 "Number of Illumina fragments in data set (zero if none)"},
		{"IlluminaAvgLen", "FLOAT", "NOT NULL", 
		 "Average length of Illumina fragments"},
		{"IlluminaSD", "FLOAT", "NOT NULL", 
		 "Standard deviation in length of Illumina sequences"},
		// Columns for Illumina sequences
		{"454Count", "INT", "NOT NULL", 
		 "Number of 454 fragments in data set (zero if none)"},
		{"454AvgLen", "FLOAT", "NOT NULL", 
		 "Average length of 454 fragments"},
		{"454SD", "FLOAT", "NOT NULL", 
		 "Standard deviation in length of 454 sequences"},
		// Additional information to verify consistency of data files if needed
		{"Timestamp", "TIMESTAMP", "NOT NULL", 
		 "Timestamp when the dataset was created (useful for verification)"},
		{"MD5CheckSum", "VARCHAR(128)", "", 
		 "MD5 checksum useful for verifying source dataset consistency & match"}
	};
	
	/**
	 * Definition for the columns in the Job table.
	 * Each column in the DataSet table is defined using a set of
	 * following 4 entries: 
	 * 
	 * <ol>
	 * 
	 * <li> The first entry is the column name. The first entry in the
	 * first row is assumed to be the primary key.</li>
	 * 
	 * <li> the second entry is MySQL data type (such as: MEDIUMINT, 
	 * MEDIUMTEXT, etc.)</li>
	 * 
	 * <li> The third entry is additional SQL attributes that are to be
	 * used when creating a table which the given column entry.</li>
	 * 
	 * <li>The fourth entry is a a brief comment about the contents of 
	 * the column and its purpose.</li>
	 */
	private static final String[][] JobTableStructure = {
		{"JobID", "INT", "NOT NULL AUTO_INCREMENT", 
		 "A unique identifier for each job in this table"},
		{"UserID", "VARCHAR(256)", "NOT NULL", 
		 "The ID of the user who added a given row",},
		{"Description", "VARCHAR(8192)", "", 
		 "User supplied description for the job"},
		{"LoadDate", "TIMESTAMP", "NOT NULL", 
		 "Timestamp when the row was added to the table"},
		{"SrcGeneID", "INT", "NOT NULL", 
		 "ID of row in DataSet table that contains source genes"},
		{"DataSetID", "INT", "NOT NULL", 
		 "ID of row in DataSet table that contains data used for job"},
		{"Assembler", "VARCHAR(256)", "NOT NULL", 
		 "Name and version of the genomic-assembler used for assembly"},
		{"CmdLine", "VARCHAR(8192)", "", 
		 "The full command line used to run assembler"},
		{"RunDate", "TIMESTAMP", "NOT NULL", 
		 "Timestamp when the job finished running"},
		{"RunHost", "VARCHAR(256)", "NOT NULL", 
		 "The name of the machine/cluster on which job was run"},
		{"CPUModel", "VARCHAR(256)", "",
		 "The CPU model of the machine/cluster if available"},
		{"NumProcesses", "SMALLINT", "NOT NULL", 
		 "The number of parallel processes/threads used for the job"},
		{"PeakMemory", "FLOAT", "NOT NULL", 
		 "The peak memory in (MiB) used by the job"},
		{"RunTime", "INT", "NOT NULL", 
		 "Wall-clock time taken for job in seconds"}, 
		{"AScore", "FLOAT", "NOT NULL", 
		 "The A-score for contigs vs. SrcGeneID entries generated by assembler"},
		{"N50", "FLOAT", "NOT NULL", 
		 "The N50 score generated for contigs"},
	};
	
	/**
	 * Generated serialization UID to keep compiler happy
	 */
	private static final long serialVersionUID = 1744039409477868365L;
}
