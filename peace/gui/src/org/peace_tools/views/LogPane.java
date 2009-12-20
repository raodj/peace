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
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.views;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.SwingConstants;

import org.peace_tools.generic.Log;
import org.peace_tools.generic.LogListener;
import org.peace_tools.generic.Utilities;

/**
 * <p>This is an abstract base class that provides the core features for
 * displaying log entries. The logging system has been developed using the
 * Model-View-Controller (MVC) pattern. This class represents the core 
 * component of the "View" that provides a graphical view of the logs
 * using the data stored in the "Model" (the Log class). The top-level
 * Frame that contains the views acts as the "Controller".</p> 
 * 
 * Specifically, this base class provides the following features:
 * 
 * <ul>
 * 
 * <li>It provides a generic place holder that is used to UserLog and
 * ProgrammerLog classes to display information, suitably formatted to fit the
 * needs of the two different types of logs.</li>
 * 
 * <li>It provides a mechanism for enabling and disabling saving the logs to a
 * given file.</li>
 * 
 * </ul>
 */
public abstract class LogPane extends JPanel 
implements ActionListener, LogListener {
	/**
	 * Constant to make the fonts in the controls in the tool bar a bit
	 * small because they look huge on Linux.
	 */
	private final int FONT_SIZE_CHANGE = -2;
	
	/**
	 * A constructor for the LogPane. It performs the initial layout of the
	 * common elements in this panel.
	 * 
	 * @param log The log from which entries are actually displayed by this
	 * log pane.
	 *    
	 * @param createLevelControl If this flag is true, then a control to change
	 * the log level is created in the tool bar. This feature is used 
	 * only by the user log.
	 */
	public LogPane(Log log, boolean createLevelControl) {
		// Setup our static field.
		this.log = log;
		// Perform the basic GUI configuration
		setMinimumSize(new Dimension(100, 100));
		setLayout(new BorderLayout());
		// Create the controls at the top of the log pane(s) that permit the
		// user to: set/change the file name. Control if data is being logged,
		// and set the log level (if the user has chosen) it. These
		// components are created via helper methods.
		JToolBar controlPanel = new JToolBar("Log ToolBar");
		controlPanel.setFloatable(false);
		// Add components to the control panel.
		addLogFileNameControls(controlPanel);
		if (createLevelControl) {
			controlPanel.addSeparator();
			addLogLevelChanger(controlPanel);
		}
		// Add tool bar to ourselves
		add(controlPanel, BorderLayout.NORTH);
		// The derived classes will add a suitable component to center.
		
		// Finally register this view with the actual log to obtain
		// updates whenever log information changes.
		log.addLogListener(this);
	}

	/**
	 * This is a refactored utility method that is used in the constructor to
	 * create the log file display area, the "Change..." button, and the toggle
	 * button to enable/disabling writing logs to the file.
	 * 
	 * @param toolBar The tool bar to which the controls are to be added.
	 * 
	 */
	private void addLogFileNameControls(JToolBar toolBar) {
		// First create the various controls associated with changing/setting
		// the log file name.
		fileNameDisplay = new JTextField(40); // displays file name.
		fileNameDisplay.setEditable(false);
		fileNameDisplay.setBackground(Color.white);
		Utilities.adjustFont(fileNameDisplay, FONT_SIZE_CHANGE, 8, -1);
		// Create toggle start/stop saving button
		saveLogsButton = new JToggleButton(Utilities
				.getIcon("images/16x16/SaveLog.png"));  
		saveLogsButton.setEnabled(false);
		saveLogsButton.setActionCommand(TOGGLE_LOG_SAVING);
		saveLogsButton.setToolTipText("Start/Stop writing logs to log file");
		Utilities.adjustFont(saveLogsButton, FONT_SIZE_CHANGE, 8, -1);
		// Create button to permit user to choose a log file name.
 		JButton changeButton = new JButton(Utilities
				.getIcon("images/16x16/ChangeLogFile.png"));
		changeButton.setActionCommand(CHANGE_LOG_FILE_NAME);
		changeButton.addActionListener(this);
		changeButton.setToolTipText("Set/Change Log file name");
		Utilities.adjustFont(changeButton, FONT_SIZE_CHANGE, 8, -1);
		// Now add the various controls to the tool bar.
		JLabel label = new JLabel("Log filename: ");
		Utilities.adjustFont(label, FONT_SIZE_CHANGE, 8, -1);
		toolBar.add(label);
		toolBar.add(fileNameDisplay);
		toolBar.addSeparator();
		toolBar.add(changeButton);
		toolBar.addSeparator();
		toolBar.add(saveLogsButton);
	}

	/**
	 * This is a utility method that is invoked just from the constructor to
	 * create a combo box that permits the user to change the log level. This
	 * field is created only for the UserLog. This method was introduced to 
	 * keep the code clutter in the constructor to a minimum.
	 * 
	 * @param toolBar The tool bar to which the controls are to be added.
	 * 
	 */
	private void addLogLevelChanger(JToolBar toolBar) {
		String[] logLevels = { "Informational (unimportant)",
				"Notice (mildy important)", "Warnings", "Errors" };
		logLevelList = new JComboBox(logLevels);
		logLevelList.setSelectedIndex(1);
		logLevelList.setActionCommand(CHANGE_LOG_LEVEL);
		logLevelList.addActionListener(this);
		Utilities.adjustFont(logLevelList, FONT_SIZE_CHANGE, 8, -1);
		logLevelList.setBackground(Color.white);
		// Create the container with suitable layout.
		toolBar.add(new JSeparator(SwingConstants.VERTICAL));
		toolBar.add(Box.createHorizontalStrut(5));
		JLabel label = new JLabel("Set logging level: ");
		Utilities.adjustFont(label, FONT_SIZE_CHANGE, 8, -1);
		toolBar.add(label);
		toolBar.add(logLevelList);
	}

	/**
	 * The event-callback handler for the "Choose..." button, the "SaveLogs"
	 * button, and the combo box. This is an aggregate handler for all of the
	 * three controls in this class.
	 * 
	 * @param event
	 *            The action event to be processed by this method.
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		// This is the action listener for the two buttons
		// and combo box for this log.
		if (CHANGE_LOG_LEVEL.equals(event.getActionCommand())) {
			// The user wants to change the logging level.
			JComboBox cb = (JComboBox) event.getSource();
			Log.LogLevel level = Log.LogLevel.values()[cb.getSelectedIndex()];
			// Set the log level in our log.
			log.setLevel(level);
		} else if (TOGGLE_LOG_SAVING.equals(event.getActionCommand())) {
			// The user wants to enable/disable writing log to file.
			log.toggleLogWriting();
		} else if (CHANGE_LOG_FILE_NAME.equals(event.getActionCommand())) {
			final JFileChooser fc = new JFileChooser();
			int returnVal = fc.showSaveDialog(this);
			fc.setDialogTitle("Select log file name");
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				// The user seems to have chosen a file. So use it.
				String fileName = fc.getSelectedFile().getAbsolutePath();
				if (!log.setFileName(fileName)) {
					// The file was not successfully set. Show an error message.
					JOptionPane.showMessageDialog(this, "Error writing logs to " +
							"the log file.\n(file: " + fileName + ")", 
							"Error writing logs",
							JOptionPane.WARNING_MESSAGE);
				}
			}
		}
	}

	/**
	 * This method implements the LogListener interface method to update
	 * the GUI to reflect the current state of the logger associated
	 * with this log pane. Refer to the documentation in the interface
	 * for additional details.
	 */
	@Override
	public void logChanged(boolean logsChanged, boolean fileNameChanged,
			boolean levelChanged) {
		// Update the default controls based on the information.
		fileNameDisplay.setText(log.getFileName() != null ? log.getFileName() : "<none set>");
		saveLogsButton.setEnabled(log.getFileName() != null);
		if (logLevelList != null) {
			logLevelList.setSelectedIndex(log.getLevel().ordinal());
		}
	}
	
	/**
	 * This instance variable holds a reference to the actual log class
	 * from which this log pane is displaying log entries. This value is
	 * set in the constructor and is never changed during the life time
	 * of this object.
	 */
	final protected Log log;
	
	/**
	 * This is an internal text field that is used to display the file to which
	 * the logs are being currently written. If the fileName is null, then this
	 * field displays "<none>" and the toggle button to save data to the file
	 * is disabled.
	 */
	private JTextField fileNameDisplay;

	/**
	 * A toggle button that can be in "on" or "off" state to indicate if logs
	 * are to be saved to the specified fileName (if it is not null) or not.
	 * This button is disabled until the user chooses a valid log file name.
	 * When in the "on" state, new logs are written to the supplied file.
	 */
	private JToggleButton saveLogsButton;

	/**
	 * A combo box to permit the user to set the level at which logs are to be
	 * saved by this logger. Typically this facility is used only by the UserLog
	 * and not by the programmer log.
	 */
	private JComboBox logLevelList = null;

	/**
	 * A string constant to ensure that the action command associated with a
	 * swing component and its action listener are consistently used without
	 * causing any issues. Just a better programming practice and a coding
	 * convention.
	 */
	private static final String CHANGE_LOG_FILE_NAME = "ChangeLogFileName";

	/**
	 * A string constant to ensure that the action command associated with a
	 * swing component and its action listener are consistently used without
	 * causing any issues. Just a better programming practice and a coding
	 * convention.
	 */
	private static final String TOGGLE_LOG_SAVING = "ToggleLogSaving";

	/**
	 * A string constant to ensure that the action command associated with a
	 * swing component and its action listener are consistently used without
	 * causing any issues. Just a better programming practice and a coding
	 * convention.
	 */
	private static final String CHANGE_LOG_LEVEL = "ChangeLogLevel";

	/**
	 * A serialization UID to keep compiler happy. 
	 */
	private static final long serialVersionUID = -2203502234634483701L;

	/** 
	 * This method overrides the default implementation in the base class
	 * to remove the log pane from the list of log listeners to receive
	 * updates.
	 */
	@Override
	protected void finalize() throws Throwable {
		// Remove ourselves from the list of interested listeners for
		// log updates.
		log.removeLogListener(this);
		super.finalize();
	}
}
