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
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import javax.swing.WindowConstants;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.PEACEProperties;
import org.peace_tools.decagon.PropertyKeys;
import org.peace_tools.generic.BorderWithComponent;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

/**
 * This class serves as an interactive page in a JobWizard.
 * This page permits the user to provide the information about the
 * local MST file where the data is to be stored. This wizard page 
 * also permits the user to select a server on which the job must 
 * be run. In addition, the wizard page permits the user to specify
 * the number of nodes and CPUs to be used for a parallel job.
 * 
 * Note that this wizard page checks to ensure that the MST file
 * does not exist, yet.
 */
public class MetasimConfigPage extends GenericWizardPage 
implements ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a text box and
	 * button to specify the MST file (this file must not yet exist
	 * as we are going to generate it from the job), a combo box to
	 * select the server to run the job (this server must have peace
	 * already installed on it), JSpinners for number of nodes and
	 * CPUs per node to run the job. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param helpInfo A piece of HTML fragment loaded from the
	 * DecagonConfigInfo.html file. This information is extracted
	 * from the data file by DecagonConfigurationManager, wrapped
	 * in a suitable GUI component, and passed-in to this method.
	 */
	public MetasimConfigPage(DecagonConfigurationManager dcm, JComponent helpInfo) {
		assert(dcm != null);
		// Setup the title(s) for this page and border
		setTitle("MetaSim", "Configure path to MultiSim");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Setup the browse button and text field where the user 
		// can browse and setup path to MultiSim install.
		metasimPath = new JTextField();
		browse      = new JButton("Browse...");
		browse.addActionListener(this);
		browse.setActionCommand("browse");
		JPanel entryPane = new JPanel(new BorderLayout(5, 0));
		entryPane.add(metasimPath, BorderLayout.CENTER);
		entryPane.add(browse, BorderLayout.EAST);
		// Since don't have much going on in this page, we can
		// use some more vertical space for the help information.
		Dimension helpSize = helpInfo.getPreferredSize();
		helpSize.height   += 50;
		helpInfo.setPreferredSize(helpSize);
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents("Set path to Metasim executable", 
				"<html>The path will be checked " +
				"when the next button is clicked</html>", 
				0, true, entryPane, Box.createVerticalStrut(5), helpInfo);

		// Create the top-level enable/disable check box
		enableMetaSim = new JCheckBox("Enable MetaSim", true);
		enableMetaSim.setFont(enableMetaSim.getFont().deriveFont(Font.BOLD));
		enableMetaSim.setActionCommand("enableMetaSim");
		enableMetaSim.addActionListener(this);
		// Setup a custom border with the check-box part of the border
		BorderWithComponent border = new BorderWithComponent(enableMetaSim, subPanel, 
				BorderFactory.createCompoundBorder(BorderFactory.createEtchedBorder(),
						BorderFactory.createEmptyBorder(5, 10, 5, 10)));
		// Set up the border to make things look good.
		subPanel.setBorder(border);
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
		// Setup initial default values for these controls
		setupInitalValues();
	}

	/**
	 * This method is invoked from the constructor to setup the default
	 * values for this page.
	 */
	private void setupInitalValues() {
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		// Obtain the enable/disable flag values
		final boolean enableFlag = props.getBooleanProperty(PropertyKeys.METASIM_ENABLED, true);
		enableInputs(enableFlag);
		enableMetaSim.setSelected(enableFlag);
		// Load path to metasim and the flags
		metasimPath.setText(props.getProperty(PropertyKeys.METASIM_PATH, ""));
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
		summary.append("<b>Use of MetaSim : <i>" + 
				(enableMetaSim.isSelected() ? "Enabled" : "Disabled") + "</i></b><br>");
		summary.append(Utilities.HTML_4SPACES);
		summary.append("Metasim executable path: " +
				(enableMetaSim.isSelected() ? metasimPath.getText() : "n/a") +
				"<br>");
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
		// Obtain the enable/disable flag values
		props.setProperty(PropertyKeys.METASIM_ENABLED, enableMetaSim.isSelected());
		// Load path to Metasim and the flags
		props.setProperty(PropertyKeys.METASIM_PATH, 
				(enableMetaSim.isSelected() ? metasimPath.getText() : ""));
	}
	
	/**
	 * Helper method to enable/disable controls in this panel. 
	 * 
	 * This is a helper method that is used to enable or disable the
	 * various pertinent inputs in this wizard page.
	 * 
	 * @param enableFlag If this flag is true, then the controls are
	 * enabled. Otherwise they are disabled.
	 */
	private void enableInputs(boolean enableFlag) {
		metasimPath.setEnabled(enableFlag);
		browse.setEnabled(enableFlag);
	}
	
	/**
	 * Method to handle clicking of "Browse" button. This method essentially
	 * enables and disables the various inputs depending on the check box.
	 * 
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		final String cmd = event.getActionCommand();
		if (browse.getActionCommand().equals(cmd)) {
			JFileChooser jfc = new JFileChooser();
			jfc.setDialogTitle("Choose target MST File directory");
			jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			if (jfc.showDialog(this, " OK ") == JFileChooser.APPROVE_OPTION) {
				// Copy the chosen file as the Metasim Path
				String wsPath = jfc.getSelectedFile().getAbsolutePath();
				metasimPath.setText(wsPath);
			}
		} else if (enableMetaSim.getActionCommand().equals(cmd)) {
			// Enable/disable inputs based on check-box status
			final boolean enableFlag = enableMetaSim.isSelected();
			enableInputs(enableFlag);
		}
	}


	
	/**
	 * This method validates the MetaSim path set in this page.
	 * 
	 * This method is invoked just before the user navigates to the
	 * adjacent page. This method checks to ensure that the path to
	 * MetaSim executable is is valid before user navigates to the next 
	 * page. The actual check is performed using an instance of the
	 * MetaSimRunner internal class.
	 * 
	 * @param dialog The wizard dialog that logically owns this page.
	 * 
	 * @param currPage The zero-based index number of the current page.
	 * 
	 * @param nextPage The zero-based index number of the page to which
	 * the user is requesting to navigate.
	 * 
	 * @return This method returns true if the user should be permitted
	 * to navigate from currPage to nextPage. 
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		if (!enableMetaSim.isSelected()) { 
			// Metasim is not selected. Nothing further to do.
			return true;
		}
		
		// Verify that metasim path has been correctly set.
		final String msPath = metasimPath.getText().trim();
		if (msPath.length() == 0) {
			// The path is not valid. Let the user know that.
			JOptionPane.showMessageDialog(this, METASIM_PATH_NOT_SET_MSG,
					"Metasim not usable", JOptionPane.WARNING_MESSAGE);
			// Do not proceed to the next page.
			return false;
		}
		MetaSimRunner msr = new MetaSimRunner(dialog);
		msr.start();
		// Let the wizard know whether to process to next page or not.
		return msr.isSuccess();
	}

	/**
	 * A helper class to run MetaSim in the background for validation.
	 * 
	 * This is a quick helper class that is used to run MetaSim in the
	 * background. MetaSim is run in the background to ensure that the
	 * GUI does not become unresponsive for a long time thereby confusing
	 * the user about the activity that is going on.  This class is
	 * created and used by {@link MetasimConfigPage#pageChanging(WizardDialog, int, int)}
	 * method.
	 *
	 */
	private class MetaSimRunner extends SwingWorker<Boolean, Boolean> {
		/**
		 * The only constructor for this class.
		 * 
		 * The constructor sets up the dialog box that is displayed by this
		 * class. The dialog is actually made visible by the caller by calling
		 * the {@link #start()} method.
		 * 
		 * @param parent The parent window for the temporary dialog box
		 * displayed by this method.
		 */
		public MetaSimRunner(WizardDialog parent) {
			progressBar = new JProgressBar();
			progressBar.setIndeterminate(true);
			progressBar.setStringPainted(true);
			progressBar.setString("Please wait...");
			// Create an informative label
			JLabel infoLabel = new JLabel(METASIM_TESTING_MSG, 
					Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
			// Set up an overall top-level container panel.
			waitPanel = new JDialog(parent, "Verifying MetaSim Path...", true);
			waitPanel.setAlwaysOnTop(true);
			waitPanel.setUndecorated(true);
			waitPanel.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
			waitPanel.setResizable(false);
			// Create a top-level content pane to contain various GUI components.
			JPanel contentPane = new JPanel(new BorderLayout(10, 10));
			contentPane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(Color.BLACK),
					BorderFactory.createCompoundBorder(BorderFactory.createRaisedBevelBorder(),
							BorderFactory.createEmptyBorder(10, 20, 10, 20))));
			// Setup sub-elements
			contentPane.add(infoLabel, BorderLayout.NORTH);
			contentPane.add(progressBar, BorderLayout.SOUTH);
			// Add content pane to the dialog
			waitPanel.setContentPane(contentPane);
		}

		/**
		 * Determine if MetaSim test was successful.
		 * 
		 * This method can be used once the runner has completed its task.
		 *  
		 * @return This method returns true if the test was successful;
		 * false otherwise.
		 */
		public boolean isSuccess() {
			boolean result = false;
			try {
				result = get();
			} catch (Exception e) {
				ProgrammerLog.log(e);
			}
			return result;
		}

		/**
		 * Method to start the actual work of the MetaSimRunner.
		 * 
		 * This method is a convenience method that must be used to 
		 * display the dialog (indicating work is underway) and
		 * actually start running MetaSim in the background.
		 */
		public void start() {
			// Finish layout.
			waitPanel.pack();
			waitPanel.setLocationRelativeTo(waitPanel.getParent()); 
			// Start background task.
			execute();
			// Make the wait dialog visible, blocking the caller.
			waitPanel.setVisible(true);
		}

		@Override
		protected Boolean doInBackground() throws Exception {
			// Try and run metasim using the supplied path and ensure
			// it is operating as expected.
			try {
				Process process     = new ProcessBuilder(metasimPath.getText(), "-h").start();
				// Read all the data into a single string from standard output
				// assuming standard error does not fill up.
				final String stdout = Utilities.readFullStream(process.getInputStream());
				// Convert the standard output to a string as well.
				final String stderr = Utilities.readFullStream(process.getErrorStream());
				// Wait for process to finish and get exit code.
				final int exitCode  = process.waitFor();
				if ((exitCode == 0) && (stdout.indexOf("MetaSim - Metagenome simulator") != -1) &&
						(stderr.length() == 0)) {
					// Successfully tested Metasim.
					success = Boolean.TRUE;
				}
				// Setup outputs from running metasim no matter what
				outputs = "-- Exit code: " + exitCode + " --\n-- Standard Output --\n" + 
				stdout + "\n-- Standard Error --\n" + stderr;
			} catch (Exception e) {
				// Save exception for processing in the done() method
				exp = e;
			}
			return success;
		}

		/** Method to report results at the end of the test.
		 * 
		 * This method is automatically invoked from the main Swing thread 
		 * once the {@link #doInBackground()} method has completed. This
		 * method hides the {@link #waitPanel} dialog and displays the
		 * results from running MetaSim.
		 */
		@Override
		protected void done() {
			// Hide the waitPanel.
			waitPanel.setVisible(false);
			waitPanel.dispose();
			// Display success/failure information.
			// Control drops here on successful case or when an exception occurs.
			final String mainMsg = (success ? "MetaSim path has been successfully set" :
			"Path to MetaSim executable is not valid.");
			String detailedMsg = "DETAILS FROM TEST:\n";
			if (exp != null) {
				detailedMsg += Utilities.toString(exp);
			} else {
				// Tag out outputs generated by running MetaSim
				detailedMsg += outputs;
			}
			JPanel fullInfo = Utilities.collapsedMessage(mainMsg, detailedMsg);
			JOptionPane.showMessageDialog(waitPanel.getParent(), fullInfo,
					"MetaSim Path Status", 
					(success ? JOptionPane.INFORMATION_MESSAGE : JOptionPane.ERROR_MESSAGE));
			super.done();
		}

		/**
		 * A simple dialog box requesting the user to wait while MetaSim is run.
		 * This dialog box provides the user with brief information about the
		 * fact that MetaSim is being run.
		 */
		private JDialog waitPanel;

		/**
		 * A roving progress bar to indicate to the user that some work
		 * is underway. This is more a visual cue and is not very meaningful.
		 * But it does provide some visual confirmation that work is being done.
		 */
		private JProgressBar progressBar;

		/**
		 * The outputs and results from running MetaSim. This value is initialized
		 * to empty string.  It is changed in the {@link #doInBackground()}
		 * method. 
		 */
		private String outputs;

		/**
		 * The exception that was generated (if any) when MetaSim was run.
		 * This is initialized to null and is later on
		 * changed if an exception occurs in the {@link #doInBackground()}
		 * method.
		 */
		private Exception exp = null;

		/**
		 * Final flag that indicates success or failure of MetaSim test.
		 * If this flag is false (default initial value) then MetaSim
		 * tests failed. Conversely, a true value indicates that the
		 * MetaSim test was successful and MetaSim path is good.
		 */
		private Boolean success = Boolean.FALSE;

		/**
		 * A static text displayed in the dialog displayed by this 
		 * class.
		 */
		private static final String METASIM_TESTING_MSG =
			"<html>" +
			"Attempting to validate path to MetaSim by trying to<br>" +
			"execute MetaSim and checking its default outputs.<br>" +
			"This test will take a few moments. Please wait..." +
			"</html>";
	}

	/**
	 * The top-level check-box that enables or disables the items in
	 * this wizard-page. If the check-box is checked, then other GUI
	 * items are enabled and permit the user to configure MetaSim path.
	 * Otherwise the use of MetaSim is disabled.
	 */
	private JCheckBox enableMetaSim;

	/**
	 * Field to read/display the install path for MetaSim
	 */
	private JTextField metasimPath;

	/**
	 * The browse button to be enabled the user to choose the
	 * directory where the file is to be stored.
	 */
	private JButton browse;

	/**
	 * A simple message that is displayed to the user when the path
	 * to Metasim is not set.
	 */
	private static final String METASIM_PATH_NOT_SET_MSG =
		"<html>" +
		"Path to Metasim executable is not set.<br/>" +
		"You need install Metasim (separately outside of PEACE)<br/>" +
		"and set the path to Metasim executable here";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1193172876744546755L;
}
