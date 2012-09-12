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

package org.peace_tools.core;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.List;
import java.util.Scanner;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.UIManager;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;

import org.peace_tools.decagon.helpers.DADXHelper;
import org.peace_tools.generic.CustomPanel;
import org.peace_tools.generic.Log.LogLevel;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Workspace;

/**
 * This class provides a pretty sophisticated dialog box that is
 * displayed initially to the user to select a work space. A work space
 * is simply a directory in which PEACE GUI stores all the files that
 * are pertinent to a given GUI session. This class is used in the
 * following manner:
 * 
 * <ol>
 * 
 * <li>The main() method creates an instance of this class and calls
 * the createTabs() method to create all the necessary tabs in this
 * dialog box. If the tab creation fails, then the GUI exits. </li>
 * 
 * <li>If this is the first time PEACE is being launched (that is
 * the default PEACE directory does not exist) then this method
 * ensures that the user accepts the license agreement, creates the
 * default working directory for PEACE.</li>
 * 
 * <li>When a user selects a valid work space this class performs
 * the following tasks:
 * <ol>
 * <li>It hides the "OK" and "Cancel" buttons and displays a progress
 * bar to indicate progress of loading the necessary files. The actual
 * loading of files is performed using a SwingWorker thread that is
 * implemented by the BackgroundLoader inner class.</li>
 * 
 * <li>It interfaces with the DECAGON modules (DADXHelper) to load the
 * necessary assembler description files used by DECAGON for further
 * processing. As the DADX files are loaded the progress bar is
 * suitably updated.</li> 
 * 
 * <li>The work space instance variable, launches the MainFrame
 * and the dialog winds up.</li>
 * </ol>
 * </li>
 *  
 * </ol>
 *
 */
public class WorkspaceChooser extends JDialog implements ActionListener {
	/**
	 * Flag to indicate if this is the first time ever the user is running
	 * PEACE GUI. In this situation the user is forced to accept our
	 * license before proceeding further.
	 */
	private boolean firstLaunch = true;

	/**
	 * Flag that tracks if the first time welcome message has already
	 * been displayed or now.
	 */
	private boolean showFirstTimeMsg = true;
	
	/**
	 * The set of tabs that are displayed in this dialog for the user
	 * to work with.
	 */
	private JTabbedPane tabs;

	/**
	 * The list of work spaces that the user can choose from or edit
	 * to start up the GUI session.
	 */
	private JComboBox wsList;

	/**
	 * The selected work space that is set after the user successfully
	 * chooses a work space.
	 */
	private String workspace = null;
	
	/**
	 * This is a progress bar that is displayed instead of the "OK" and "Cancel"
	 * buttons after the user presses the OK button in this dialog. The 
	 * progress bar is used display the progress of loading various DECAGON
	 * data files and finally the actual workspace file. These operations
	 * can take several seconds and the visual indication to the user is
	 * critical in this situation. The actual object is created in 
	 * the {@link #createButtonProgBarPanel()}.
	 */
	private JProgressBar progBar = null;
	
	/**
	 * A simple label that is used to display some short additional information
	 * to the user. This label is displayed just above the {@link #progBarInfo}.
	 * The actual object is created in the {@link #createButtonProgBarPanel()}
	 * method.
	 */
	private JLabel progBarInfo = null;
	
	/**
	 * This panel manager is used to display either buttons or the progress
	 * bar. This panel and its layout manager is created in the 
	 * {@link #createButtonProgBarPanel()}
	 * and is used to switch from buttons to progress bar in the
	 * {@link #okAction()} method.
	 * 
	 */
	private JPanel btnProgPanel = null;
	
	/**
	 * An internal helper class to load various files in the background.
	 * 
	 * This helper class has been introduced to enable loading various
	 * DECAGON and PEACE workspace XML files in the background. The files
	 * are loaded in the background to ensure that the GUI is responsive
	 * and does not appear to have locked up during this slightly longer
	 * operation. This class is created from the {@link WorkspaceChooser#okAction()}
	 * method.
	 */
	private class BackgroundLoader extends SwingWorker<String, String> {
		/**
		 * The workspace directory path from where the workspace data is to
		 * be loaded.
		 */
		private final File wsPath;
		/**
		 * Flag to indicate if the directory path is empty and a default
		 * PEACE workspace file is to be created.
		 */
		private final boolean createDefWorkspace;
		
		/**
		 * The constructor for this inner class.
		 * 
		 * The constructor merely saves the parameters in instance variables
		 * to be used in various methods in this class.
		 * 
		 * @param wsPath The workspace directory path from where the workspace data is to
		 * be loaded.
		 * @param createDefWorkspace Flag to indicate if the directory path is empty and a default
		 * PEACE workspace file is to be created.
		 */
		public BackgroundLoader(final File wsPath, final boolean createDefWorkspace) {
			this.wsPath = wsPath;
			this.createDefWorkspace = createDefWorkspace;
		}
		@Override
		protected String doInBackground() throws Exception {
			return useWorkspace(this, createDefWorkspace, wsPath);
		}
		@Override
		protected void process(List<String> entries) {
			WorkspaceChooser.this.updateProgress(entries);
		}
		@Override
		protected void done() {
			String result = null;
			try {
				// Get overall status from background thread.
				result = get();
			} catch (Exception exp) {
				// If we have trouble getting status, then assume
				// the status failed.
				result = Utilities.toString(exp);
			}
			// Let helper method in main class to do some more 
			// work to keep this inner class small and streamlined.
			WorkspaceChooser.this.done(result, wsPath, createDefWorkspace);
		}
		/**
		 * Convenience method to report initial step count or progress.
		 * 
		 * This method provides a convenient public method to report
		 * progress via the protected {@link #publish(String...)} method
		 * in this class.
		 * 
		 * @param entry The progress to be reported via the main Swing GUI
		 * thread. This parameter cannot be null.
		 */
		public void updateProgress(String entry) {
			this.publish(entry);
		}
	}
	
	/** The constructor.
	 * 
	 * The constructor essentially sets up the default properties of
	 * this frame and logs the top image for display. The core setup
	 * of various tabs is done when the createTabs() method is called.
	 */
	public WorkspaceChooser(JFrame owner, PEACE peace) {
		super(owner, "PEACE: Choose Workspace", true);
		this.peace = peace; // save for call back later on.
		setLayout(new BorderLayout(0, 0));
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		this.getContentPane().setBackground(Color.white);
		setResizable(false);
		// Setup icon in title bar.
		this.setIconImage(Utilities.getIcon("images/16x16/PEACE.png").getImage());
		// Create a CustomPanel with the peace logo in it. This panel
		// will contain other components.
		cp = new CustomPanel(new BorderLayout(0, 0));
		cp.setBackground(Color.white);
		cp.setImage("images/peace_blue_header.png");
		// Use the image width to setup empty border for logo spacer
		JLabel imgHolder = new JLabel(Utilities.getIcon("images/peace_blue_header.png"));
		// Use the logo's width as the preferred width for this dialog.
		Dimension size = imgHolder.getPreferredSize();
		size.height = 400;
		setPreferredSize(size);
		setSize(size);
		// Add logo and tabs to the main dialog
		cp.add(Box.createVerticalStrut(100), BorderLayout.NORTH);
		// Check and setup the firstLaunch flag suitably.
		File defDir = new File(Utilities.getDefaultDirectory());
		firstLaunch = !(defDir.isDirectory() && defDir.exists() &&
				defDir.canRead() && defDir.canWrite());
		// Add the custom panel to this dialog
		add(cp, BorderLayout.CENTER);
	}

	/**
	 * Helper method to setup the tabs to be displayed in this frame.
	 * This method creates the various tabs and configures the 
	 * tabs for "First time Launch" as well.
	 * 
	 * @return This method returns true if the data in the tabs and
	 * license information was populated successfully.  On errors,
	 * this method returns false (after displaying suitable error
	 * messages).
	 */
	public boolean createTabs() {
		// Create the tabs in the work space chooser to contain
		// the actual work space chooser dialog and the license info.
		// Maybe we can add a tab for help to include the HTML help file?
		tabs = new JTabbedPane();
		JComponent chooserDialog = createChooserDialog(firstLaunch);
		// Create label to draw the title with a bit of border to
		// make the tab a bit larger than usual (for fancy)
		JLabel tabTitle = new JLabel("Workspace", Utilities.getIcon("images/16x16/DataSet.png"),
				JLabel.LEFT);
		tabTitle.setBorder(new EmptyBorder(4, 0, 1, 0));
		tabs.add(chooserDialog);
		tabs.setTabComponentAt(0, tabTitle);
		try {
			tabs.addTab("License", 
					Utilities.getIcon("images/16x16/GPL.png"),
					createLicenceTab(firstLaunch));
		} catch (Exception e) {
			JPanel msg = Utilities.collapsedMessage("<html>" +
					"Unable to load license data for PEACE.<br>" +
					"This indicates a problem with your JAR file.<br>" +
					"Try to download a fresh install " +
					"of PEACE and try again.</html>", 
					Utilities.toString(e));
			// Error occurred when loading the license file. This is no
			// good. Show an error message and bail out.
			JOptionPane.showMessageDialog(this, msg, "Initialization error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		// Add tabs to the custom panel
		cp.add(tabs, BorderLayout.CENTER);

		// Check and setup tabs for first launch
		if (firstLaunch) {
			// during first launch only the license tab is enabled and
			// in focus. During other launches the work space tab is the
			// the one where the user starts.
			tabs.setEnabledAt(0, false);
			tabs.setSelectedIndex(1);
			addWindowFocusListener(new WindowAdapter() {
				public void windowGainedFocus(WindowEvent evt) {
					if (showFirstTimeMsg) {
						showFirstTimeMsg = false;
						JOptionPane.showMessageDialog(tabs, FirstTimeMsg, 
								"Welcome!", JOptionPane.INFORMATION_MESSAGE);
					}					
				}
			});
		}
		// Everything went well.
		return true;
	}

	/**
	 * Obtain the work space that was successfully selected by
	 * the user. If a work space was not successfully selected
	 * then this method returns null.
	 * 
	 * @return The workspace selected by the user. If a valid workspace
	 * has not yet been selected, then this method returns null.
	 */
	public String getWorkspace() {
		return workspace;
	}
	
	/**
	 * Helper method to streamline code in {@link #createChooserDialog(boolean)} method.
	 * 
	 * This method is a helper method that is used to create a compound panel that
	 * can display one of the following two pieces of information:
	 * <ul>
	 * <li>A panel containing "OK" and "Cancel" buttons which is initially displayed</li>
	 * <li>A panel containing a progress bar to indicate progress in loading files. 
	 * This panel is displayed only after the user clicks on the "OK" button.</li>
	 * </ul>
	 *  
	 * This method is invoked only once from the createChooserDialog button.
	 * 
	 * @return A compound panel with two sub-panels to ease switching their display.
	 */
	private JPanel createButtonProgBarPanel() {
		// Create the ok and cancel buttons at the bottom.
		JButton okButton = Utilities.createButton(null, "    OK    ", "ok", 
				this, "Use sepecified workspace to launch main window", true);
		// Cancel button.
		JButton cancelButton = Utilities.createButton(null, " Cancel ", "cancel", 
				this, "Exit without performing any further operations", true);
		// Put ok and cancel button into a horizontal layout with		
		Box buttonPanel = Box.createHorizontalBox();
		buttonPanel.add(Box.createHorizontalGlue()); // right align
		buttonPanel.add(okButton);
		buttonPanel.add(Box.createHorizontalStrut(10));
		buttonPanel.add(cancelButton);
		buttonPanel.add(Box.createHorizontalStrut(10));
		
		// Next create a panel with an empty title label and a progress bar.
		// The range of the progress bar and the titles are set later on.
		progBar     = new JProgressBar(SwingConstants.HORIZONTAL);
		progBar.setIndeterminate(true);
		progBar.setStringPainted(true);
		progBar.setString("Please wait...");
		progBarInfo = new JLabel("Loading ...", null, JLabel.HORIZONTAL);
		// Place progress bar below label.
		Box progPanel = Box.createVerticalBox();
		progPanel.add(Box.createVerticalGlue()); // center align along with trailing glue
		progPanel.add(progBarInfo);
		progBarInfo.setAlignmentX(0);
		progPanel.add(Box.createVerticalStrut(5));
		progPanel.add(progBar);
		progBar.setAlignmentX(0);
		progPanel.add(Box.createVerticalStrut(5));
		progPanel.add(Box.createVerticalGlue()); // needed for center vertical alignment
		
		// Finally pack the button and progress panels into a main panel.
		btnProgPanel = new JPanel(new CardLayout());
		btnProgPanel.add(buttonPanel, "buttonPanel");
		btnProgPanel.add(progPanel,   "progressPanel");
		return btnProgPanel;
	}
	
	/**
	 * This is a helper method that creates the work space chooser 
	 * dialog that prompts the user to choose a work space for this
	 * session of the GUI. This method places all the components into
	 * a suitable panel and returns the panel back to the caller.
	 * Most of the code is dedicated to just laying out the various 
	 * components correctly in the dialog.
	 * 
	 * @param firstTime This flag indicates if PEACE is being run
	 * for the first time ever (or not).
	 * 
	 * @return The panel that contains the work space chooser dialog
	 * components.
	 */
	private JComponent createChooserDialog(boolean firstTime) {
		// Create two labels to display some basic prompt for inputs.
		JLabel info = 
			new JLabel("<html>PEACE stores your workspace information "    +
					"and other data files<br/>in a directory. This directory " +
					"is called a <i>workspace</i>.<br/>Choose a workspace " +
					"folder to use for this GUI session.</html>");
		info.setIcon(UIManager.getIcon("OptionPane.informationIcon"));
		Utilities.adjustFont(info, 1, 10, -1);
		info.setBorder(new EmptyBorder(20, 20, 30, 20));
		
		// Create browse button.
		JButton browse = new JButton("Browse...");
		browse.setActionCommand("browse");
		browse.addActionListener(this);
		browse.setAlignmentY(0); // align to top

		// Create the combo-box with recent work spaces and the browse
		// button to choose a work space in a sub-panel.
		wsList = new JComboBox();
		wsList.setEditable(true);
		Dimension maxSize = wsList.getPreferredSize();
		// Need some annoying logic to get the layout correct in 
		// windows and Linux. Don't know what Mac does with it
		maxSize.width = getPreferredSize().width - 
			browse.getPreferredSize().width - 90;
		wsList.setPreferredSize(maxSize);
		wsList.setMaximumSize(maxSize);
		wsList.setAlignmentY(0); // align to top
		// Load the last used work spaces from the .workspace_list
		// file.
		loadWorkspaceList(firstLaunch, wsList);

		// Modify Sub-panel to put fileList to left and browse button to right.
		Box horizPanel = Box.createHorizontalBox();
		horizPanel.add(Box.createHorizontalStrut(10));
		horizPanel.add(wsList);
		horizPanel.add(Box.createHorizontalStrut(10));
		horizPanel.add(browse);
		horizPanel.add(Box.createHorizontalStrut(10));
		// Put ok and cancel button into a horizontal layout with
		// some version information to the left.
		String version = "<html>" + Version.GUI_VERSION + "</html>";
		version = version.replaceAll("\n", "<br>"); // Convert newline to HTML <br>
		JLabel versionInfo = new JLabel(version, 
				Utilities.getIcon("images/24x24/PEACE.png"), JLabel.LEFT);
		JPanel buttonProgPanel = createButtonProgBarPanel();
		// Put ok and cancel button into a horizontal layout with
		// some version information to the left.
		Box buttonPanel = Box.createHorizontalBox();
		buttonPanel.add(versionInfo);
		buttonPanel.add(buttonProgPanel);
		
		// Put the above components into a suitable
		// panel with vertical layout.
		Box vertPanel = Box.createVerticalBox();
		info.setAlignmentX(0.0f);
		vertPanel.add(info);
		JLabel prompt = new JLabel("  Select Workspace:");
		prompt.setAlignmentX(0.0f);
		vertPanel.add(prompt);
		horizPanel.setAlignmentX(0.0f);
		vertPanel.add(horizPanel);
		buttonPanel.setAlignmentX(0.0f);
		vertPanel.add(Box.createVerticalGlue());
		vertPanel.add(buttonPanel);
		// Place the vertPanel in a panel.
		JPanel panel = new JPanel(new BorderLayout());
		panel.add(vertPanel, BorderLayout.CENTER);
		panel.setBorder(new CompoundBorder(new EtchedBorder(EtchedBorder.LOWERED),
				new EmptyBorder(5, 20, 5, 20)));
		return panel;
	}

	/**
	 * This is a helper method that is used to create the license
	 * tab in the work space chooser. This method creates an 
	 * accept button at the bottom of the license tab if the
	 * firstLaunch instance variable is "true"
	 * 
	 * @param firstTime If this parameter is true, then this method
	 * creates an accept button at the bottom of the license panel. In
	 * addition, it also adds an introductory message to the license
	 * information. 
	 * 
	 * @return A panel containing the license/copyright information
	 * to be displayed to the user.
	 */
	private JComponent createLicenceTab(boolean firstTime) throws Exception {
		JTextArea space = new JTextArea();
		space.setEditable(false);
		// Load a set the license in this space.
		space.setText(Utilities.readSmallTextFile("installFiles/License.txt"));
		space.setCaretPosition(0);
		// Wrap the license space in a scroll pane to let users
		// scroll through the license information.
		JScrollPane jsp = new JScrollPane(space);
		jsp.setBorder(new EtchedBorder(EtchedBorder.RAISED));

		// If this is the first launch then display an accept button
		// at the bottom.
		JPanel container = new JPanel(new BorderLayout());
		if (firstTime) {
			space.insert(FirstTimeMsg + Separator, 0);
			JButton accept = new JButton("  Accept  ");
			accept.setActionCommand("accept");
			accept.addActionListener(this);
			JPanel buttonPanel = new JPanel();
			buttonPanel.add(accept);
			container.add(buttonPanel, BorderLayout.SOUTH);
		}
		// Add scroll pane that has license text in it to container.
		container.add(jsp, BorderLayout.CENTER);
		// Return the container
		return container;
	}

	/**
	 * This is a helper method that is used to load the list of previously
	 * used work space paths from the global property named 
	 * "WorkspaceChooser.WorkspaceList". The properties are already
	 * loaded from the properties file in the default working directory. 
	 * If no previous work spaces were found, then this method adds a 
	 * suggested default work space for the user.
	 * 
	 * @param firstTime If this flag is true then it indicates that the
	 * PEACE program is being launched for the first time.
	 * 
	 * @param list The combo box to be populated with the list of work 
	 * space paths.
	 */
	private void loadWorkspaceList(boolean firstTime, JComboBox list) {
		// If this is not the first time running, load data from the
		// "WorkspaceChooser.WorkspaceList" property. Otherwise we use
		// a default.
		final String defWSPath = Utilities.getDefaultDirectory() + 
			File.separator + "workspace1";
		final String workspaceList = PEACEProperties.get().getProperty(WSLIST_PROPERTY_NAME, defWSPath);
		Scanner wsFile     = new Scanner(workspaceList);
		// Read line-by-line from file and add it to combo box
		while (wsFile.hasNext()) {
			String wsPath = wsFile.nextLine().trim();
			if (wsPath.length() < 1) {
				// Ignore empty lines.
				continue;
			}
			// Add items.
			wsList.addItem(wsPath);
		}
		// We will never have an empty list.
		assert(wsList.getItemCount() > 0 );
		// Select the first item in the list as the default work space
		wsList.setSelectedIndex(0);
	}

	/**
	 * This is a helper method that is used to save the current list of
	 * work space paths into the global properties.
	 * 
	 * This method updates the list of files in the global properties and
	 * saves the properties to disk. 
	 * 
	 * @param list The combo box whose data is to be saved as the list of
	 * workspaces.
	 */
	private void saveWorkspaceList(JComboBox list) {
		StringBuilder sb = new StringBuilder(1024);
		for(int idx = 0; (idx < list.getItemCount()); idx++) {
			sb.append(list.getItemAt(idx));
			sb.append('\n');
		}
		final String workspaceList = sb.toString();
		// Update properties and save the list.
		PEACEProperties.get().setProperty(WSLIST_PROPERTY_NAME, workspaceList);
		PEACEProperties.get().save(list, true);
	}

	/**
	 * This method handles the actions associated with various controls
	 * in the dialog. Specifically it handles the actions associated
	 * with the following components: the "Accept" button (shown only
	 * the first time the GUI is launched to accept the license), the
	 * "Browse..." button, the "OK" button and the "Cancel" button.
	 * 
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		String command = event.getActionCommand();
		if ("accept".equals(command)) {
			// Use helper method to do the various actions.
			accept(event);
		}

		if ("ok".equals(command)) {
			// The user must have chosen a work space. Check if it exists
			// or create it otherwise. Close this dialog.
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			okAction();
			setCursor(Cursor.getDefaultCursor());
		}
		
		// The browse action is sufficiently small to be inlined here.
		if ("browse".equals(command)) {
			// Create and display a file chooser for user to browse
			JFileChooser jfc = new JFileChooser();
			jfc.setDialogTitle("Choose workspace directory");
			jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			if (jfc.showDialog(this, "Use As Workspace") == JFileChooser.APPROVE_OPTION) {
				// Copy the chosen directory to the work space combo box.
				String wsPath = jfc.getSelectedFile().getAbsolutePath();
				// Set the selected item in this path.
				wsList.setSelectedItem(wsPath);
			}
		}
		
		// The cancel action is sufficiently small to be in-lined here.
		if ("cancel".equals(command)) {
			this.setVisible(false);
		}
	}

	/**
	 * This is a refactored (to keep code clutter to a minimum in the
	 * actionPerformed()  method) method that is called only from
	 * the actionPerformed() method whenever the user clicks on the
	 * "Accept" button in the license tab. This button is displayed
	 * only when the GUI is launched for the first time. Accordingly,
	 * this method performs the following action(s):
	 * 
	 * <ol>
	 * 
	 * <li>Creates a default working directory. If errors occur then
	 * this method displays an error and returns immediately.</li>
	 * 
	 * <li>It hides the "Accept" button (as the user has already
	 * accepted the license there is no need to display it).</li>
	 * 
	 * <li>It activates the "Workspace" tab so user can now choose
	 * a workspace.</li>
	 * 
	 * <li>It informs the user that the default working directory has
	 * been created.</li>
	 * 
	 * <li>Switches tabs to the "Workspace" tab so user can choose the
	 * workspace.</li>
	 * 
	 * </ol>
	 * 
	 * @param event The action event that triggered this action.
	 */
	private void accept(ActionEvent event) {
		// The user has accepted the license. First create default
		// directory. 
		String path = Utilities.getDefaultDirectory();
		File dir = new File(path);
		if (!dir.mkdir()) {
			// Unable to create the working folder.
			JOptionPane.showMessageDialog(this, 
					"Unable to create the default working directory " +
					"for PEACE.\nEnsure the following directory does "  +
					"not exist but can be created:\n" +
					"Directory: " + path,
					"Error creating default directory", 
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		// OK, the path was successfully created. 
		// First, hide the accept button
		JButton button = (JButton) event.getSource();
		button.getParent().setVisible(false);
		//  Next, enable the work space tab.
		tabs.setEnabledAt(0, true);
		// Show message that default work directory was created.
		JOptionPane.showMessageDialog(this, 
				"The default working directory for PEACE was created.\n" +
				"(Directory: " + path + ")",
				"Default directory created", 
				JOptionPane.INFORMATION_MESSAGE);
		// Now switch tabs over to the work space tab.
		tabs.setSelectedIndex(0);
	}
	
	/**
	 * This is a refactored (to keep code clutter to a minimum in the
	 * actionPerformed()  method) method that is called only from
	 * the actionPerformed() method whenever the user clicks on the
	 * "OK" button in the workspace tab. This method performs the
	 * following operations.
	 * 
	 * <ol>
	 * 
	 * <li>Checks and creates the selected work space directory. If the
	 * directory already exists checks are made to ensure that the
	 * directory is actually a directory and is readable.</li>
	 * 
	 * <li>If errors occur during directory creation or verification
	 * then a suitable error message is displayed and the method
	 * exists immediately.</li>
	 * 
	 * <li>The list of recently used workspaces is updated and the
	 * list is saved.</li>
	 * 
	 * <li>The method sets up the workspace instance variable and
	 * disposes the dialog to return control back to the main method.</li>
	 * 
	 * </ol>
	 * 	 
	 */
	private void okAction() {
		// Flag to indicate if user selected path is new or not 
		boolean createDefWorkSpace = false;
		// The user must have chosen a work space. Check if it exists
		// or create it otherwise.
		final File wsPath = new File((String) wsList.getSelectedItem());
		if (!wsPath.exists()) {
			// Try and create the working directory
			if (!wsPath.mkdir()) {
				JOptionPane.showMessageDialog(this,
						"Unable to create the specified workspace directory.\n" +
						"(Directory: " + wsPath + ")\n" + 
						"Choose a different workspace directory.",
						"Unable to create Workspace", JOptionPane.ERROR_MESSAGE);
				// Use early exit to avoid deeply nested if-else
				return;
			}
			// Flag to indicate path is new to create default workspace below.
			createDefWorkSpace = true;
		} else {
			// The path exists. Ensure it is usable.
			if (!wsPath.isDirectory() || !wsPath.canRead()) {
				JOptionPane.showMessageDialog(this,
						"The workspace directory is not usable.\n" +
						"(Directory: " + wsPath + ")\n" +
						"Choose a different workspace directory.",
						"Invalid Workspace", JOptionPane.ERROR_MESSAGE);
				// Use early exit to avoid deeply nested if-else
				return;
			}
			// Check to see if default work space file exists.
			String filename = Workspace.getWorkspaceFile(wsPath.getAbsolutePath());
			File   wsDataFile = new File(filename);
			
			if (!wsDataFile.exists()) {
				int option = JOptionPane.showConfirmDialog(this, 
						"The directory already exists but does not have " +
						"valid workspace metadata.\n" +
						"(directory: " + wsPath.getAbsolutePath() + ")\n" +
						"PEACE will need to create a default workspace " +
						"metadata.\n" +
						"Do you want to create a default workspace?",
						"Checking workspace", JOptionPane.YES_NO_OPTION);
				if (option == JOptionPane.NO_OPTION) {
					// Let user re-choose a folder in this case.
					return;
				}
				// Make a default work space in the directory.
				createDefWorkSpace = true;
			} 
		}
		// The loading of the workspace is not a short process.
		// This is done in a background thread to ensure GUI remains
		// response and we are able to show progress information.
		CardLayout cl = (CardLayout) btnProgPanel.getLayout();
		cl.show(btnProgPanel, "progressPanel");
		// call useWorkspace() method from a separate background thread.
		BackgroundLoader bl = new BackgroundLoader(wsPath, createDefWorkSpace);
		bl.execute();
	}
	
	private void updateProgress(List<String> entries) {
		if ((entries == null) || (entries.isEmpty())) {
			// No entries to process. Nothing further to be done.
			return;
		}
		// Check and handle initialization situation first.
		int startEntry = 0;
		if (progBar.isIndeterminate()) {
			// This is the first update entry from our background thread.
			// This value is an integer that indicates the number of
			// operations to be performed. This is used to set the 
			// range for the progress bar that displays progress.
			int numOperations = Integer.parseInt(entries.get(0));
			progBar.setIndeterminate(false);
			progBar.setMaximum(numOperations);
			progBar.setMinimum(0);
			progBar.setValue(0);
			startEntry++; // Skip over first entry in for-loop below.
		}
		// Check and process any additional entries in the list by
		// setting current status to the last entry
		if (startEntry < entries.size()) {
			final int lastEntry = entries.size() - 1;
			progBarInfo.setText(entries.get(lastEntry));
			progBar.setValue(progBar.getValue() + lastEntry + 1 - startEntry);
		} 
	}
	
	private void done(final String result, 
			final File wsPath, final boolean createDefWorkSpace) {
		// First get overall result from the background worker thread.
		boolean success = Boolean.TRUE.toString().equals(result);
		// If data was successfully loaded then move forward with
		// launching the main frame
		if (success) {
			// Dispose this dialog returning control to main.
			setVisible(false);
			dispose();
		} else {
			// Error loading workspace information.
			// Re-enable OK/Cancel button and bail out
			CardLayout cl = (CardLayout) btnProgPanel.getLayout();
			cl.show(btnProgPanel, "buttonPanel");
			progBar.setIndeterminate(true);
			// Error creating/using workspace. Report error.
			final String subMsg = (createDefWorkSpace ?
				"Error creating initial work space in directory" :
				"The chosen directory has invalid PEACE workspace metadata");
			JPanel msg = Utilities.collapsedMessage("<html>" + subMsg +
					"<br/>(Directory: " + wsPath.getAbsolutePath() + ")<br/>" +
					" or there was an error loading DECAGON assembler description.<br/><br/>" +					
					"PEACE cannot use the above directory and operate correctly.<br/>" +
					"Pleace choose a different workspace folder.</html>", 
					result);
			JOptionPane.showMessageDialog(this, msg,
					"Cannot use workspace", JOptionPane.ERROR_MESSAGE);
		}
	}
	
	private String useWorkspace(BackgroundLoader loader,
			final boolean createDefWorkSpace, final File wsPath) {
		// Set standard set of DECAGON Assembler Description XML (DADX) files 
		// to be loaded.
		final List<String> dadxList = DADXHelper.getDADXFilesToLoad();
		// Let the GUI know this background method has done some work.
		loader.updateProgress("" + (dadxList.size() + 3));
		// Now load each DADX file one at a time from the list and report progress
		for(int i = 0; (i < dadxList.size()); i++) {
			try {
				// Create helper and get it to load DADX file and register
				// it with DADXSummary class.
				UserLog.log(LogLevel.INFO, "WorkspaceChooser", "Loading DADX file " + dadxList.get(i));
				// Extract just the file name from the full path to the file for progress bar status
				File tmpFile = new File(dadxList.get(i));
				loader.updateProgress("Loading: " + tmpFile.getName());
				DADXHelper helper = new DADXHelper();
				helper.unmarshal(dadxList.get(i), true, false);
			} catch (Exception exp) {
				// Error loading DADX file. Bail out from background processing
				// with an exception message.
				return "Error loading DADX file: " + dadxList.get(i) + "\n" + 
					Utilities.toString(exp);
			}
		}
		// Now load the actual PEACE workspace file.
		try {
			// Update progress
			loader.updateProgress("Loading workspace...");
			// The value for the createDefWorkSpace flag is set in
			// the okAction() method.
			if (createDefWorkSpace) {
				// Create default workspace in the given path.
				Workspace.createDefault(wsPath.getAbsolutePath());
			} else {
				// Load the existing workspace definition from the given
				// directory.
				Workspace.useWorkspace(wsPath.getAbsolutePath());
			}
		} catch (Exception exp) {
			// Error creating/using workspace. Report error.
			return Utilities.toString(exp);
		}
		// The work space is good and is usable. First save the
		// recently used work space information.
		if (wsList.getSelectedIndex() != -1) {
			// Remove existing item to avoid duplicates below.
			wsList.removeItemAt(wsList.getSelectedIndex());
		}
		// Workspace is successfully read for further use.
		// Make the selected item the first item in the list
		// for the next launch.
		wsList.insertItemAt(wsPath.getAbsolutePath(), 0);
		wsList.setSelectedIndex(0);
		// Now set the workspace instance variable.
		workspace = wsPath.getAbsolutePath();
		// Let the main PEACE class know it can now launch the main frame.
		loader.updateProgress("Launching Mainframe...");
		peace.launchMainFrame(workspace, firstLaunch);
		// save the list for future use at a later time.
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				saveWorkspaceList(wsList);
			}
		});
		// Everything went well.
		loader.updateProgress("Done.");
		return Boolean.TRUE.toString();
	}
	
	/**
	 * Reference to the main PEACE class. The reference is used to
	 * create the main frame after the user has selected a valid
	 * work space.
	 */
	private final PEACE peace;
	
	/**
	 * A custom panel that actually has the PEACE logo as the 
	 * background image. This panel is used to hold all the 
	 * components just to make the GUI look nice.
	 */
	private CustomPanel cp;
	
	/**
	 * The generated serial version ID to keep the Java compiler
	 * happy about serializing this class. 
	 */
	private static final long serialVersionUID = 5779507970297217323L;

	/**
	 * A static piece of text that is prepended to the license text
	 * when PEACE is launched for the first time to provide the user
	 * with some additional information.
	 */
	private static final String FirstTimeMsg = 
		"It appears that this is the first time you are running PEACE GUI.\n" +
		"Thank you for trying out PEACE. Prior to using PEACE you need to:\n" +
		"    1. Read and accept the license below.\n" +
		"    2. Allow PEACE GUI to create a default working folder.";
	
	/**
	 * A simple separator to make the license text look pretty.
	 */
	private static final String Separator = 
		"\n\n---------------------------[ License ]----------------------------\n\n";
	
	/**
	 * The name of the property that is used to save the workspace list.
	 * 
	 * This string indicate the name of the property (in the global list of
	 * properties) that is used to save the newline delimited list of
	 * workspaces that have been used by the user. 
	 */
	private static final String WSLIST_PROPERTY_NAME = "WorkspaceChooser.WorkspaceList";
}
