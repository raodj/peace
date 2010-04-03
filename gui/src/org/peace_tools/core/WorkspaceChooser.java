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
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
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
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.UIManager;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;

import org.peace_tools.generic.CustomPanel;
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
 * <li>When a user selects a valid work space this method sets up
 * the work space instance variable and the dialog winds up.</li>
 * 
 * <li>The main method uses the work space directory for further
 * use to create the main frame.</li>
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
		// Create the ok and cancel buttons at the bottom.
		JButton okButton = new JButton("    OK    ");
		okButton.addActionListener(this);
		okButton.setActionCommand("ok");
		// Cancel button.
		JButton cancelButton = new JButton(" Cancel ");
		cancelButton.addActionListener(this);
		cancelButton.setActionCommand("cancel");
		// Put ok and cancel button into a horizontal layout with
		// some version information to the left.
		String version = "<html>" + Version.GUI_VERSION + "</html>";
		version = version.replaceAll("\n", "<br>"); // Convert newline to HTML <br>
		JLabel versionInfo = new JLabel(version, 
				Utilities.getIcon("images/24x24/PEACE.png"), JLabel.LEFT);
		Box buttonPanel = Box.createHorizontalBox();
		buttonPanel.add(versionInfo);
		buttonPanel.add(Box.createHorizontalGlue()); // right align
		buttonPanel.add(okButton);
		buttonPanel.add(Box.createHorizontalStrut(10));
		buttonPanel.add(cancelButton);
		buttonPanel.add(Box.createHorizontalStrut(10));
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
	 * used work space paths from the default ".workspace_list" file in the
	 * default working directory. If no previous work spaces were found, then
	 * this method adds a suggested default work space for the user.
	 * 
	 * @param firstTime If this flag is true then it indicates that the
	 * PEACE program is being launched for the first time.
	 * 
	 * @param list The combo box to be populated with the list of work 
	 * space paths.
	 */
	private void loadWorkspaceList(boolean firstTime, JComboBox list) {
		// If this is not the first time running, load data from the
		// ".workspace_list" file.
		try {
			String wsFileName  = Utilities.getDefaultDirectory() + "/.workspace_list";
			Scanner wsFile     = new Scanner(new FileReader(wsFileName));
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
			wsFile.close();
		} catch (Exception e) {
			// Do nothing if the work space file could not be read.
		}
		// Add a default work space entry to the list.
		if (wsList.getItemCount() == 0) {
			String wsPath = Utilities.getDefaultDirectory() + 
				File.separator + "workspace1";
			wsList.addItem(wsPath);
		}
		// Select the first item in the list as the default work space
		wsList.setSelectedIndex(0);
	}

	/**
	 * This is a helper method that is used to save the current list of
	 * work space paths to the work space list file. 
	 * 
	 * @param list The combo box whose data is to be saved.
	 */
	private void saveWorkspaceList(JComboBox list) {
		try {
			String wsFileName  = Utilities.getDefaultDirectory() + "/.workspace_list";
			PrintWriter printer = new PrintWriter(new FileWriter(wsFileName));
			for(int idx = 0; (idx < list.getItemCount()); idx++) {
				printer.println(list.getItemAt(idx));
			}
			printer.close();
		} catch (Exception e) {
			// Ignore errors on saving work space list
		}
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
		boolean makeDefWS = false;
		// The user must have chosen a work space. Check if it exists
		// or create it otherwise.
		File wsPath = new File((String) wsList.getSelectedItem());
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
			makeDefWS = true;
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
				makeDefWS = true;
			} else {
				// A work space file exists.
				try {
					// Load and check to ensure it is valid.
					Workspace.useWorkspace(wsPath.getAbsolutePath());
				} catch (Exception e) {
					// Workspace metadata is invalid. User must chose a
					// different directory.
					JPanel msg = Utilities.collapsedMessage("<html>" +
							"The chosen directory has invalid PEACE workspace " +
							"metadata:<br>" +
							"Directory: " + wsPath.getAbsolutePath() + "<br><br>" + 
							"PEACE cannot use the above directory and operate correctly.<br>" +
							"Pleace choose a different workspace folder.</html>", 
							Utilities.toString(e));
					JOptionPane.showMessageDialog(this, msg,
						"Invalid workspace directory", JOptionPane.ERROR_MESSAGE);
					// Bail out.
					return;
				}
			}
		}
		// If a new work space needs to be created do it now.
		if (makeDefWS) {
			try {
				Workspace.createDefault(wsPath.getAbsolutePath());
			} catch (Exception e) {
				// O!o! error creating default work space. Report error
				// and let user chose another directory.
				JPanel msg = Utilities.collapsedMessage("<html>" +
						"Error creating initial work space in directory.<br>" +
						"Directory: " + wsPath.getAbsolutePath() + "<br><br>" +
						"PEACE cannot use the above directory and operate correctly.<br>" +
						"Pleace choose a different workspace folder.</html>", 
						Utilities.toString(e));
				JOptionPane.showMessageDialog(this, msg,
						"Cannot create workspace", JOptionPane.ERROR_MESSAGE);
				// bail out
				return;
			}
		}
		
		// The work space is good and is usable. First save the
		// recently used work space information.
		if (wsList.getSelectedIndex() != -1) {
			// Remove existing item to avoid duplicates below.
			wsList.removeItemAt(wsList.getSelectedIndex());
		}
		// Make the selected item the first item in the list
		// for the next launch.
		wsList.insertItemAt(wsPath.getAbsolutePath(), 0);
		// save the list for future use
		saveWorkspaceList(wsList);
		// Now set the workspace instance variable.
		workspace = wsPath.getAbsolutePath();
		// Let the main PEACE class know it can now launch the main frame.
		peace.launchMainFrame(workspace, firstLaunch);
		// Dispose this dialog returning control to main.
		setVisible(false);
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
	private static String FirstTimeMsg = 
		"It appears that this is the first time you are running PEACE GUI.\n" +
		"Thank you for trying out PEACE. Prior to using PEACE you need to:\n" +
		"    1. Read and accept the license below.\n" +
		"    2. Allow PEACE GUI to create a default working folder.";
	
	/**
	 * A simple separator to make the license text look pretty.
	 */
	private static String Separator = 
		"\n\n---------------------------[ License ]----------------------------\n\n";
}
