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

package org.peace_tools.views;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.Date;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextPane;
import javax.swing.plaf.ComponentUI;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultStyledDocument;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

/**
 * A common base class shared by JobDetailsView and ServerJobsView.
 *
 * This is a common base class that is shared between two very similar
 * DetailView classes, namely JobDetailsView and ServerDetailsView.
 * The base class houses several of the commonly used features of these
 * two detail viewer classes, such as:
 * 
 * <ul>
 * 
 * <li>createTopPanel: Helper method to layout the top panel that 
 * consists of job ID/user ID following by the server name from 
 * where the details are being displayed. The second row has a
 * refresh button, progress bar, following by time stamp when the
 * information was last refreshed.</li>
 * 
 * <li>createDocument: A helper method to create a DefaultStyledDocument 
 * that can be used to display information with specific styles</li>
 * 
 *  <li>
 * 
 * </ul>
 */
abstract class DetailView extends JPanel 
implements ActionListener, Runnable {
	/**
	 * The constructor.
	 * 
	 * The constructor initializes the various instance variables and
	 * creates the different tabs and panels.
	 *  
	 * @param idTitle The label or title to be displayed for the 
	 * ID value. This value is "Job ID" or "User ID". 
	 *
	 * @param idValue The actual value to be displayed for the label.
	 * This is typically a job ID or a user ID depending on the type
	 * of the detail view.
	 * 
	 * @param server The server from where the details are to be 
	 * displayed.
	 */
	public DetailView(String idTitle, String idValue, Server server) {
		// Save reference to the server for later use
		this.server = server;
		// Setup the panels
		this.setOpaque(false);
		setLayout(new BorderLayout(0, 0));
		add(createTopPanel(idTitle, idValue), BorderLayout.NORTH);
	}

	/**
	 * Method to intercept clicks on "Refresh" button.
	 * 
	 * This method is invoked whenever the user clicks on the
	 * "Refresh" button (or from the constructor). This method
	 * launches a background thread that performs the task of
	 * obtaining and populating job information from the remote
	 * server.
	 * 
	 * <p><b>Note:</b>  The derived class must implement the 
	 * Runnable.run() method.</p>
	 * 
	 * @param arg0 The action event associated with a button
	 * click. This event is not really used. So you can call this
	 * method with null as parameter.
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		// First disable the refresh button.
		refreshButton.setEnabled(false);
		// Start background thread.
		Thread t = new Thread(this);
		t.start();
	}
		
	/**
	 * Helper method to update information in a document.
	 * 
	 * This is a helper method that is used to set/update/reset 
	 * information in documents. This method is specifically used
	 * to update information in tabs such as: outputs, error, and 
	 * scripts tabs.
	 * 
	 * @param doc The document whose information is to be cleared
	 * and then updated.
	 * 
	 * @param title The title for the document.
	 * 
	 * @param text The actual content to be displayed in the document.
	 * 
	 * @param textStyle A predefined style to be set for the text.
	 * 
	 * @throws BadLocationException This method simply exposes the
	 * exception generated by the document.
	 */
	protected void updateDocument(DefaultStyledDocument doc, String title, 
			String text, String textStyle) throws BadLocationException {
		// First clear out all the information
		doc.replace(0, doc.getLength(), "", null);
		doc.insertString(doc.getLength(), title + "\n", 
				doc.getStyle("title"));
		doc.insertString(doc.getLength(), text, doc.getStyle(textStyle)); 
	}
	
	/**
	 * Helper method to run the job runner on the server and return
	 * outputs.
	 * 
	 * This method is invoked from the run method to actually run the
	 * job runner on the remote machine. This method execute the job
	 * runner on the remote machine with the given parameter as the 
	 * command line parameter.
	 * 
	 * @param session The server session to be used to execute the 
	 * command. The command should include the complete path in it
	 * without the suffix. This method adds a ".sh" or ".bat"
	 * suffix to the command before executing it.
	 * 
	 * @param command The command to be run to via the specified
	 * session. 
	 * 
	 * @param option The command line option to be passed to the
	 * job runner script as the command line parameter. If a command
	 * does not have any options then set the option to an empty
	 * string ("").
	 * 
	 * @param streams The standard output and error streams that were
	 * returned by running the job runnser script with the specified
	 * option.
	 */
	protected void runCommand(ServerSession session, String command, 
			String option, String[] streams) throws Exception {
		Server.OSType os       = session.getOSType();
		final String extension = (Server.OSType.WINDOWS.equals(os)) ? "bat" : "sh"; 
		String remoteCmd = command + "." + extension + " " + option;
		// Run command and get stdin and stdout
		streams[0] = "";
		streams[1] = "";
		ProgrammerLog.log("Attempting to run command '" + remoteCmd + 
				"' on server: " + server.getName() + "\n");
		int exitCode = session.exec(remoteCmd, streams);
		ProgrammerLog.log("Exit code: '" + exitCode); 
		if (exitCode != 0) {
			throw new IOException("The remote command " + remoteCmd + "\n" +
					"terminated with exit code: " + exitCode + ".\n" +
					"The error message reported was:\n" + streams[1] + "\n.");
		}
	}
	
	/**
	 * Helper method to check and create styles.
	 * 
	 * This helper method is invoked from the createOutputDocs()
	 * method. This method creates the set of styles that are
	 * used to display information in the styled documents
	 * associated with this class. The styles are stored in a shared
	 * style context class. The styles are created only if the
	 * style context class does not exist. Therefore calling
	 * this method over and over has no side effects. Currently
	 * this method creates the following styles: stdout, stderr,
	 * scripts, title, subtitle.
	 * 
	 */
	private static synchronized void createStyles() {
		if (sc == null) {
			sc = new StyleContext();
			// Setup some standard styles we are going to use.
			Style style = sc.addStyle("stdout", null);
			style.addAttribute(StyleConstants.Foreground, Color.blue);
			style.addAttribute(StyleConstants.FontFamily, Font.MONOSPACED);
			style = sc.addStyle("stderr", null);
			style.addAttribute(StyleConstants.Foreground, Color.red);
			style.addAttribute(StyleConstants.FontFamily, Font.MONOSPACED);
			style = sc.addStyle("scripts", null);
			style.addAttribute(StyleConstants.Foreground, Color.green.darker());
			style.addAttribute(StyleConstants.FontFamily, Font.MONOSPACED);
			style = sc.addStyle("info", null);
			style.addAttribute(StyleConstants.Foreground, Color.black);
			style = sc.addStyle("title", null);
			style.addAttribute(StyleConstants.FontSize, new Integer(18));
			style = sc.addStyle("subtitle", null);
			style.addAttribute(StyleConstants.FontSize, new Integer(14));
		}
	}
	
	/**
	 * Helper method to create styled documents.
	 * 
	 * This helper method is invoked from the constructor to create
	 * the four tabs to display various information about a given
	 * job. This method creates four StyledDocuments with custom
	 * styles tagged with the names "stdout", "stderr", "scripts",
	 * "info", "title", and "subtitle". The styled documents
	 * are placed within scroll panes and added to different tabs. 
	 *
	 * @param outputs The array that must contain the various styled
	 * documents to be created by this method.
	 * 
	 * @param TabNames The titles for the various tabs to be created
	 * by this method. The number of titles must match the list of
	 * outputs -- that is, outputs.length == TabNames.length.
	 * 
	 * @return This method returns the tabbed pane that contains all
	 * the styled documents where the outputs have been added.
	 */
	protected JTabbedPane createOutputDocs(DefaultStyledDocument outputs[], 
			final String TabNames[]) {
		// Create our document with a suitable style context.
		createStyles();
		// Create the tabs
		JTabbedPane tabPane = new JTabbedPane(JTabbedPane.TOP);
		// Create the documents and add tabs
		for(int i = 0; (i < outputs.length); i++) {
			outputs[i] = new DefaultStyledDocument(sc);
			// Create the text pane to display output logs
			JTextPane textPane = new JTextPane(outputs[i]) {
				private static final long serialVersionUID = -2728430057905161506L;
				@Override
				public boolean getScrollableTracksViewportWidth() {
					Component parent = getParent();
					ComponentUI ui = getUI();
					return parent != null ? (ui.getPreferredSize(this).width <= parent
							.getSize().width) : true;
				}
			};
			// Add the text pane to a scroll pane to handle large outputs
			JScrollPane jsp = new JScrollPane(textPane);
			tabPane.add(jsp, TabNames[i]);
		}
		// Return the newly created tab
		return tabPane;
	}
	
	/**
	 * Helper method to create the controls displayed at the top.
	 * 
	 * This is a refactored method that was introduced to streamline
	 * the constructor. This helper method performs the task of 
	 * laying out the controls at the top of the job details
	 * tab. The creation and layout of all the controls is done
	 * using a GridBagLayout as it seems to work out the best
	 * in this situation.
	 * 
	 * @param idTitle The label or title to be displayed for the 
	 * ID value. This value is "Job ID" or "User ID". 
	 *
	 * @param idValue The actual value to be displayed for the label.
	 * This is typically a job ID or a user ID depending on the type
	 * of the detail view.
	 * 
	 * @return The panel to be displayed at the top of this view
	 * with all the necessary components.
	 */
	private JPanel createTopPanel(String idTitle, String idValue) {
		GridBagLayout gbLayout = new GridBagLayout();
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.fill               = GridBagConstraints.HORIZONTAL;
		// Create the top panel.
		JPanel firstRow        = new JPanel(gbLayout);
		// The list of labels we are going to create on first row
		String labels[] = {idTitle, idValue, "  Server: ", server.getName()};
		for(int i = 0; (i < labels.length); i++) {
			JLabel label = new JLabel(labels[i]);
			gbc.weightx  = 0;
			if (i % 2 != 0) {
				// Set a border to make the label look like a field
				label.setBorder(BorderFactory.createCompoundBorder(
						BorderFactory.createEtchedBorder(), 
						BorderFactory.createEmptyBorder(3, 4, 4, 4)));
				label.setBackground(Color.white);
				label.setOpaque(true);
				gbc.weightx = 1;
			}
			if (i == labels.length - 1) {
				gbc.gridwidth = GridBagConstraints.REMAINDER;
			}
			gbLayout.setConstraints(label, gbc);
			firstRow.add(label);
		}
		
		// Create components for the second row.
		refreshButton = Utilities.createButton("images/16x16/Refresh.png",
				"Refresh", "Refresh", this, "Refresh the job data", false);
		progressBar = new JProgressBar(0, 7);
		Dimension size = progressBar.getPreferredSize();
		size.width     = 100;
		progressBar.setMinimumSize(size);
		progressBar.setBorder(BorderFactory.createEmptyBorder(1, 20, 1, 20));
		progressBar.setStringPainted(true);
		progressBar.setOpaque(false);
		
		JLabel tsLabel = new JLabel("  Last update: ");
		Date logTime = new Date();
		timestamp      = new JLabel(logTime.toString());
		timestamp.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createEtchedBorder(), 
				BorderFactory.createEmptyBorder(3, 4, 4, 4)));
		timestamp.setBackground(Color.white);
		timestamp.setOpaque(true);
		
		// Layout all the components in the second row.
		gbLayout = new GridBagLayout();
		gbc      = new GridBagConstraints();
		JPanel secondRow = new JPanel(gbLayout);
		gbLayout.setConstraints(refreshButton, gbc);
		secondRow.add(refreshButton);
		gbc.weightx   = 1.0;
		gbc.fill      = GridBagConstraints.HORIZONTAL;
		gbLayout.setConstraints(progressBar, gbc);
		secondRow.add(progressBar);
		gbc.weightx   = 0;
		gbc.fill      = GridBagConstraints.NONE;
		gbLayout.setConstraints(tsLabel, gbc);
		secondRow.add(tsLabel);
		gbc.fill      = GridBagConstraints.HORIZONTAL;
		gbc.gridwidth = GridBagConstraints.REMAINDER;
		gbLayout.setConstraints(timestamp, gbc);
		secondRow.add(timestamp);
		
		// Create a top-panel with the two rows
		JPanel topPanel = new JPanel(new BorderLayout(2, 3));
		topPanel.setBorder(BorderFactory.createEmptyBorder(2, 6, 2, 6));
		topPanel.add(firstRow,  BorderLayout.NORTH);
		topPanel.add(secondRow, BorderLayout.SOUTH);
		// Make the panels opaque to make them blend in with the tab
		topPanel.setOpaque(false);
		firstRow.setOpaque(false);
		secondRow.setOpaque(false);
		
		return topPanel;
	}

	/**
	 * The refresh button that the user can click to reload the
	 * job information.
	 */
	JButton refreshButton;
	
	/**
	 * The progress bar that is placed in indeterminate mode to
	 * indicate to the user that the job data is being download.
	 */
	JProgressBar progressBar;
	
	/**
	 * This label is used to display the time stamp when the 
	 * job information was last updated.
	 */
	JLabel timestamp;
	
	/**
	 * The server from where the details are bing displayed. The 
	 * server reference is set when this class is created. This entry
	 * is used to connect to the server to obtain information.
	 */
	final protected Server server;

	/**
	 * A shared style context that contains the various styles that
	 * are used to layout information in the DefaultStyledDocument
	 * objects created and used by this class. The StyleContext is
	 * created the first time createOutputDocs()  method is called
	 */
	private static StyleContext sc = null;
	
	/**
	 * A generated serialization UID just to keep the Java compiler
	 * happy. 
	 */
	private static final long serialVersionUID = 6192661377313095929L;

}
