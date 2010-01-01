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

package org.peace_tools.core.dataset;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ProgressMonitorInputStream;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.peace_tools.data.DataStore;
import org.peace_tools.data.ESTList;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;

/**
 * This class serves as the second interactive page in a the Data Set
 * Wizard. This page permits the user to enter information about the
 * top-level EST file that forms the primary input for this data set.
 * The EST file must be in FASTA file format in order for PEACE
 * to successfully process it.  In order to ensure validity of the 
 * primary FASTA file, this wizard page verifies the file when
 * the "Next" button is clicked (and before proceeding to the 
 * next page). The verification is done on a separate thread so
 * that the GUI does not appear to hang (as loading FASTA files
 * can be a time consuming task).
 */
public class ESTInfoWizardPage extends GenericWizardPage 
implements Runnable, ActionListener {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a text area for
	 * description, a text field for install directory.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * @param dataSet The in-memory object that contains the
	 * core information about the data set being edited.
	 * @param lockPath Flag to indicate if this page must permit the
	 * user to edit the EST file path.
	 */
	public ESTInfoWizardPage(DataSetWizard wizard, DataSet dataSet,
			boolean lockPath) {
		this.wizard   = wizard;
		this.dataSet  = dataSet;
		this.lockPath = lockPath; 
		// Setup the title(s) for this page and border
		setTitle("EST Information", 
		"Set the primary input EST FASTA file for this data set");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Put the install path text file and browse button into a 
		// single horizontal panel.
		browse = Utilities.createButton(null, " Browse ", 
				"Browse", this, "Browse local file system", false);
		estFile = new JTextField();
		JPanel horizBox = new JPanel(new BorderLayout(10, 0));
		horizBox.add(estFile, BorderLayout.CENTER);
		horizBox.add(browse, BorderLayout.EAST);
		// Create a labeled component.
		JComponent dirBox =
			Utilities.createLabeledComponents("Enter EST FASTA File Name:",
					"(FASTA file must exist and must be valid)",
					0, false, horizBox);		
		// Create panels with description, install folder, and
		// polling time.
		description = new JTextArea(4, 4);
		JScrollPane jsp = new JScrollPane(description);
		jsp.setMinimumSize(description.getPreferredSize());
		JComponent descBox = 
			Utilities.createLabeledComponents("Description for EST file:",
					"(This is for your reference & can be anything)", 0, 
					false, jsp); 
		// Pack the options along with a pretty icon and label into 
		// a sub-panel.
		JComponent infoLabel = null;
		// Set up additional informational panel if the install
		// directory is editable for verification.
		if (!lockPath) {
			// Create the informational panel.
			infoLabel = createInfoPanel();
		} else {
			estFile.setEnabled(false);
			browse.setEnabled(false);
		}
		// Create a sub panel with all the components
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
				dirBox, Box.createVerticalStrut(10),
				descBox, Box.createVerticalStrut(10),
				infoLabel);
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.NORTH);
	}

	/**
	 * This is a refactored utility method to create informational
	 * labels. This method was introduced to keep the code clutter in
	 * the constructor to a minimum. 
	 * 
	 * @return The label to be used as the information panel.
	 */
	private JComponent createInfoPanel() {
		// Let the user know the remote directory will be validated when
		// they click the "Next>" button.
		JLabel info = new JLabel("<html>EST file will be verified when " + 
				"the<br>Next button is clicked</html>", 
				Utilities.getIcon("images/16x16/Information.png"),
				JLabel.LEFT);
		return info;
	}

	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the data being displayed in
	 * the GUI fields from the data stored in the in-memory
	 * Server object.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// First update the necessary information.
		description.setText(dataSet.getDescription());
		estFile.setText(dataSet.getPath());
		estFile.setEnabled(!lockPath);
		// Enable browse button only for local installs
		browse.setEnabled(!lockPath);
	}

	/**
	 * Method to intercept page change and trigger FASTA validation.
	 * 
	 * This method is invoked by the wizard when the user clicks the
	 * "Next" button. This method verifies that the FASTA file is
	 * valid. The actual detailed validation is done from a separate
	 * thread so that the GUI does not appear to hang.
	 * 
	 * @param dialog The wizard dialog which logically owns this page.
	 * @param currPage The current zero-based page number.
	 * @param nextPage The zero-based index number of the next page that
	 * this wizard will display.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		if (nextPage < currPage) {
			// The user want's to go back. That's OK.
			return true;
		}
		// Save the information for future reference.
		dataSet.setDescription(description.getText());
		if (lockPath) {
			// The EST path was uneditable. So assume it is good.
			return true;
		}
		// Normalize and do basic validation of file name.
		String estFileName = estFile.getText().trim();
		// Check if the specified EST file even exists and is readable.
		File est = new File(estFileName);
		if (!est.exists() || !est.isFile() || !est.canRead()) {
			JOptionPane.showMessageDialog(wizard, "<html>The specified EST file " +
					"does not exist (or is not readable).<br/>" +
					"Pleace specifiy a different source EST FASTA file.</html>",
					"Invalid EST File", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		// OK, the file at least exists.
		dataSet.setPath(est.getAbsolutePath());
		// Next let's do a more
		// through check on the FASTA file. This is done on a 
		// separate thread to ensure that the GUI does not lock up.

		// First Disable any changes to this page for now.
		Utilities.setEnabled(this, false);
		// Disable the buttons in in the wizard while this operation
		// is underway as we are going to spin off a separate thread
		// for this operation.
		dialog.setButtonStatus(0, 0, -1);
		// Launch the file verification thread.
		Thread verifThread = new Thread(this);
		// Flag this thread as a daemon because some I/O issues can
		// cause threads to hang around for a long time.
		verifThread.setDaemon(true);
		wizard.addThread(verifThread); // Track threads
		verifThread.start();

		// Prohibit page change for now.
		return false;
	}

	/**
	 * This is a helper method that is used to validate a given FASTA
	 * file. This method verifies the file is valid and then tries to
	 * load the EST data. If the data is loaded successfully this 
	 * method returns true. This method was primarily introduced to
	 * keep the code more readable.
	 */
	@Override
	public void run() {
		Exception exp = null;
		try {
			String estFileName = dataSet.getPath();
			// Before we proceed further, we verify that the specified
			// FASTA file has valid content.
			File file = new File(estFileName);
			if (!file.exists()) {
				throw new IOException("File " + estFileName + " was not found.");
			}
			// Check to ensure we have sufficient memory
			DataStore.get().memoryCheck(file, wizard);
			// Create input stream to read the file.
			FileInputStream fis = new FileInputStream(file);
			// Create a progress monitor in case the FASTA file is large
			String msg = "Verifying FASTA File: " + file.getName();
			ProgressMonitorInputStream pmis =
				new ProgressMonitorInputStream(this, msg, fis);
			pmis.getProgressMonitor().setNote("Please wait...");
			// Load the FASTA file.
			ESTList ests = ESTList.loadESTs(estFileName, pmis);
			if (ests.getESTs().size() < 1) {
				throw new Exception("The FASTA file did not have any ESTs");
			}
			// EST data was loaded successfully. That is good.
			wizard.setESTList(ests);
			// Do the page change in from the main AWT thread.
		} catch (Exception e) {
			exp = e;
		}
		// Show the results in the main AWT thread.
		final Exception error = exp;
		final JPanel page     = this;
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				// First re-enable the buttons
				wizard.setButtonStatus(1, 1, -1);
				// In addition re-enable the input fields so user
				// can change it if desired.
				Utilities.setEnabled(page, true);
				// Disable controls that should not be enabled
				estFile.setEnabled(!lockPath);
				browse.setEnabled(!lockPath);
				// Report error (if any)
				if (error != null) {
					JPanel msg = Utilities.collapsedMessage(INVALID_FASTA_MSG, 
							Utilities.toString(error));
					JOptionPane.showMessageDialog(wizard, msg, 
							"Invalid FASTA File", JOptionPane.ERROR_MESSAGE);
				} else {
					// Onto the next wizard page...
					wizard.changePage(wizard.getCurrentPage(), 
						wizard.getCurrentPage() + 1);
				}
			}
		});

		// Finally, remove the thread from the list of workers.
		wizard.removeThread(Thread.currentThread());
	}

	/**
	 * This method handles the tasks performed when the "Browse"
	 * button is clicked. It displays a JFileChooser and permits
	 * the user to select a local folder.
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		JFileChooser jfc = new JFileChooser();
		jfc.setDialogTitle("Choose EST File");
		jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		if (jfc.showDialog(this, "Use FASTA File") == JFileChooser.APPROVE_OPTION) {
			// Copy the chosen directory to the work space combo box.
			String wsPath = jfc.getSelectedFile().getAbsolutePath();
			// Set the selected item in this path.
			estFile.setText(wsPath);
		}
	}

	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final DataSetWizard wizard;

	/**
	 * Information about the actual data set entry being edited. This
	 * reference is set when this class is instantiated and is never
	 * changed during the life time.
	 */
	private final DataSet dataSet;

	/**
	 * Field to read/display the install path where PEACE is (or has
	 * been) installed on the given server.
	 */
	private JTextField estFile;

	/**
	 * Field to read and edit a brief description about the server.
	 * This information can be anything the user desires and is
	 * meaningful only to the user.
	 */
	private JTextArea description;

	/**
	 * If this flag is set to true, then this wizard permits the
	 * user only to change the install path on the server.
	 */
	private boolean lockPath;

	/**
	 * The browse button to be enabled only for local installs.
	 */
	private JButton browse;

	/**
	 * Just a simple error message that is displayed when connection
	 * to the remote host could not be successfully established.
	 */
	private static final String INVALID_FASTA_MSG = "<html>" + 
		"The specified EST file is not a valid FASTA file.<br/>" +
		"Please choose a different EST file.</html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752144L;
}
