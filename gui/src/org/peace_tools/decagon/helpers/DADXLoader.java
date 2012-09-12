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

package org.peace_tools.decagon.helpers;

import java.awt.Component;
import java.io.File;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.peace_tools.generic.BackgroundTask;
import org.peace_tools.generic.Log.LogLevel;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;

/**
 * Helper method to load a DADX file in the background while displaying
 * progress in the foreground.
 *
 * This is a convenience class that can be used by GUI components to load
 * a given DADX file. This class provides the user with a dialog box that
 * shows the progress of loading the DADX file. If errors occurs when 
 * loading files this class displays a suitable error message. The actual
 * unmarshalling is done via call to {@link DADXHelper#unmarshal(String)}
 * method.
 *
 */
public class DADXLoader implements BackgroundTask.UserTask {
	/**
	 * The names of the DADX file to be loaded.
	 * This value is initialized in the constructor.
	 */
	private final String dadxFileNames[];

	/**
	 * The DADXHelper object that is used to actually load the XML and
	 * unmarshall it into an in-memory data structure for further processing.
	 * This value is returned by the {@link #loadDADX()} method.
	 */
	final DADXHelper dadxHelper = new DADXHelper();

	/**
	 * Create a loader to load DECAGON Assembler Details XML (DADX) file. 
	 * 
	 * This constructor must be used to create an instance of this class
	 * to load assembler details from a given DADX files. The constructor 
	 * does not perform the task of loading the files.
	 * It merely setups the GUI components to perform the task. The actual
	 * loading is done when the {@link #loadDADX()} method is called.
	 * 
	 * @param dadxFileName An array of file names of the DADX files to 
	 * be unmarshalled. This value cannot be null and this array must 
	 * contain at least one entry. The paths can be to an internal resource
	 * file that is available in PEACE jar file or an external file on disk.
	 * 
	 * @see {@link org.peace_tools.decagon.DecagonMenuHelper#runAssembler(String)}
	 */
	public DADXLoader(final String dadxFileNames[]) {
		this.dadxFileNames = dadxFileNames;
	}
	
	/**
	 * Method to load (or unmarshall) the specified DADX file.
	 * 
	 * This method must be used to trigger the actual process of loading
	 * assembler details from the given DADX file. This method displays
	 * a dialog in the foreground and performs the actual task of loading
	 * the DADX file in the background. This method displays any errors
	 * that may occur during loading. If the DADX data is loaded
	 * successfully, then this method does not require any specific
	 * interactions from the user.  
	 * 
	 * @param parent The parent window that should be used to display
	 * dialog boxes.
	 * 
	 * @return If the DADX file is <b>not</b>loaded successfully then this
	 * method returns null. If parameters are loaded successfully then
	 * this method returns the DADXHelper object that can be used to 
	 * perform further operations using the loaded data. 
	 */
	public DADXHelper loadDADX(Component parent) {
		// Create the background task object to perform database creation in 
		// the background.
		final JLabel infoLabel = new JLabel("Loading assembler details. Please wait...",
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		// Create the background task (that shows a dialog box)
		BackgroundTask bt = new BackgroundTask(this, infoLabel,
				true, dadxFileNames, true, false, true, 1, false, false, false, null);
		// Kick start the background operations
		File tmpFile = new File(dadxFileNames[0]);
		bt.start(true, parent, "Loading " + tmpFile.getName(), 
				Utilities.getIcon("images/16x16/Import.png"));
		// Wait for the list to be loaded.
		try {
			// The following get() call waits for background task to finish.
			if (bt.get() != null) {
				// An error occurred! Return invalid parameter list.
				return null;
			}
			// Parameters were successfully loaded.
			return dadxHelper;
		} catch (Exception e) {
			// Some unforeseen exception occurred.
			ProgrammerLog.log(e);
			// Display error to the user.
			final String htmlErrMsg = "<html>Error loading assembler description from XML file.<br>"
					+ "<dl><dt><b>Error:</b></dt><dd>" + Utilities.wrapStringToHTML(e.getMessage(), 65)
					+ "</dd></dl></html>";
				JPanel msgPanel = Utilities.collapsedMessage(htmlErrMsg, Utilities.toString(e));
			// Display the message
			JOptionPane.showMessageDialog(parent, msgPanel, "Error Loading Parameters",
					JOptionPane.ERROR_MESSAGE);
		}
		// Return invalid parameter list to indicate some problem occurred.
		return null;
	}
	
	@Override
	public void run(BackgroundTask bTask) throws Exception {
		// Load the DADX filed via DADX helper.
		for(String fileName: dadxFileNames) {
			dadxHelper.unmarshal(fileName, false, true);
			// Cut log to indicate some feed back to user 
			UserLog.log(LogLevel.NOTICE, "DADXHelper", 
					"Successfully loaded data from DADX file: " + fileName);
			// Let user know we have made some progress.
			bTask.updateProgress();
		}
	}

	@Override
	public void done(BackgroundTask bTask) {
		// Check and display an error message only if error occurs.
		// We intentionally don't display a success message here to ensure
		// smooth user experience.
		if (bTask.getResult() != null) {
			bTask.showMessage("", "Error loading assembler details from DADX file");
		}
		// In either case we want to close the dialog box after this one.
		bTask.disposeDialog();
	}

	@Override
	public void dialogClosed(BackgroundTask bTask) {
		// Nothing to be done here for now.
	}
}
