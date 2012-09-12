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
import java.util.ArrayList;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.peace_tools.generic.BackgroundTask;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;

/**
 * Helper method to load a set of parameters from XML files.
 *
 * This is a convenience class that can be used by GUI components to load
 * a given set of parameter definitions from XML files. This class
 * provides the user with a dialog box that shows the progress of
 * loading the parameters. If errors occurs when loading files this
 * class displays a suitable error message.
 *
 */
public class ParameterSetLoader implements BackgroundTask.UserTask {
	/**
	 * The names of the file from where the parameters are to be loaded.
	 * This value is initialized in the constructor. This array should not
	 * be null or an empty array.
	 */
	private final String[] paramFileNames;

	/**
	 * The list of parameters that were loaded by this class. This value is
	 * returned by the {@link #loadParameterSets()} method.
	 */
	final ArrayList<ParameterSetGUIHelper> paramList = new ArrayList<ParameterSetGUIHelper>();

	/**
	 * Create a loader to load one or more parameter sets from XML file(s).
	 * 
	 * This constructor must be used to create an instance of this class
	 * to load a one or more parameter sets from a given set of XML files(s).
	 * The constructor does not perform the task of loading the files.
	 * It merely setups the GUI components to perform the task. The actual
	 * loading is done when the {@link #loadParameterSets()} method is called.
	 * 
	 * @param paramFileNames The file names of the parameter files from where
	 * the data is to be loaded. This array must have at least one entry in it.
	 * The data is loaded by calling {@link ParameterSetGUIHelper#unmarshal(String)}.
	 */
	public ParameterSetLoader(final String[] paramFileNames) {
		this.paramFileNames = paramFileNames;
	}
	
	/**
	 * Convenience method to load the parameters from the XML file(s).
	 * 
	 * This method must be used to trigger the actual process of loading
	 * parameters from the given list of XML files. This method displays
	 * a dialog in the foreground and performs the actual task of loading
	 * parameters in the background. This method displays any errors
	 * that may occur during parameter loading. If parameters are loaded
	 * successfully, then this method does not require any specific
	 * interactions from the user (other than informing the user about
	 * progress when loading large parameter sets).  
	 * 
	 * @param parent The parent window that should be used to display
	 * dialog boxes.
	 * 
	 * @return If the parameters are to be loaded successfully then this
	 * method returns null. If parameters are loaded successfully then
	 * this method returns an array list containing the loaded parameters
	 * for further use by GUI. 
	 * 
	 * <b>NOTE:</b> It is guaranteed that the parameters will be returned in 
	 * the same order in which the files were specified when this class
	 * was instantiated.
	 */
	public ArrayList<ParameterSetGUIHelper> loadParameterSets(Component parent) {
		// Create the background task object to perform database creation in 
		// the background.
		final JLabel infoLabel = new JLabel("Loading parameter definitions. Please wait...",
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		// Create step labels based on the list of file names to load.
		final String[] stepLabels = new String[paramFileNames.length];
		for(int step = 0; (step < paramFileNames.length); step++) {
			stepLabels[step] = "Loading: " + paramFileNames[step];
		}
		// Create the background task (that shows a dialog box)
		BackgroundTask bt = new BackgroundTask(this, infoLabel,
				true, stepLabels, true, true, true, 1, false, false, false, null);
		// Kick start the background operations
		bt.start(true, parent, "Loading Parameter Definitions", 
				Utilities.getIcon("images/16x16/Import.png"));
		// Wait for the list to be loaded.
		try {
			// The following get() call waits for background task to finish.
			if (bt.get() != null) {
				// An error occurred! Return invalid parameter list.
				return null;
			}
			// Parameters were successfully loaded.
			return paramList;
		} catch (Exception e) {
			// Some unforeseen exception occurred.
			ProgrammerLog.log(e);
			// Display error to the user.
			final String htmlErrMsg = "<html>Error loading parameter definitions from XML file.<br>"
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
		// Load parameters one at a time from the given list
		for(int param = 0; (param < paramFileNames.length); param++) {
			ParameterSetGUIHelper psgh = new ParameterSetGUIHelper();
			// Load data from the given XML file. This call can generate
			// several different types of exceptions.
			psgh.unmarshal(paramFileNames[param]);
			// This parameter was successfully loaded. Add it to outgoing list.
			paramList.add(psgh);
			ProgrammerLog.log("ParameterSetLoader successfully loaded parameter file: " + paramFileNames[param]);
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
			bTask.showMessage("", "Error loading parameters from an XML file");
		}
		// In either case we want to close the dialog box after this one.
		bTask.disposeDialog();
	}

	@Override
	public void dialogClosed(BackgroundTask bTask) {
		// Nothing to be done here for now.
	}
}
