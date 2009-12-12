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
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import org.peace_tools.generic.ProgrammerLog;

/**
 * This log pane customizes the features of the LogPane to display the 
 * programmer logs in a free form manner. This log can be used by programmers
 * to log information and troubleshoot problems. 
 * The log can contain any message that the system wishes to log to the user.
 * 
 */
public class ProgrammerLogPane extends LogPane {
	/**
	 * A constructor for the ProgrammerLogPane. It sets up the JTextArea 
	 * the log is written to.
	 * 
	 * @param logLevel
	 *            The Log level to display in the JTextArea
	 * @param name
	 *            The name of the file to also log to
	 */
	public ProgrammerLogPane() {
		// The base class sets up the core panel and the tool bar at top.
		super(ProgrammerLog.getLog(), false);
		// Setup the JTextArea that is used to display programmer log
		logArea = new JTextArea();
		// Setup some properties.
		// Place the log in a scroll pane so that logs can be scrolled
		JScrollPane scroller = new JScrollPane(logArea);
		add(scroller, BorderLayout.CENTER);
		// Initially update the data in programmer log
		ProgrammerLog pLog = ProgrammerLog.getLog();
		logArea.setText(pLog.getLogEntries().toString());
	}

	/**
	 * This method overrides the default implementation in the base 
	 * class to update the programmer logs displayed by this component.
	 * It still calls the base class to have the default operations
	 * performed.
	 */
	@Override
	public void logChanged(boolean logsChanged, boolean fileNameChanged,
			boolean levelChanged) {
		// Let the base class set up the basic information.
		super.logChanged(logsChanged, fileNameChanged, levelChanged);
		// Update the text in the log if logs have changed
		if (logsChanged) {
			ProgrammerLog pLog = ProgrammerLog.getLog();
			logArea.setText(pLog.getLogEntries().toString());
		}
	}
	
	/**
	 * The table to which we write our log into. This is the visual component
	 * that provides the necessary display to the user in a tabular form.
	 */
	private JTextArea logArea;
		
	/**
	 * The general serial version UID to enable serialization of this class as
	 * per Java requirements.
	 */
	private static final long serialVersionUID = 2213371948161753574L;
}
