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
import java.awt.Component;
import javax.swing.Icon;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import org.peace_tools.generic.Log;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;

/**
 * LogPane is a simple way for the system to present the user with useful
 * information. The LogPane contains a JTable which shows a "log." The log can
 * contain any message that the system wishes to log to the user. This class
 * extends the common features provided by the LogPane base class to display log
 * entries in a more user-friendly form.
 * 
 * <p>
 * <b>Note:</b> Logs for programming aids and debugging (that would be
 * pertinent for a programmer to view and understand) must be generated in the
 * free form ProgrammerLog and not in the UserLog.
 * </p>
 */
public class UserLogPane extends LogPane {
	/**
	 * A constructor for the LogPane. It sets up the JTextArea the log is
	 * written to.
	 */
	public UserLogPane() {
		// The base class sets up the core panel and the tool bar at top.
		super(UserLog.getLog(), true);
		// Setup the JTable that is used to display user logs
		log = new JTable(UserLog.getLog().getLogEntries()) {
			private static final long serialVersionUID = 1430052270991586572L;
			
			@Override
			public boolean isCellEditable(int row, int column) {
				return false;
			}
			@Override
			public TableCellRenderer getCellRenderer(int row, int column) {
				if (column == 0) {
					return TimeStampRenderer; 
				}
				return super.getCellRenderer(row, column);
			}
		};
		// Set some of the table properties
		log.setBorder(null);
		log.setShowHorizontalLines(true);
		log.setFillsViewportHeight(true);
		log.setDragEnabled(false);
		log.setDropTarget(null);
		log.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		// Setup some properties.
		TableColumnModel tcm = log.getColumnModel();
		tcm.getColumn(0).setPreferredWidth(75);
		tcm.getColumn(1).setPreferredWidth(75);
		tcm.getColumn(2).setPreferredWidth(600);
		// Place the log in a scroll pane so that logs can be scrolled
		JScrollPane scroller = new JScrollPane(log);
		scroller.getViewport().setBackground(log.getBackground());
		add(scroller, BorderLayout.CENTER);
	}

	/**
	 * The table to which we write our log into. This is the visual component
	 * that provides the necessary display to the user in a tabular form.
	 */
	private JTable log;
	
	/**
	 * The set of icons/images used by the StatusTimeStamp class to
	 * provide an aggregate view of both status and time stamps for
	 * every entry in the user log.
	 */
	private static Icon LevelIcons[] = {
		Utilities.getIcon("images/16x16/Information.png"),
		Utilities.getIcon("images/16x16/Information.png"),
		Utilities.getIcon("images/16x16/Warning.png"),
		Utilities.getIcon("images/16x16/Error.png"),
	};

	private class StatusTimeStamp extends JLabel
    implements TableCellRenderer {
		/**
		 * The general serial version UID to enable serialization of this class as
		 * per Java requirements.
		 */
		private static final long serialVersionUID = 5654108539980884223L;
		
		/** Default constructor.
		 * 
		 * The default constructor merely sets the default font to be used
		 * by the label to be non-bold.
		 */
		public StatusTimeStamp() {
			// Make the font non-bold.
			Utilities.adjustFont(this, 0, 8, -1);
			setOpaque(true);
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			String data = (String) value;
			// The first character is the log level code. Use that to 
			// determine the ICON to set and rest of the data (other than
			// the first character) is the time stamp.
			setText(data.substring(1));
			Log.LogLevel level = Log.decode(data.charAt(0));
			setIcon(LevelIcons[level.ordinal()]);
			setBackground(isSelected ? table.getSelectionBackground() : table.getBackground());
			return this;
		}
	}

	/**
	 * An instance of the cell renderer that is used to render the
	 * first column in the table with a suitable status icon.
	 */
	private final StatusTimeStamp TimeStampRenderer = new StatusTimeStamp();
	
	/**
	 * The general serial version UID to enable serialization of this class as
	 * per Java requirements.
	 */
	private static final long serialVersionUID = 2213371948161753574L;
}
