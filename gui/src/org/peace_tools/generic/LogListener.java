package org.peace_tools.generic;
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

import java.util.EventListener;

/**
 * Interface for "Views" that are interested in receiving log changes.
 * 
 * This interface provides a mechanism for the Logging system to notify all
 * interested "View"s (Views as in the Model-View-Controller pattern) about
 * changes in logging information. Typically the view classes are child classes
 * of the LogPane class. Whenever new logs are cut, or if the log file name is
 * changed (by a view) then all views that implement this interface and register
 * themselves with the corresponding Log class will receive a call back
 * notification informing about change in the log information.
 */
public interface LogListener extends EventListener {
	/**
	 * This is the only call back method in this interface. This method is
	 * invoked whenever the information associated with the log changes or
	 * if the file name information changes.
	 * 
	 * @param logsChanged This parameter is true if log entries were added
	 * or removed from the logs. Typically this is the common scenario for
	 * the call back to occur.
	 * 
	 * @param fileNameChanged This parameter is true if the file name into
	 * which logs are written was changed by the user (via some view).
	 * 
	 * @param levelChanged This parameter is true if the current logging
	 * level set for the log was changed by the user (via some view).
	 */
	public void logChanged(final boolean logsChanged,
			final boolean fileNameChanged, final boolean levelChanged);
}
