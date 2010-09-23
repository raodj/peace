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
 * Interface for views (GUI components) dialog that are interested in 
 * receiving "Find" requests from the FindDialog.
 * 
 * This interface provides a straightforward mechanism for dispatching find
 * (or search) requests from the generic find dialog. Each GUI component that
 * is interested in receiving find events must implement this interface and
 * register itself with the find dialog to receive search requests.
 */
public interface FindListener extends EventListener {
	/**
	 * This is the only call back method in this interface. This method is
	 * invoked whenever the user clicks the "Find" button and all the information
	 * that the user has entered is deemed valid.
	 * 
	 * @param event The find event that contains all the necessary information
	 * required to execute a find operation.
	 * 
	 * @return This method must return true if the find operation actually
	 * found information. Otherwise it returns false.
	 */
	public boolean find(FindEvent event);
}
