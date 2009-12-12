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

package org.peace_tools.workspace;

/**
 * A simple interface class that is used to report changes occuring to
 * the currently active, global Workspace. The changes can be inserts (or
 * addition of a new entry), updates (change in status or information 
 * associated with an entry), or deletes (removal of existing entries).
 * The events are generated for various entries in a workspace, 
 * including: data set, mst data, cluster data, job, and server entries.
 */
public interface WorkspaceListener {
	/**
	 * This method is invoked on all workspace listeners registered with
	 * the current Workspace. The event contains the necessary information
	 * to report the change that has occured to a workspace.
	 *  
	 * @param event The event that contains the information regarding the
	 * change that has occured to the workspace.
	 */
	public void workspaceChanged(WorkspaceEvent event);
}
