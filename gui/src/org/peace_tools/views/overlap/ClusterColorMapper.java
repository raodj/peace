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

package org.peace_tools.views.overlap;

import java.awt.Color;


/**
 * A simple interface to map a clusterID to a Color object.
 * 
 * This interface was introduced to streamline the process of 
 * creating, managing, and using colors (defined by the user)
 * for fragments that are part of different clusters. This 
 * interface provides methods to set and get colors given a 
 * cluster ID. Currently only the OverlapView implements this
 * interface. The OverlapPanel uses the colors while the 
 * ClusterList class changes them.
 */
public interface ClusterColorMapper {
	/**
	 * Return the Color associated with a given cluster ID.
	 * 
	 * This method provides the currently used Color for a given
	 * cluster. 
	 *  
	 * @param clusterID The ID of the cluster whose current color
	 * is to be returned.
	 *  
	 * @return If the cluster ID is valid then this method returns
	 * a valid Color object. Otherwise it returns null. 
	 */
	public Color getColor(final Integer clusterID);

	/**
	 * Set a color to be used for rendering fragments in the given cluster.
	 * 
	 * This method must be used to set the color to be used for
	 * rendering the fragments in this cluster. Note that this method
	 * merely sets the color but does not trigger a refresh or an update.
	 * 
	 * @param clusterID The cluster ID whose color is to be set.
	 * @param color The color to be set for the given cluster ID.
	 */
	public void setColor(final Integer clusterID, final Color color);
}
