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

package org.peace_tools.generic;

import javax.swing.Icon;
import javax.swing.tree.DefaultMutableTreeNode;

import org.peace_tools.views.PropertiesTreeMaker;

/**
 * A simple tree node that also permits an icon to be set.
 * 
 * This is an custom class that is used in the properties tree to
 * create a tree node with a predefined icon. Setting the icon
 * in the tree node makes developing a suitable cell renderer
 * easier.
 * 
 * @see PropertiesTreeMaker
 */
public class IconTreeNode extends DefaultMutableTreeNode {
	/**
	 * The icon to be used for this mutable tree node. This value
	 * is typically set in the constructor.
	 */
	private Icon icon = null;
	
	/**
	 * Constructor to create a tree node with an icon.
	 * 
	 * This constructor permits the creation of a tree node with a 
	 * suitable icon.
	 * 
	 * @param object The object to be associated with this tree node. This
	 * value is the user object set in the DefaultMutableTreeNode class.
	 * 
	 * @param icon The icon (if any) to be associated with this tree node.
	 * This value can be null.
	 */
	public IconTreeNode(Object object, Icon icon) {
		super(object);
		this.icon = icon;
	}
	
	/**
	 * Obtain the icon set for this tree node.
	 * 
	 * This method must be used to obtain the icon set for this tree node.
	 * 
	 * @return The icon set for this tree node. If a valid icon has not
	 * been set then this method returns null.
	 */
	public Icon getIcon() {
		return icon;
	}
	
	/**
	 * A generated serialization UID (just to keep the compiler happy).
	 */
	private static final long serialVersionUID = 6654295035253838886L;
}
