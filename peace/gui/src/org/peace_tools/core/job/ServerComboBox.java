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

package org.peace_tools.core.job;

import java.awt.Color;
import java.awt.Component;
import java.util.ArrayList;

import javax.swing.ImageIcon;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.ListCellRenderer;

import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.ServerList;
import org.peace_tools.workspace.Workspace;

/**
 * A GUI sub-component used permit user to select an pre-configured 
 * server entry for further use.
 * 
 * This component was primarily introduced to permit various wizards
 * to share common functionality of permitting a user to select
 * a server. This GUI component provides a customized JComboBox
 * that can is populated with Server entries. The server entries can
 * be filtered/restricted to set of a servers with specific properties
 * (like those that have PEACE, EAST, or DECAGON installed). In addition,
 * this class also provides a custom cell renderer that provides a more
 * detailed information about the server entries.
 * 
 */
public class ServerComboBox extends JComboBox implements ListCellRenderer {
	/** The default and only constructor.
	 * 
	 * This constructor is very straightforward and sets various instance
	 * variables to their default initial value.
	 */
	public ServerComboBox() {
		cellRenderer = new JLabel();
		setName("ServerComboBox");
		setBackground(Color.white);
		setRenderer(this);
	}
	
	/**
	 * Helper method to build/set the list of valid servers for listing in
	 * this combo-box.
	 * 
	 * This method is typically invoked from the wizard page (that logically
	 * owns this object). In addition, this method is also used by some classes
	 * to ensure that the workspace contains the necessary entries for their 
	 * correct operation (that this why this method is static).
	 *   
	 * This method essentially updates the list of server entries
	 * in the supplied combo-box (which is typically <code>this</code>) that 
	 * meet the given requirements (based on the values of the parameters). 
	 * 
	 * @param needPEACE If this flag is set to true, then only server entries
	 * that have a usable PEACE install are presented to the user as valid choices.
	 * 
	 * @param needEAST If this flag is set to true, then only server entries
	 * that have a usable EAST install are presented to the user as valid choices.
	 * 
	 * @return This method returns the number of entries that are finally present
	 * in the supplied server list.
	 */
	public static int getServerList(final boolean needPEACE, final boolean needEAST,
			JComboBox serverList) {
		final int srvrChoiceType = (needPEACE ? 0 : 1) + (needEAST ? 2 : 0); 
		// Clear any previous entries
		serverList.removeAllItems();
		// If we don't have a workspace entry to work with get out
		if (Workspace.get() == null) {
			return 0;
		}
		// Populate the combo-box with server entries from work space.
		ServerList list = Workspace.get().getServerList();
		ArrayList<Server> servers = list.getServers();
		// On add servers that are a successful PEACE install.
		for(Server srvr: servers) {
			boolean validServer = false;
			switch (srvrChoiceType) {			
			case 2: 
			case 3: validServer = srvr.hasEASTInstalled();
			break;
			default: validServer = true;
			}
			if ((!Server.ServerStatusType.GOOD.equals(srvr.getStatus()) && 
					!Server.ServerStatusType.CONNECT_FAILED.equals(srvr.getStatus())) ||
					!validServer) {
				// Server is not usable as its status is not good or we are
				// currently unable to connect to it.
				continue;
			}			
			serverList.addItem(srvr);
		}
		return serverList.getItemCount();
	}

	/**
	 * This method is called by core Java GUI system whenever it
	 * needs to render a data set or MST entry for the user to
	 * choose from the {@link #dataSetList} combo box. This
	 * method appropriately updates the {@link #cellRenderer}
	 * label (created in the constructor) and return the
	 * label back. See JavaDoc on this method for additional
	 * details.
	 */
	@Override
	public Component getListCellRendererComponent(JList list, Object value,
			int index, boolean isSelected, boolean cellHasFocus) {
		// Set the default background color depending on selection
		cellRenderer.setOpaque(true);
        if (isSelected) {
            cellRenderer.setBackground(list.getSelectionBackground());
            cellRenderer.setForeground(list.getSelectionForeground());
        } else {
        	cellRenderer.setBackground(Color.white);
        	cellRenderer.setForeground(list.getForeground());
        }
        if (value != null) {
        	// Set the text to be rendered based on the current server.
        	Server srvr = (Server) value;
        	String srvrInfo = "<html><b>" + srvr.getName() + "</b><br/>" +
        	"<font size=\"-2\">" + Utilities.trim(srvr.getDescription(), 50) + 
        	"</font></html>";
        	cellRenderer.setText(srvrInfo);
        	//Set the icon to make things look pretty and meaningful
        	int iconIndex = 0 + (srvr.hasEASTInstalled() ? 1 : 0) + (srvr.hasDECAGONInstalled() ? 1 : 0);
        	cellRenderer.setIcon(ServerIcons[iconIndex]);
        } else {
        	// We don't have a valid list item to render
        	cellRenderer.setText("");
        	cellRenderer.setIcon(null);
        }
        return cellRenderer;
	}
	
	/**
	 * This label is used by this class to render the Server entries
	 * presented by this class to the user. This cell renderer is created
	 * once (in the constructor) and reused each time the 
	 * {@link #getListCellRendererComponent(JList, Object, int, boolean, boolean)}
	 * method is called (by JComboBox parent class)
	 */
	private final JLabel cellRenderer;
	
	/**
	 * The list of icons that are displayed along with various servers listed
	 * by this class. The icons are primarily used to provide some quick
	 * additional information to the user in a pretty way.
	 */
	private static final ImageIcon ServerIcons[] = {
		Utilities.getIcon("images/24x24/Server.png"), // PEACE only server
		Utilities.getIcon("images/24x24/Server.png"), // PEACE + EAST server
		Utilities.getIcon("images/24x24/Server.png")  // PEACE + EAST + DECAGON server
	};
	
	/**
	 * A generated serialization GUID (primarily to keep compiler happy
	 * and prevent a warning from being generated)
	 */
	private static final long serialVersionUID = -5425868665130124990L;
}
