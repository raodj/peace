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

package org.peace_tools.decagon;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.AbstractButton;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TreeSelectionListener;

import org.peace_tools.core.AbstractMenuHelper;
import org.peace_tools.core.MainFrame;
import org.peace_tools.core.PEACEProperties;
import org.peace_tools.decagon.configMgr.DecagonConfigurationManager;
import org.peace_tools.decagon.djc.DecagonJobCreator;
import org.peace_tools.decagon.helpers.DADXHelper;
import org.peace_tools.decagon.helpers.DADXLoader;
import org.peace_tools.decagon.helpers.DADXSummary;
import org.peace_tools.decagon.sdg.SyntheticDataSetGenerator;
import org.peace_tools.generic.Log.LogLevel;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;

/**
 * The DECAGON menu helper for PEACE-GUI.
 * 
 * This class encapsulates the code related to the operations performed by
 * various menu items in the "DECAGON" menu. This class is typically created once
 * from the MainFrame.createMenus() method. The primary motivation for
 * introducing sub-menu handler classes is to improve code organization and
 * minimize code clutter. Note that this class is essentially an event
 * handler that is set on the various menu items created by this class. 
 * This helper class also provides a  createDecagonMenu() method that 
 * actually creates the "DECAGON" menu.
 */
public class DecagonMenuHelper extends AbstractMenuHelper
implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "Help" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public DecagonMenuHelper(MainFrame mainFrame) {
		super(HelperType.DECAGON_MENU, mainFrame);
		
	}

	/**
	 * A static constant string used for action command prefixes for few
	 * menu entries (pertaining to launching the wizard for a given
	 * genomic assembler) used in this class. 
	 */
	private static final String RUN_ASSEMBLER_CMD = "runAssembler";
	
	/**
	 * Helper method to create the DECAGON menu.
	 * 
	 * <p>This method performs the actual task of creating the "DECAGON" menu. This
	 * method has been introduced to organize all the methods related to the
	 * "DECAGON" menu into a single class. This method is invoked from the MainMenu
	 * class to create the "Server" menu.<p/>
	 * 
	 * <p>This method tracks menu items whose status is to be modified in the context sensitive
	 * list maintained in the base class. These menu items are not really
	 * context sensitive but more configuration sensitive. But we are essentially
	 * reusing the feature available in the base class correctly.<p>
	 * 
	 * @param toolbar The toolbar to which frequently used shortcuts can
	 * be added typically in the form of icons. If the toolbar is null, then
	 * shortcuts are not added.
	 */
	public JMenu createDecagonMenu(JToolBar toolbar) {
		// Create just the menu entry.
		JMenu decagonMenu = new JMenu("DECAGON  ");
		// Create various menu items.
		assert( decagonMenu != null );
		// First create entries for various known assemblers whose summary
		// information is cached.
		JMenu asmList = (JMenu) Utilities.createMenuItem(Utilities.MenuItemKind.SUB_MENU_ITEM, 
				"Run a genomic assembler", "Assembles a selected data set using a given assembler",
				null, null, "images/24x24/EAST.png", null, true, false);
		DADXSummary.addMenuItems(asmList, RUN_ASSEMBLER_CMD, this);
		decagonMenu.add(asmList);
		
	    // Next create and add the new file creation menu options.
		JMenuItem multiSimMI = getMenuItem(ActionType.GENERATE_READS_VIA_MULTISIM, true);
		decagonMenu.add(multiSimMI);
		this.contextItemList.add(multiSimMI);
		decagonMenu.addSeparator();
		
		// Create a menu to configure the various DECAGON options
		decagonMenu.add(getMenuItem(ActionType.CONFIGURE_DECAGON, true));
		
		// Add tool bar entry to add a server entry.
		if (toolbar != null) {
			// Add some of the tools that we anticipate users to work
			// with frequently to the tool bar.
		}
		
		// Enable/disable menu items based on current configuration.
		updateMenuItems();
		return decagonMenu;
	}
	
	@Override
	public ActionListener getActionListener() {
		return this;
	}

	@Override
	public JMenuItem getMenuItem(ActionType actionType, boolean mainMenu) {
		int index = actionType.ordinal() - ActionType.GENERATE_READS_VIA_MULTISIM.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainMenu ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";
		// Create and return the main menu item
		JMenuItem item =
			Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, MenuTitles[index],
				(mainMenu ? MenuSubTitles[index] : null),
				ActionCmds[index], this, IconPath, 
				null, true, false);
		return item;
	}
	
	@Override
	public AbstractButton getTool(ActionType actionType, boolean mainToolBar) {
		int index = actionType.ordinal() - ActionType.CONFIGURE_DECAGON.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainToolBar ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";

		AbstractButton item =
			Utilities.createToolButton(IconPath, null, ActionCmds[index], this, 
				MenuSubTitles[index], true);
		if (index > 0) {
			// Track context sensitive entries
			contextItemList.add(item);
			item.setEnabled(false);
		}
		return item;
	}

	@Override
	public void actionPerformed(ActionEvent event) {
		String cmd = event.getActionCommand();
		if (ActionCmds[0].equals(cmd)) {
			// Synthetic DataSet Generator
			 SyntheticDataSetGenerator sdg = SyntheticDataSetGenerator.create(mainFrame);
			 if (sdg != null) {
				 Utilities.centerPanel(mainFrame, sdg);
				 sdg.showWizard("http://www.peace-tools.org/downloads/manual.pdf#page=30");
			 }
		} else if (ActionCmds[1].equals(cmd)) {
			// Configure option...
			DecagonConfigurationManager dcm = new DecagonConfigurationManager(mainFrame);
			 Utilities.centerPanel(mainFrame, dcm);
			 dcm.showWizard("http://www.peace-tools.org/downloads/manual.pdf#page=30");
		} else if ((cmd != null) && cmd.startsWith(RUN_ASSEMBLER_CMD)) {
			runAssembler(cmd.substring(RUN_ASSEMBLER_CMD.length() + 1), true);
		}
	}
	
	/**
	 * Convenience method to load a DECAGON Assembler Description XML (DADX) 
	 * file and optionally the DECAGON Analyzer DADX file.  
	 * 
	 * This method is a convenience method that can be used to load a
	 * given DADX file. This method can optionally also load the DADX 
	 * file that provides the DECAGON analyzer definition. The DECAGON
	 * analyzer runs the analyzer pipeline that generates scores and
	 * statistics regarding the quality of contigs generated by the
	 * assembler pipeline.  
	 * 
	 * @param assemblerName The name of the assembler whose DADX file
	 * description is to be loaded by this method. The assembler name
	 * must match the name of an earlier DADX file that was loaded
	 * during PEACE start up.
	 * 
	 * @param addDecAnalyzer A boolean flag that indicates if the 
	 * additional DECAGON analyzer DADX file is to be loaded.
	 */
	public void runAssembler(final String assemblerName, boolean addDecAnalyzer) {
		DADXSummary summary = DADXSummary.getSummary(assemblerName);
		if (summary == null) {
			JOptionPane.showMessageDialog(mainFrame, 
					"Unable to find information for assembler " + assemblerName, 
					"Invalid Assembler", JOptionPane.ERROR_MESSAGE);
			return;
		}
		// Use the summary information for the assembler to load full details
		// about the assembler and DECAGON analyzer job for use by the wizard.
		String fileNames[] = null;
		if (addDecAnalyzer) {
			final String AnalyzerFileName = "decagon/DECAGONAnalyzer.dadx";
			fileNames = new String[]{summary.getDadxFilePath(), AnalyzerFileName};
		} else {
			fileNames = new String[]{summary.getDadxFilePath()};
		}
		DADXLoader dadxLoader = new DADXLoader(fileNames);
		DADXHelper dadxHelper = dadxLoader.loadDADX(mainFrame);
		if (dadxHelper != null) {
			UserLog.log(LogLevel.NOTICE, "DECAGON", "Launching wizard for assembler " + assemblerName);
			// Create the wizard (if various checks pass)
			DecagonJobCreator djc = DecagonJobCreator.create(dadxHelper, mainFrame);
			if (djc != null) {
				djc.showWizard("http://www.peace-tools.org/downloads/manual.pdf#page=30");
			}
		}
	}
	
	/**
	 * Helper method to enable/disable menu items depending on the current DECAGON
	 * configuration settings.
	 * 
	 * This method is typically called from {@link DecagonConfigurationManager#done(boolean)}
	 * method to enable or disable menu items based on the current 
	 * DECAGON configuration settings. This method uses corresponding entries
	 * in the global PEACEProperties to enable/disable menu items.
	 */
	public void updateMenuItems() {
		// Shortcut to global PEACE properties.
		// Short cut to the current global properties
		PEACEProperties props = PEACEProperties.get();
		// Enable/disable Metasim use to generate synthetic data sets
		final boolean enableMetaSim = props.getBooleanProperty(PropertyKeys.METASIM_ENABLED, false);
		setEnabled(ActionCmds[0], enableMetaSim);		
	}
	
	@Override
	public ListSelectionListener getListSelectionListener(JTable table) {
		return null;
	}
	
	@Override
	public TreeSelectionListener getTreeSelectionListener(JTree tree) {
		return null;
	}
	
	/**
	 * The strings for each menu title created by this helper. The list of
	 * values are organized in the same order as the ordinal values of 
	 * the {@link AbstractMenuHelper#ActionType} enumeration. Preserving 
	 * the order is important.
	 */
	private static final String MenuTitles[] = {
		"Generate Synthetic DataSet",
		"Configure Properties..."
	};
	
	/**
	 * The icon file names for each menu title created by this helper. 
	 * The list of values are organized in the same order as the ordinal 
	 * values of the ActionType enumeration. Preserving the order is
	 * important. Note that a prefix directory (such as: images/16x16)
	 * and a suffix extension (.png) is added when tools or menu items
	 * are created.
	 */
	private static final String IconNames[] = {
		"ServerAdd", 
		"Baton",
	};
	
	/**
	 * The strings for the action commands generated by the various menu items
	 * in the help menu. The list of values are organized in the same order as
	 * the ordinal values of the ActionType enumeration. Preserving the order 
	 * is important.
	 */
	private static final String ActionCmds[] = {
		"RunMultiSim", 
		"Configure",
	};
	
	/**
	 * The various sub menu titles that are used in the main menu. The
	 * sub menu titles are used to provide the user with a bit more 
	 * verbose description on the action that will be performed by a given
	 * menu item. The list of values are organized in the same order as 
	 *  the ordinal values of the ActionType enumeration. Preserving the 
	 *  order is important.
	 */
	private static final String MenuSubTitles[] = {
		"Generate synthetic dataset (via MultiSim) for evaluating assemblers", 
		"Configure various DECAGON features to streamline operations",
	};
}
