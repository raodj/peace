package org.peace_tools.decagon.helpers;

import java.awt.Image;
import java.awt.event.ActionListener;
import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JMenu;
import javax.swing.JMenuItem;

import org.peace_tools.decagon.jaxb.AssemblerDetails;
import org.peace_tools.generic.Log.LogLevel;
import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Utilities;

/**
 * Class to hold essential information for use by the GUI regarding
 * DECAGON Assembler Description XML (DADX) file.
 * 
 * This class serves the following two purposes for managing summary
 * information about DADX files (the summary information is frequently
 * used various components constituting DECAGON):
 * 
 * <ol>
 * 
 * <li>First it serves as a registry where the DADX files that were
 * loaded at startup time are summarized and managed as a list. The
 * list is used for creating menus and other operations.</li>
 * 
 * <li>This class serves to encapsulate the essential information from
 * DADX files that are frequenty used by various GUI components.</li>
 * 
 * </ol>
 * 
 */
public class DADXSummary {
	/**
	 * The full path to the DADX file from where the summary information
	 * for this object was loaded,.
	 */
	private final String dadxFilePath;
	
	/**
	 * A short name for the assembler to be displayed in main menu
	 * items.
	 */
	final private String name;
	
	/**
	 * A long form description for the assembler. The description may
	 * have HTML elements in it.
	 */
	final private String description;
	
	/**
	 * A 24x24 icon for the assembler. If the DADX file did not explicitly
	 * specify an icon, then this icon is a default icon supplied by
	 * DECAGON.  
	 */
	final private Icon largeIcon;
	
	/**
	 * A 16x16 icon for the assembler. If the DADX file did not explicitly
	 * specify an icon, then this icon is a default icon supplied by
	 * DECAGON.  
	 */
	final private Icon smallIcon;
	
	private static final HashMap<String, DADXSummary> assemblerList =
		new HashMap<String, DADXSummary>();
	
	/**
	 * The constructor.
	 * 
	 * The constructor initializes the core information managed by this
	 * class using the details about a given assembler. The details
	 * about an assembler are typically loaded from a DADX file 
	 * (see {@link DADXHelper#unmarshal(String, boolean)} for example).
	 * 
	 * @param filePath The path to the source DADX file from where the
	 * assembler description was loaded. This file name can be null.
	 * 
	 * @param src The AssemblerDetails object from where the summary
	 * information managed by this class are to be loaded. This parameter
	 * cannot be null.
	 */
	public DADXSummary(String dadxFilePath, AssemblerDetails src) {
		this.dadxFilePath = dadxFilePath;
		this.name         = src.getName();
		this.description  = src.getDescription();
		this.largeIcon    = getIcon(dadxFilePath, src.getIcon(), "images/24x24/EAST.png", 24);
		this.smallIcon    = getIcon(dadxFilePath, src.getIcon(), "images/16x16/EAST.png", 16);
	}
	
	/**
	 * Helper method to load an icon for an assembler and scale it to desired size.
	 * 
	 * This is a helper method that is invoked from the constructor to load
	 * an icon for a given assembler. This method attempts to load an icon
	 * specified in the DADX file (if any). Otherwise it uses a default
	 * assembler icon specified in the defIconFileName. If the icon is not
	 * of the specified size, then this method derives a suitably scaled
	 * icon from the original icon image.
	 * 
	 * @param dadxFilePath The path to the original DADX file (if any) which
	 * is used as a reference location to search for the icon file.
	 * 
	 * @param iconFileName The name of the icon file specified in the DADX file.
	 * This value can be null in which case the default icon is used.
	 * 
	 * @param defIconFileName The name of the default icon file to be used.
	 * 
	 * @param size The required size of the icon.
	 * 
	 * @return The icon object to be used.
	 */
	private Icon getIcon(final String dadxFilePath, 
			final String iconFileName, final String defIconFileName, int size) {
		URL       iconURL = null; // Path to the icon file.
		ImageIcon icon    = null; // The actual icon
		try {
			if (iconFileName != null) {
				// Obtain the fully resolved path to the icon
				File dadxFile = new File(dadxFilePath);
				iconURL = DADXHelper.resolveFilePath(iconFileName, dadxFile.getParent());
			}
		} catch (MalformedURLException e) {
			// Can't locate the user specified icon.
			UserLog.log(LogLevel.WARNING, "DADXSummary", "Unable to locate icon file: " + iconFileName);
		}
		// Load user-specified icon or the default icon.
		if (iconURL != null) {
			// Try and load user-specified icon.
			icon = Utilities.getIcon(iconURL.getPath());
		} else {
			// User specified icon not available. Load default.
			icon = Utilities.getIcon(defIconFileName);
		}
		// If icon was loaded scale it to the desired size.
		if (icon != null) {
			if ((icon.getIconHeight() != size) || (icon.getIconWidth() != size)) {
				// The icon needs to be scaled to the desired size.
				icon = new ImageIcon(icon.getImage().getScaledInstance(size, 
						size, Image.SCALE_SMOOTH));
			}
		}
		
		// Return the icon to be cached in the object.
		return icon;
	}
	
	/**
	 * Convenience method to create a main menu entry.
	 * 
	 * This method is invoked from DecagonMenuHelper to create entries
	 * for the various registered assemblers (entries in 
	 * {@link DADXSummary#assemblerList} hash map) in a given main menu. All of
	 * the entries are added with a specified action command to be processed
	 * by a given action listener.
	 * 
	 * @param menu The JMenu to which entries for the various assemblers
	 * are to be added.
	 * 
	 * @param actionCommandPrefxi The action command prefix string value to be 
	 * used with each menu item created by this method. The final action
	 * command is of the form actionCommandPrefix + ":" + name of assembler.
	 * (example: <code>runAssembler: EAST-1.0</code>).
	 * 
	 * @param listener The listener to be set for processing GUI events
	 * fired when a menu option is chosen.
	 */
	public static void addMenuItems(JMenu menu, String actionCommandPrefix, 
			ActionListener listener) {
		for(Entry<String, DADXSummary> entry: assemblerList.entrySet()) {
			menu.add(entry.getValue().getMenuItem(actionCommandPrefix, listener, true));
		}
	}
	
	/**
	 * Convenience method to create menu item using the summary information
	 * stored in this object.
	 * 
	 * @param actionCommandPrefxi The action command prefix string value to be 
	 * used with each menu item created by this method. The final action
	 * command is of the form actionCommandPrefix + ":" + name of assembler.
	 * (example: <code>runAssembler: EAST-1.0</code>).
	 * 
	 * @param listener The listener to be set for processing GUI events
	 * fired when a menu option is chosen.
	 * 
	 * @param mainMenu If this flag is true then a larger ICON is chosen.
	 * Otherwise a smaller icon is used in the menu entry.
	 * 
	 * @return This method returns a menu item created using the summary
	 * information stored in this object.
	 */
	public JMenuItem getMenuItem(String actionCommandPrefix, ActionListener listener, 
			boolean mainMenu) {
		// Create and return the main menu item
		JMenuItem item =
			Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, name,
				(mainMenu ? description : null),
				actionCommandPrefix + ":" + name, listener, null, null, true, false);
		item.setIcon(mainMenu ? largeIcon : smallIcon);
		item.setName(name);
		if (mainMenu) {
			item.setIconTextGap(6); // same as in Utilities.createMenuItem()
		}
		return item;
	}
	
	/**
	 * Obtain the small ICON associated with this assembler description.
	 * 
	 * This method can be used to obtain the a 16x16 icon associated
	 * with a given assembler's description. The icon is cached by
	 * this object and consequently is returned very fast (it is <b>not</b>
	 * loaded from disk each time this method is called).
	 * 
	 * @return The icon to be associated with the assembler's description.
	 * The returned object is never null.
	 */
	public Icon getSmallIcon() {
		return smallIcon;
	}
	
	/**
	 * Method to register assembler details for use by other parts of DECAGON.
	 * 
	 * This method must be used to "register" an assembler for use by other
	 * parts of DECAGON. This method essentially creates a DADXSummary 
	 * object and adds it to the list of entries in {@link #assemblerList}
	 * hash map for future reference. The {@link DADXSummary#name} field is
	 * used as the key. 
	 * 
	 * @param dadxFilePath The full path to the DADX file from where the
	 * assembler information was loaded. 
	 * 
	 * @param src The unmarshalled object containing the information from
	 * the DADX file. This information is used to create a DADXSummary
	 * object to be added to {@link #assemblerList}.
	 */
	public static void register(String dadxFilePath, AssemblerDetails src) {
		DADXSummary ds = new DADXSummary(dadxFilePath, src);
		assemblerList.put(ds.name, ds);
	}

	/**
	 * Obtain the summary information associated with a given genomic-assembler.
	 * 
	 * This method can be used to obtain the summary information about a
	 * given genomic assembler. This method is typically used by
	 * the {@link org.peace_tools.decagon.DecagonMenuHelper#runAssembler(String)} 
	 * method to obtain summary information about  given assembler in order
	 * to launch the appropriate wizard.
	 * 
	 * @param name The name of the assembler whose summary information is
	 * required. This string is the same name that was specified in the
	 * DADX file.
	 * 
	 * @return If the summary information was found then this method returns
	 * a valid summary object. Otherwise this method returns null.
	 */
	public static DADXSummary getSummary(final String name) {
		return assemblerList.get(name);
	}
	
	/**
	 * Obtain the path to the DADX file that contains the full assembler
	 * details.
	 * 
	 * @return The path to the DADX file from where the data summarized
	 * in this class was loaded. 
	 */
	public String getDadxFilePath() {
		return dadxFilePath;
	}
}
