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

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.net.URLConnection;
import java.util.List;

import javax.swing.AbstractButton;
import javax.swing.Box;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.KeyStroke;
import javax.swing.ProgressMonitor;
import javax.swing.SwingUtilities;

/**
 * This class acts as a helper class and contains STATIC methods (only) that
 * acts as helpers to handle some of the common, often performed tasks in the
 * presentation manager
 */
public class Utilities {
	/**
	 * The default prefix that is added to the path of the files to be loaded by
	 * the various methods in this class. This prefix path enables organizing
	 * and collating files appropriately.
	 */
	public static final String PATH_PREFIX = "../../../";

	/**
	 * A simple string that has 4 HTML spaces. This string is typically used
	 * to indent HTML content a bit to organize them better.
	 */
	public static final String HTML_4SPACES = "&nbsp; &nbsp; &nbsp; &nbsp;";
	
	/** Constants that can be supplied to the createMenuItem method 
	 * to create different types of menu items to be added to a menu
	 */
	public enum MenuItemKind {
		/** Creates a routine menu item with an icon, main text, and sub-text */
		MENU_ITEM, 
		/** Create a check-box type menu item that contains a small check-box
		 * that the user can select/unselect.
		 */
		CHECK_BOX_ITEM, 
		/**
		 * Create a radio button item that can be used to select one of 
		 * several options. Refer to the standard Javadoc on 
		 * {@link JMenuItem} for details on how to use radio button item.
		 */
		RADIO_BUTTON_ITEM, 
		/**
		 * This is a sub-menu item that can contain other sub-menus and
		 * menu items. It looks similar to a menu-item except it can 
		 * contain sub-items.
		 */
		SUB_MENU_ITEM
	};
	
	/**
	 * This method is a convenience method that provides a simple mechanism for
	 * creating a JMenuItem, given the parameters.
	 * 
	 * @param itemKind
	 *            This is a enumeration value which can be MENU_ITEM,
	 *            CHECK_BOX_ITEM, RADIO_BUTTON_ITEM, or SUB_MENU_ITEM
	 * 
	 * @param itemTitle
	 *            String representing the menu option. Note: The string may
	 *            contain a "_" (underscore) character preceeding the letter
	 *            that is suppoed to be used as an accelerator for this item.
	 *            However, if an explicit key stroke is specified (parameter 3
	 *            to this method call) then this underscore will be ignored. In
	 *            either case the "_" will be stripped out of the option! The
	 *            string cannot be "null"
	 * 
	 * @param command
	 *            The string representing the action command that will be
	 *            generated when this menu option is slected!
	 * @param al
	 *            The ActionListener to register with the Menu Item.
	 * 
	 * @param iconFileName
	 *            String representing the file name of the icon (*.jpg or *.gif)
	 *            to be assiocated with this menu option. If the name is "null"
	 *            then no icon is associates with this menu.
	 * 
	 * @param shortCut
	 *            Keystroke representing the shortcut for this menu item
	 * 
	 * @param enable
	 *            Boolean value that indicates if this menu item should be
	 *            enabled by default. "true" implies menu item is enabled and
	 *            "false" implies the menu item is disabled.
	 * 
	 * @param auxFlag
	 *            A boolean value that is used to select or unselect a
	 *            CHECK_BOX_ITEM. Refer to the documentation of CheckBoxMenuItem
	 *            in Java SDK API for further information
	 * 
	 * @return The newly created JMenuItem or NULL on error!
	 */
	public static JMenuItem createMenuItem(MenuItemKind itemKind, String itemTitle,
			String command, ActionListener al, String iconFileName,
			KeyStroke shortCut, boolean enable, boolean auxFlag) {

		JMenuItem menuItem = null;
		Icon icon = null;
		int mnemonic = -1;

		// / First we check if the itemTitle has an "_" and make a note
		// / of the mnemonic character to be used and rip it out of
		// / the menu title
		int underScorePosition = itemTitle.indexOf("_");
		if (underScorePosition != -1) {
			// There was an underscore in the string title!
			// let's make a note of the character after it and use the
			// character to make the accelerated key stroke!

			// First remove the "_" from the string! We do this by
			// concatenating two substrings around the "_" together.
			String newTitle = null;
			newTitle = itemTitle.substring(0, underScorePosition)
			+ itemTitle.substring(underScorePosition + 1, itemTitle
					.length());
			itemTitle = newTitle;

			// Now, the actual key to be used for acceleration will be
			// at the index indicated by newTitle! We use it to
			// create the new keystroke.
			String accKey = newTitle.substring(underScorePosition,
					underScorePosition + 1);
			accKey.toUpperCase();
			mnemonic = KeyStroke.getKeyStroke(accKey).getKeyCode();
		}

		// Next check for any icon specifications associated with this
		// menu item. If so we create a new icon here and hook the icon
		// to the menu item in the constructor of the menu item!
		if (iconFileName != null) {
			// okay we have an icon file name. Using the file name we
			// will create an icon!
			icon = getIcon(iconFileName);
		}

		switch (itemKind) {
		case CHECK_BOX_ITEM:
			menuItem = new JCheckBoxMenuItem(itemTitle, icon, auxFlag);
			break;

		case RADIO_BUTTON_ITEM:
			menuItem = new JRadioButtonMenuItem(itemTitle, icon, auxFlag);
			break;

		case SUB_MENU_ITEM:
			menuItem = new JMenu(itemTitle);
			menuItem.setIcon(icon);
			break;
			
		case MENU_ITEM:
		default:
			menuItem = new JMenuItem(itemTitle, icon);
		}

		menuItem.setEnabled(enable);

		if (command != null) {
			menuItem.setActionCommand(command);
			menuItem.addActionListener(al);
			menuItem.setAccelerator(shortCut);
			menuItem.setMnemonic(mnemonic);
		}

		return menuItem;
	}

	/**
	 * Helper method to break messages (using HTML &gt;br/&lt;) at word
	 * boundaries so that they do not exceed a given length.
	 * 
	 * This method is a helper method that is used to break a given 
	 * string into a multi-line HMTL text. If the source string is
	 * already in HTML then your mileage with this method will vary.
	 * Multiple lines of HTML are generated using &gt;br/&lt; HTML
	 * directive at word boundaries. This method does not add any
	 * other HTML decorations or tags to the text.
	 * 
	 * @param srcStr The source string to be converted to a multi-line
	 * HTML text. The source string is directly returned by this method
	 * if: (i) source string is null; (ii) it is an empty string; or
	 * (iii) its length is already less than WrapLen.
	 *  
	 * @param wrapLen The maximum acceptable length (in number of
	 * characters) for each line in the returned HTML text. This 
	 * value must be at least 1.
	 * 
	 * @return The source string converted to a multi-line HTML text
	 * by inserting &gt;br/&lt; at appropriate word boundaries such
	 * that each line of HTML text does not exceed WrapLen characters.
	 */
	public static String wrapStringToHTML(String srcStr, int wrapLen) {
		assert(wrapLen > 0);
		// Check to ensure we have something to work with.
		if ((srcStr == null) || (srcStr.length() < wrapLen)) {
			// Nothing further to be done.
			return srcStr;
		}
		// Fold the subTitle message so that it does not exceed given 
		// number of characters and fold at word boundaries.
		StringBuilder wrappedTitle = new StringBuilder(srcStr.length() + 128);
		while (srcStr.length() > 0) {
			if (wrappedTitle.length() > 0) {
				wrappedTitle.append("<br/>");
			}
			// See if subTitle is longer than wrapLen chars.
			String words = srcStr.substring(0, Math.min(wrapLen, srcStr.length()));
			// Handle the line suitably
			if (words.length() >= wrapLen) {
				// Go back to previous word from the end but not to
				// beginning of string (otherwise we will strip off nothing).
				int lastSpc = Math.max(words.lastIndexOf(' '), 1);
				words = words.substring(0, lastSpc);
				// Save the rest of the line
				srcStr = srcStr.substring(lastSpc + 1);
			} else {
				// Used up all the words
				srcStr = "";
			}
			// Tag on the words we extracted out
			wrappedTitle.append(words);
		}
		return wrappedTitle.toString();
	}
	
	/**
	 * This method is a convenience method that provides a simple mechanism for
	 * creating a JMenuItem, given the parameters.
	 * 
	 * <p><b>Note</b>: For created sub-menu's the return value can be 
	 * safely type cast to a JMenu object.</p>
	 * 
	 * @param itemKind
	 *            This is a enumeration value which can be MENU_ITEM,
	 *            CHECK_BOX_ITEM, or RADIO_BUTTON_ITEM
	 * 
	 * @param itemTitle
	 *            String representing the menu option. Note: The string may
	 *            contain a "_" (underscore) character preceding the letter
	 *            that is supposed to be used as an accelerator for this item.
	 *            However, if an explicit key stroke is specified (parameter 3
	 *            to this method call) then this underscore will be ignored. In
	 *            either case the "_" will be stripped out of the option! The
	 *            string cannot be "null"
	 * 
	 * @param subTitle A sub title string to be used for this menu item. If this
	 * string is null then this method calls the overloaded method that does not
	 * take a subTitle as parameter.
	 * 
	 * @param command
	 *            The string representing the action command that will be
	 *            generated when this menu option is selected!
	 * @param al
	 *            The ActionListener to register with the Menu Item.
	 * 
	 * @param iconFileName
	 *            String representing the file name of the icon (*.jpg or *.gif)
	 *            to be associated with this menu option. If the name is "null"
	 *            then no icon is associates with this menu.
	 * 
	 * @param shortCut
	 *            Keystroke representing the shortcut for this menu item
	 * 
	 * @param enable
	 *            Boolean value that indicates if this menu item should be
	 *            enabled by default. "true" implies menu item is enabled and
	 *            "false" implies the menu item is disabled.
	 * 
	 * @param auxFlag
	 *            A boolean value that is used to select or unselect a
	 *            CHECK_BOX_ITEM. Refer to the documentation of CheckBoxMenuItem
	 *            in Java SDK API for further information
	 * 
	 * @return The newly created JMenuItem or NULL on error!
	 */
	public static JMenuItem createMenuItem(MenuItemKind itemKind, String itemTitle,
			String subTitle, String command, ActionListener al, 
			String iconFileName, KeyStroke shortCut, boolean enable, 
			boolean auxFlag) {
		// Check to ensure we have a valid subTitle to work with
		if (subTitle == null) {
			// Use the overloaded method to create default menu item.
			return createMenuItem(itemKind, itemTitle, command, al, iconFileName, 
					shortCut, enable, auxFlag);
		}
		// Create the menu item with all the given data.
		JMenuItem menuItem = null;
		Icon icon = null;
		int mnemonic = -1;

		// First we check if the itemTitle has an "_" and make a note
		// of the mnemonic character to be used and rip it out of
		// the menu title
		int underScorePosition = itemTitle.indexOf("_");
		if (underScorePosition != -1) {
			// There was an underscore in the string title!
			// let's make a note of the character after it and use the
			// character to make the accelerated key stroke!

			// First remove the "_" from the string! We do this by
			// concatenating two substrings around the "_" together.
			String newTitle = null;
			newTitle = itemTitle.substring(0, underScorePosition)
			+ itemTitle.substring(underScorePosition + 1, itemTitle
					.length());
			itemTitle = newTitle;

			// Now, the actual key to be used for acceleration will be
			// at the index indicated by newTitle! We use it to
			// create the new keystroke.
			String accKey = newTitle.substring(underScorePosition,
					underScorePosition + 1);
			accKey.toUpperCase();
			mnemonic = KeyStroke.getKeyStroke(accKey).getKeyCode();
		}

		// Next check for any icon specifications associated with this
		// menu item. If so we create a new icon here and hook the icon
		// to the menu item in the constructor of the menu item!
		if (iconFileName != null) {
			// okay we have an icon file name. Using the file name we
			// will create an icon!
			icon = getIcon(iconFileName);
		}
		// Fold the subTitle message so that it does not exceed 45 
		// characters and fold at word boundaries.
		final String wrappedTitle = wrapStringToHTML(subTitle, 45);
		// Construct a full HTML item text.
		final String fullItem = "<html><b>" + itemTitle + "</b><br>" +
		"<font size=\"-2\">" + wrappedTitle + "</font></html>";
		// Create suitable menu item.
		switch (itemKind) {
		case CHECK_BOX_ITEM:
			menuItem = new JCheckBoxMenuItem(fullItem, icon, auxFlag);
			break;

		case RADIO_BUTTON_ITEM:
			menuItem = new JRadioButtonMenuItem(fullItem, icon, auxFlag);
			break;

		case SUB_MENU_ITEM:
			menuItem = new JMenu(fullItem);
			menuItem.setIcon(icon);
			break;
		case MENU_ITEM:
		default:
			menuItem = new JMenuItem(fullItem, icon);
		}
		// Set gap to make things better
		if (icon != null) {
			menuItem.setIconTextGap(6);
		}
		// Setup the command string
		if (command != null) {
			menuItem.setActionCommand(command);
			menuItem.addActionListener(al);
			menuItem.setAccelerator(shortCut);
			menuItem.setMnemonic(mnemonic);
		}
		// Disable the menu item if needed.
		if (!enable) {
			setEnabled(menuItem, enable);
		}
		// Let the caller have it.
		return menuItem;
	}

	/**
	 * Utility method to update menu item text that uses HTML.
	 * 
	 * For some reason when menu items use HTML for formatting the 
	 * enabling/disabling of the text does not operate correctly 
	 * because HTML ignores default color settings. This method
	 * introduces colors into the menu items to make them appear
	 * enabled/disabled.
	 * 
	 * <p><b>Note:</b>Calling this method with menu items whose text is not
	 * HTML causes no side effects. So it is safe to call this method
	 * with any menu item.</p>
	 * 
	 * @param menuItem The menu item that must be enabled or disabled.
	 * @param enabled If this parameter is true, then the menu item is
	 * enabled. Otherwise it is disabled.
	 */
	public static void setEnabled(JMenuItem menuItem, boolean enabled) {
		menuItem.setEnabled(enabled);
		String text = menuItem.getText();
		// Reset text item.
		menuItem.setText(enableHTMLText(text, enabled));
	}
	
	/**
	 * Utility method to modify HTML text to look enabled or disabled.
	 * 
	 * For some reason when labels or menu items use HTML for 
	 * formatting the enabling/disabling of the text does not 
	 * operate correctly because HTML ignores default color settings. 
	 * This method introduces colors into the HTML text to make them 
	 * appear enabled/disabled.
	 * 
	 * <p><b>Note:</b>Calling this method with text that is not HTML causes no 
	 * side effects. So it is safe to call this method with any text.</p>
	 *  
	 * @param text The message whose text is to be modified to make
	 * it appear enabled to disabled.
	 * 
	 * @param enabled If this parameter is true, then the text is
	 * updated to appear enabled. Otherwise it is set to appear disabled.
	 * @return The modified string. If the original message wasn't HTML
	 * then this method simply returns the original text.
	 */
	public static String enableHTMLText(String text, boolean enabled) {
		if (!text.startsWith("<html>")) {
			// Plain menu item. nothing to do
			return text;
		}
		// Add color to html depending on enabled status. First strip
		// out leading and trailing <html> and </html>
		text = text.substring(6, text.length() - 7);
		// Remove any font color specifications if present
		if (text.startsWith("<font color")) {
			text = text.substring(text.indexOf('>') + 1, text.length() - 7);
		}
		// Redo HTML text based on status.
		if (!enabled) {
			text = "<font color='gray'>" + text + "</font>";
		}
		text = "<html>" + text + "</html>";
		// Return the modified text.
		return text;
	}
	
	/**
	 * This method is a utility method to create generic buttons.
	 * 
	 * @param iconFileName
	 *            A string representing the icon file (gif/jpeg) that is to be
	 *            associated with this button
	 * @param title
	 *            A string representing the title for this button
	 * @param command
	 *            The command string to be generated when this button is clicked
	 * @param al
	 *            The action listener for this button. When the button is
	 *            clicked the action listener gets triggered by the java system
	 *            where the processing for the button click is performed
	 * @param toolTip
	 *            The string to be associated with the tool tip associated with
	 *            this button
	 * @param enable
	 *            A boolean that indicates if the button must be enabled by
	 *            default ("true" => button enabled. "false" => button disabled)
	 * @return A button that meets the specifications
	 */
	public static JButton createButton(String iconFileName, String title,
			String command, ActionListener al, String toolTip, boolean enable) {
		JButton button = null;

		if (iconFileName != null) {
			button = new JButton(title, Utilities.getIcon(iconFileName));
		} else {
			button = new JButton(title);
		}

		button.setToolTipText(toolTip);
		button.getAccessibleContext().setAccessibleName(title);
		button.setActionCommand(command);
		button.addActionListener(al);
		button.setEnabled(enable);
		button.setRequestFocusEnabled(false);
		// button.setMargin(new Insets(1, 1, 1, 1));
		// button.setBorder(new BevelBorder(BevelBorder.RAISED));
		return button;
	}

	/**
	 * This method is a utility method to create tool bar buttons.
	 * This method is synonymous to calling createButton followed
	 * by makeToolBarButton.
	 * 
	 * @param iconFileName
	 *            A string representing the icon file (gif/jpeg) that is to be
	 *            associated with this button
	 * @param title
	 *            A string representing the title for this button
	 * @param command
	 *            The command string to be generated when this button is clicked
	 * @param al
	 *            The action listener for this button. When the button is
	 *            clicked the action listener gets triggered by the java system
	 *            where the processing for the button click is performed
	 * @param toolTip
	 *            The string to be associated with the tool tip associated with
	 *            this button
	 * @param enable
	 *            A boolean that indicates if the button must be enabled by
	 *            default ("true" => button enabled. "false" => button disabled)
	 * @return A button that meets the specifications
	 */
	public static JButton createToolButton(String iconFileName, String title,
			String command, ActionListener al, String toolTip, boolean enable) {
		// Create button using helper.
		JButton button = createButton(iconFileName, title, command, 
				al, toolTip, enable);
		// Make button a tool bar button.
		makeToolBarButton(button, true);
		// Let the caller have it back.
		return button;
	}

	/**
	 * Method to add mouse adapter to make button like a tool bar button.
	 * 
	 * This method can be used to add a custom mouse adapter to a given
	 * button to make it behave like a tool bar button. The mouse adapter
	 * displays the background and border when a mouse enters the button.
	 * Similarly the mouse adapter clears out the background color and border
	 * when the mouse exits.
	 *  
	 * @param button The button whose properties have to be modified.
	 * 
	 * @param fillContentArea If this parameter is true, then filling 
	 * of the content area is enabled when the mouse rolls over a 
	 * button.
	 */
	public static void makeToolBarButton(AbstractButton button, 
			final boolean fillContentArea) {
		// Turn off the border and fill by default if button is not toggled
		if (!button.isSelected()) {
			button.setBorderPainted(false);
			button.setContentAreaFilled(false);
		}
		// Add a mouse adapter to the button.
		button.addMouseListener(new MouseAdapter() {
			/**
			 * This method overrides the default implementation to 
			 * show border when mouse pointer moves onto the area of
			 * the button.
			 * 
			 * @param e The mouse event from where the actual button
			 * component is obtained.
			 */
			public void mouseEntered(MouseEvent e) {
				Component component = e.getComponent();
				if (component instanceof AbstractButton) {
					AbstractButton button = (AbstractButton) component;
					button.setBorderPainted(true);
					button.setContentAreaFilled(fillContentArea);
				}
			}

			/**
			 * This method overrides the default implementation to 
			 * hide the border when mouse pointer moves onto the area of
			 * the button.
			 * 
			 * @param e The mouse event from where the actual button
			 * component is obtained.
			 */
			public void mouseExited(MouseEvent e) {
				Component component = e.getComponent();
				if (component instanceof AbstractButton) {
					AbstractButton button = (AbstractButton) component;
					// If toggle buttons are depressed they look better with border
					// and content area filled-in.
					if (!button.isSelected()) {
						button.setBorderPainted(false);
						button.setContentAreaFilled(false);
					}
				}
			}
		});
	}

	/**
	 * This method is a helper method that takes the alignment for a
	 * {@link java.awt.FlowLayout}, a label to identify the components, and the
	 * actual JComponents themselves and places them in a JPanel with a
	 * FlowLayout. The panel is useful for situations that require multiple
	 * panels to create a vertical layout of many components. (Note: FlowLayout
	 * performs only a horizontal layout)
	 * 
	 * @param alignment
	 *            The Alignment, one of {@link java.awt.FlowLayout#LEFT},
	 *            {@link java.awt.FlowLayout#RIGHT},
	 *            {@link java.awt.FlowLayout#CENTER},
	 *            {@link java.awt.FlowLayout#LEADING} or
	 *            {@link java.awt.FlowLayout#TRAILING}
	 * @param label
	 *            The label to put before the components (gets converted into a
	 *            JLabel, set to null for none)
	 * @param components
	 *            One or more comma seperated components to place in the JPanel.
	 * @return A new JPanel containing all of the specified components
	 */
	public static JPanel createSubPanel(int alignment, String label,
			JComponent... components) {
		JPanel container = new JPanel(new FlowLayout(alignment));
		if (label != null) {
			container.add(new JLabel(label));
		}
		for (int index = 0; (index < components.length); index++) {
			container.add(components[index]);
		}
		return container;
	}
	
	/**
	 * This method is a helper method to lay out components vertically.
	 * 
	 * This is a helper method that is used to layout components in 
	 * a vertical manner. In addition, this method also creates a
	 * label at the top of the components (if a non-null string 
	 * is specified).
	 * 
	 * @param label
	 *            The label to put before the components (gets converted into a
	 *            JLabel, set to null for none)
	 *            
	 * @param subLabel A sub-label to be placed below the main label. This
	 * label is displayed with a smaller font if it is not null. 
	 * 
	 * @param textBoxHeightDelta If this parameter is non-zero then it adds this
	 * value to the preferred height of any JTextField components in the
	 * component list. This feature is useful to ensure that text fields are
	 * not too small (they happen to be in GTK)
	 * 
	 * @param addEndSpacer If this flag is true, then this method adds a
	 * trailing empty component that will stretch to take up any available
	 * vertical space.
	 * 
	 * @param components
	 *            One or more comma separated components to place in the JPanel.
	 *            
	 * @return A new JPanel containing all of the specified components
	 */
	public static JPanel createLabeledComponents(String label,
			String subLabel,
			int textBoxHeightDelta, boolean addEndSpacer, 
			Component... components) {
		GridBagLayout gbLayout = new GridBagLayout();
		JPanel container       = new JPanel(gbLayout);
		// The same constraints for all components.
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.weightx            = 1.0; // Take all horizontal space
        gbc.fill               = GridBagConstraints.HORIZONTAL;
        gbc.gridwidth          = GridBagConstraints.REMAINDER;
		// Create the label.
		if (label != null) {
			JLabel title = new JLabel(label);
			gbLayout.setConstraints(title, gbc);
			container.add(title);
		}
		if (subLabel != null) {
			JLabel subTitle = new JLabel(subLabel);
			Utilities.adjustFont(subTitle, -2, 8, -1);
			gbLayout.setConstraints(subTitle, gbc);
			container.add(subTitle);
		}
		// Apply the same grid bag constraints to all the components.
		for (int index = 0; (index < components.length); index++) {
			if (components[index] == null) {
				continue;
			}
			if (textBoxHeightDelta > 0) {
				if ((components[index] instanceof JTextField) ||
					(components[index] instanceof JSpinner)) {
					Dimension size = components[index].getPreferredSize();
					size.height   += textBoxHeightDelta;
					components[index].setPreferredSize(size);
				}
			}
			// Setup constraints and add component to container
			gbLayout.setConstraints(components[index], gbc);
			container.add(components[index]);
		}
		// Add a spacer if requested
		if (addEndSpacer) {
			Component spacer = Box.createVerticalGlue();
			// Reset some properties set earlier.
			gbc.weighty    = 1.0;
			gbc.fill       = GridBagConstraints.BOTH;
			gbc.gridheight = GridBagConstraints.REMAINDER;
			// Add the vertical spacer.
			gbLayout.setConstraints(spacer, gbc);
			container.add(spacer);
		}
		// return a new JPanel that contains the components
		return container;
	}

	/**
	 * This method is a quick little helper method to pull image (icons) using
	 * getResource(String).
	 * 
	 * @param path
	 *            The path of the image (relative)
	 * @return An image icon built from that path
	 */
	public static ImageIcon getIcon(String path) {
		ImageIcon icon = null;
		try {
			// We need to check two different paths to handle difference between
			// running via a Jar file and running in development environment (Eclipse).
			// Try the Jar file choice first so that regular runs go faster.
			URL uri = Utilities.class.getResource("/" + path);
			if (uri == null) {
				uri = Utilities.class.getResource(PATH_PREFIX + path);
			}
			// Finally simply try the path
			if (uri == null) {
				File tempFile = new File(path);
				if (tempFile.exists() && tempFile.canRead()) {
					uri = tempFile.toURI().toURL();
				}
			}
			if (uri == null) {
				throw new IOException("The image '" + path + "' was not found."); 
			}
			icon = new ImageIcon(uri);
		} catch (Exception e) {
			// Seems like the image doesn't extist!
			System.out.println("Exception encountered while loading: " + path);
		}
		return icon;
	}

	/**
	 * A helper method to determine the extension of a file.
	 * 
	 * @param f The file whose extension needs to be determined
	 * @return A string containing the extension
	 */
	public static String getExtension(File f) {
		String ext = null;
		String s = f.getName();
		int i = s.lastIndexOf('.');

		if (i > 0 && i < s.length() - 1) {
			ext = s.substring(i + 1).toLowerCase();
		}
		return ext;
	}

	/**
	 * This is a utility method to convert a pair of bytes into a single
	 * integer.
	 * 
	 * @param msb
	 *            The most significant byte of the number.
	 * @param lsb
	 *            The least significant byte of the number.
	 * @return The integer value generated by combining the msb and lsb values.
	 */
	public static int toInteger(byte msb, byte lsb) {
		int msbValue = ((msb < 0) ? (256 + msb) : msb);
		int lsbValue = ((lsb < 0) ? (256 + lsb) : lsb);
		return (msbValue * 256) + lsbValue;
	}

	/**
	 * Returns a byte array containing the two's-complement representation of
	 * the integer.<br>
	 * The byte array will be in big-endian byte-order with a fixes length of 4
	 * (the least significant byte is in the 4th element).<br>
	 * <br>
	 * <b>Example:</b><br>
	 * <code>intToByteArray(258)</code> will return { 0, 0, 1, 2 },<br>
	 * <code>BigInteger.valueOf(258).toByteArray()</code> returns { 1, 2 }.
	 * 
	 * @param integer
	 *            The integer to be converted.
	 * @return The byte array of length 4.
	 */
	public static byte[] integerToByteArray(final int integer) {
		int byteNum = (40 - Integer.numberOfLeadingZeros(integer < 0 ? ~integer
				: integer)) / 8;
		byte[] byteArray = new byte[4];

		for (int n = 0; n < byteNum; n++)
			byteArray[3 - n] = (byte) (integer >>> (n * 8));

		return (byteArray);
	}

	/**
	 * Utility method to change the font size set for a given component.
	 * 
	 * @param c The component whose current font size is to be changed.
	 * 
	 * @param sizeChange The change in font size. This value gets added
	 * to the component's current font size to determine the new font
	 * size. Therefore, positive values will increase the font size while
	 * negative values will decrease the font size.
	 * 
	 * @param minSize The minimum font size below which the font value must
	 * never drop.
	 * 
	 * @param bold If this value is 1, then the font is made bold. If this
	 * value is -1, the font is made normal. If this value is 0 then the
	 * weight of the font is unchanged. 
	 */
	public static void adjustFont(Component c, int sizeChange, int minSize,
			int bold) {
		Font  currFont = c.getFont();
		float newSize  = currFont.getSize() + sizeChange;
		if (newSize < minSize) {
			newSize = minSize;
		}
		Font newFont   = currFont.deriveFont(newSize);
		if (bold != 0) {
			newFont = newFont.deriveFont((bold == 1) ? Font.BOLD : Font.PLAIN);
		}
		c.setFont(newFont);
	}

	/**
	 * Helper method for download a file from a given URL.
	 * 
	 * This method is a helper method that can be used to download a
	 * file from a given URL and save it under a given local file.
	 * This method provides an optional feature of updating progress
	 * via a progress monitor dialog. 
	 * 
	 * @param pm The progress monitor to be updated as data is downloaded.
	 * @param address The URL from where the data is to be downloaded/copied.
	 * @param localFileName The name of the local file (with optional path) where
	 * the data is to be saved.
	 * @exception Exception This method throws various exception on errors.
	 */
	public static void download(ProgressMonitor pm, String address,
			String localFileName) throws Exception {
		OutputStream out = null;
		URLConnection conn = null;
		InputStream in = null;
		try {
			URL url = new URL(address);
			out = new BufferedOutputStream(new FileOutputStream(localFileName));
			conn = url.openConnection();
			in = conn.getInputStream();
			byte[] buffer = new byte[1024];
			int numRead;
			long numWritten = 0;
			while ((numRead = in.read(buffer)) != -1) {
				if ((pm != null) && (pm.isCanceled() == true)) {
					// Delete the file and break out of here
					out.close();
					File del = new File(localFileName);
					del.delete();
					out = null;
					break;
				}
				out.write(buffer, 0, numRead);
				numWritten += numRead;
				if (pm != null) {
					pm.setProgress(numRead);
				}
			}
			// System.out.println(localFileName + "\t" + numWritten);
		} finally {
			if (in != null) {
				in.close();
			}
			if (out != null) {
				out.close();
			}
		}
	}

	/**
	 * This is a helper method that can be used to obtain an input stream for a resource.
	 * 
	 * This method is a convenience method that performs best effort to obtain a valid
	 * input stream. Specifically it handles the following two cases:
	 * 
	 * <ul>
	 * 
	 * <li>First it attempts to load the file from the root path by prepending a "/" to
	 * the resourceName (which is typically a file name). This handles the case where the
	 * resource is to be loaded from a Jar file.</li>
	 * 
	 * <li>If the above step fails, it attempts to obtain the resource by prepending the
	 * PATH_PREFIX constant defined in this class to the resource name.</li>
	 * 
	 * <li>If the above step fails, it assumes the path is an external absolute
	 * path and attempts to open the file.</li>
	 * 
	 * </ul>
	 * 
	 * @param resourceName The relative path/name of the resource to which an input stream
	 * is desired.
	 * 
	 * @return A valid input stream to the requested resource.
	 * 
	 * @throws IOException This method throws an exception if a valid input stream could not
	 * be created to the specified resource.
	 */
	public static InputStream getStream(String resourceName) throws IOException {
		InputStream is = Utilities.class.getResourceAsStream("/" + resourceName);
		if (is == null) {
			is = Utilities.class.getResourceAsStream(PATH_PREFIX + resourceName);
		}
		if (is == null) {
			File tmpFile = new File(resourceName);
			if (tmpFile.exists() && tmpFile.canRead()) {
				is = new FileInputStream(resourceName);
			}
		}
		if (is == null) {
			throw new IOException("Requested resource '" + resourceName + "' was not found.");
		}
		return is;
	}

	/**
	 * This is a helper method to read small text files (like license
	 * information) into a string.
	 * 
	 * This is a helper method that can read a text file into a string.
	 * This method is meant to read and work with small text files,
	 * particularly for display purposes in the GUI. This method should
	 * not be used to load large files as it does use a lot of memory
	 * to hold the text data. This method read files either from a 
	 * Jar file or from the file system.
	 *
	 * @param fileName The name of the file to be loaded. Prefer to use
	 * relative path names so that the file can be loaded immaterial of
	 * whether it is from a jar or from the local file system.
	 * 
	 * @return A string containing the entire contents of the given text
	 * file. All new lines and other special characters (if any) are 
	 * preserved.
	 * 
	 * @throws Exception This method throws an exception if errors occur
	 * during reading the file.
	 */
	public static String readSmallTextFile(String fileName) throws Exception{

		// Obtain input stream to the file. On errors exception
		// will be thrown.

		// We need to check two different paths to handle difference between
		// running via a Jar file and running in development environment (Eclipse).
		// Try the Jar file choice first so that regular runs go faster.
		InputStream is = Utilities.getStream(fileName);
		return readFullStream(is);
	}

	/**
	 * This is a helper method to read small text data from a stream
	 * into a String.
	 * 
	 * This is a helper method that can read textual data into a string.
	 * This method is meant to read and work with small files,
	 * particularly for display purposes in the GUI. This method should
	 * not be used to load large files as it does use a lot of memory
	 * to hold the text data. This method read files either from a 
	 * given stream.
	 *
	 * @param is The input stream from where the data is to be read.
	 * 
	 * @return A string containing the entire contents of the given text
	 * file. All new lines and other special characters (if any) are 
	 * preserved. If the stream did not have any characters then this 
	 * method returns an empty string.
	 * 
	 * @throws Exception This method throws an exception if errors occur
	 * during reading the stream.
	 */
	public static String readFullStream(InputStream is) throws Exception{
		byte buffer[]  = new byte[4096];
		int  bytesRead = -1;
		StringBuilder fileData = new StringBuilder();
		// Read data out of the stream in 4K chunks and append to
		// the file data.
		while ((bytesRead = is.read(buffer)) > 0) {
			// Append bytes read to buffer
			fileData.append(new String(buffer, 0, bytesRead));
		}
		// Close the stream.
		is.close();
		// Return the buffer as a string.
		return fileData.toString();
	}

	/**
	 * Helper method to set the preferred and maximum size of a component.
	 * 
	 * This is a helper method that can be used to adjust the preferred and
	 * maximum sizes of a component.
	 * 
	 * @param c The component whose sizes are to be adjusted.
	 * @param xDelta The change in width with respect to current width.
	 * @param yDelta The change in height with respect to current height.
	 */
	public static void adjustDimension(Component c, int xDelta, int yDelta) {
		Dimension size = c.getPreferredSize();
		size.width    += xDelta;
		size.height   += yDelta;
		c.setPreferredSize(size);
		c.setMaximumSize(size);
	}

	/**
	 * Helper method to recursively enable/disable a component and its children.
	 * 
	 * This is a helper method that can be used to recursively enable or disable
	 * a component and its children. This is a feature that is not provided
	 * by swing. In addition, this method appropriately adds color decorations
	 * to JLabel's that have HTML text in them.
	 * 
	 * @param parent The parent component whose status is to be changed.
	 * @param status If true, the components are enabled otherwise the component is
	 * disabled.
	 */
	public static void setEnabled(Container parent, boolean status) {
		setEnabled(parent, status, null);
	}

	/**
	 * Helper method to recursively enable/disable a component and its children.
	 * 
	 * This is a helper method that can be used to recursively enable or disable
	 * a component and its children. This is a feature that is not provided
	 * by swing. In addition, this method appropriately adds color decorations
	 * to JLabel's that have HTML text in them.
	 * 
	 * @param parent The parent component whose status is to be changed.
	 * @param status If true, the components are enabled otherwise the component is
	 * disabled.
	 * @param ignore This component is ignored and it is not enabled or disabled.
	 * This is helpful in cases where a child component toggles other controls
	 * to be enabled or disabled.
	 */
	public static void setEnabled(Container parent, boolean status, final Component ignore) {
		parent.setEnabled(status);
		final int ChildCount = parent.getComponentCount();
		for(int childIdx = 0; (childIdx < ChildCount); childIdx++) {
			Component child = parent.getComponent(childIdx);
			if (child == ignore) {
				// This is component to ignore
				continue;
			}
			if (child instanceof Container) {
				setEnabled((Container) child, status, ignore);
			} else {
				child.setEnabled(status);
			}
			if (child instanceof JLabel) {
				JLabel label     = (JLabel) child;
				String labelText = label.getText();
				if ((labelText != null) && (labelText.startsWith("<html>"))) {
					// Update the HTML colors to indicate enabled/disabled status
					// as this is not automatically done by Java
					label.setText(enableHTMLText(labelText, status));
				}
			}
		}
	}
	
	/**
	 * Helper method to convert the stack trace in an exception to a 
	 * string.
	 * 
	 * @param e The exception to whose stack trace and information
	 * is desired.
	 * 
	 * @return A string containing the exception's stack trace.
	 */
	public static String toString(Throwable e) {
		StringWriter strStream = new StringWriter();
		PrintWriter  writer    = new PrintWriter(strStream);
		e.printStackTrace(writer);
		return strStream.getBuffer().toString();
	}

	/**
	 * Helper method to convert a given list to a string.
	 * 
	 * This is a convenience method that can be used to convert a given
	 * list to a string. This method converts each element in the list
	 * to a string by calling the {@link #toString()} method on each
	 * entry if the entry is not null. If the list has some null objects
	 * then the literal string "<null>" is used to represent the entry.
	 * There is no way to tell the difference between a null object and
	 * an object whose {@link #toString()} method returned the literal
	 * string "<null>".
	 * 
	 * @param list The list to be converted to a string. This parameter
	 * cannot be null but can be an empty list.
	 * 
	 * @param delimiter The delimiter to be placed between successive
	 * entries in the list. If this delimiter is an empty string then
	 * no delimiter will be present between elements. Note that the 
	 * delimiter cannot be null.
	 * 
	 * @return A string containing the objects in the list separated
	 * by the given delimiter string.
	 */
	public static String toString(final List<? extends Object> list, 
			final String delimiter) {
		StringBuilder sb = new StringBuilder();
		String delim     = ""; // Changed below after first iteration.
		for(int i = 0; (i < list.size()); i++) {
			final Object o   = list.get(i);
			final String val = (o != null) ? o.toString() : "<null>";
			sb.append(delim);
			delim = delimiter; // Change for subsequent entries.
			sb.append(val);
		}
		return sb.toString();
	}
	
	/**
	 * Helper method to create a collapsed pane with message.
	 * 
	 * This is a convenience utility method that creates a panel that
	 * contains two pieces of information. The first one is a "message"
	 * that is placed within a JLabel to be displayed to the user.
	 * This information is constantly visible. The second parameter
	 * "details", is placed within a JTextArea (inside a scroll pane)
	 * that is initially not visible. The text area is made visible
	 * only when the user clicks on a "Details" button that is created
	 * by this method. The complete JPanel can then be placed within
	 * other dialogs (such as JOptionPane.showMessageDialog()) to
	 * provide additional information to the user in a form that does
	 * not overwhelm the user with information. 
	 * 
	 * @param message The message that is to be constantly displayed
	 * to the user via a JLabel.
	 * 
	 * @param details The extra information that will be placed within
	 * a JTextArea that is hidden (or shown) depending on the user's
	 * choice (indicated by clicking on the details button)
	 * 
	 * @return This method returns the JPanel containing a collapsed
	 * details box with the details.
	 */
	public static JPanel collapsedMessage(String message, String details) {
		return collapsedMessage(message, details, true);
	}

	/**
	 * Helper method to create a collapsed pane with message.
	 * 
	 * <p>This is a convenience utility method that creates a panel that
	 * contains two pieces of information. The first one is a "message"
	 * that is placed within a JLabel to be displayed to the user.
	 * This information is constantly visible. The second parameter
	 * "details", is placed within a JTextArea (inside a scroll pane)
	 * that is initially not visible. The text area is made visible
	 * only when the user clicks on a "Details" button that is created
	 * by this method. The complete JPanel can then be placed within
	 * other dialogs (such as JOptionPane.showMessageDialog()) to
	 * provide additional information to the user in a form that does
	 * not overwhelm the user with information.</p> 
	 * 
	 * <p><b>Note:</b>Use this method sparingly and only in circumstances
	 * in which you are absolutely sure that details are to be shown or
	 * hidden. When in doubt, prefer to use the overloaded
	 * {{@link #collapsedMessage(String, String)} method instead.</p>
	 * 
	 * @param message The message that is to be constantly displayed
	 * to the user via a JLabel.
	 * 
	 * @param details The extra information that will be placed within
	 * a JTextArea that is hidden (or shown) depending on the user's
	 * choice (indicated by clicking on the details button)
	 * 
	 * @param showDetails If this flag is true then the details are
	 * visible by default. If this flag is false, then the details are
	 * not visible by default.
	 * 
	 * @return This method returns the JPanel containing a collapsed
	 * details box with the details.
	 */
	public static JPanel collapsedMessage(String message, String details, boolean showDetails) {
		JPanel container = new JPanel(new BorderLayout(5, 5));
		JLabel info      = new JLabel(message);
		Dimension maxSize= info.getPreferredSize();
		maxSize.width     = Math.max(550, maxSize.width);
		info.setMinimumSize(maxSize);
		info.setPreferredSize(maxSize);
		container.add(info, BorderLayout.CENTER);
		// Use the preferred message dimension to setup the
		// dimensions of the collapsible panel.
		final JTextArea msg   = new JTextArea(details);
		final JScrollPane jsp = new JScrollPane(msg);
		jsp.setVisible(showDetails);
		// Setup the maximum scroll pane size so that it looks good.
		maxSize        = info.getPreferredSize();
		maxSize.height = 100;
		maxSize.width  = Math.max(550, maxSize.width);
		jsp.setPreferredSize(maxSize);
		jsp.setMaximumSize(maxSize);
		jsp.setMinimumSize(maxSize);
		// The simple details button with an icon.
		final JToggleButton detailBtn = new JToggleButton("Details  ", Utilities.getIcon("images/16x16/MoreDetails.png"));
		detailBtn.setSelectedIcon(Utilities.getIcon("images/16x16/LessDetails.png"));
		detailBtn.setBorder(null);
		detailBtn.setContentAreaFilled(false);
		detailBtn.setSelected(true);
		detailBtn.setFocusable(false);
		detailBtn.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				// Details button was clicked. Show or hide the message depending on status.
				jsp.setVisible(detailBtn.isSelected());
				// Get top-level parent to validate its layout again.
				SwingUtilities.getWindowAncestor(jsp).pack();
			}
		});
		// Put button and a horizontal line in a suitable sizer.
		Box btnBox = Box.createHorizontalBox();
		btnBox.add(detailBtn);
		btnBox.add(Box.createHorizontalGlue());
		// Add a JSeparator adjacent to button to make it pretty.
		JPanel subBox = new JPanel(new BorderLayout(0, 0));
		subBox.add(btnBox, BorderLayout.NORTH);
		subBox.add(new JSeparator(), BorderLayout.CENTER);
		subBox.add(jsp, BorderLayout.SOUTH);
		// Add the subBox to the main container.
		container.add(subBox, BorderLayout.SOUTH);
		return container;
	}
	
	/**
	 * Helper method to determine the default working directory for
	 * PEACE. This method returns a suitable working directory PATH
	 * depending on the operating system.
	 * 
	 * @return The default working directory for PEACE.
	 */
	public static String getDefaultDirectory() {
		String homeDir  = System.getProperty("user.home");
		String peaceDIR = homeDir + File.separator + "PEACE";
		return peaceDIR;
	}

	/**
	 * Helper method to center a given child on a parent.
	 * 
	 * This method is a helper method that can be used to center a 
	 * given child component on a parent component.
	 * 
	 * @param parent The parent component on which the child component
	 * is to be centered.
	 * @param child The child component to be centered.
	 */
	public static void centerPanel(Component parent, Window child) {
		// Get Java to put the top-left corner in center.
		child.setLocationRelativeTo(parent); 
		if (parent != null) {
			// Adjust top a bit more if parent is not null
			Point tlc      = child.getLocation();
			Dimension size = child.getPreferredSize();
			tlc.x         -= (size.width  / 2);
			tlc.y         -= (size.height / 2);
			// Ensure x and y values are not smaller than 0
			tlc.x = Math.max(0, tlc.x);
			tlc.y = Math.max(0, tlc.y);
			// Set the location back.
			child.setLocation(tlc);
		}
	}
	
	/**
	 * Helper method to trim a string to a given length.
	 * 
	 * @param src The source string (typically a description) to be
	 * trimmed to a given length.
	 * 
	 * @param maxLen The maximum length of the return string. If the
	 * src description is longer than this value, then trailing
	 * characters after this length are dropped.
	 * 
	 * @return The first part of src string such that the length of
	 * the returned string is less than maxLen. 
	 */
	public static String trim(String src, int maxLen) {
		if (src.length() > maxLen) {
			return src.substring(0, maxLen - 2) + "...";
		} 
		return src;
	}
	
	/**
	 * Helper method to try and generate a unique file name.
	 * 
	 * This method is a convenience method that can be used to generate
	 * a unique file name such that a file with the same name does not 
	 * exist in the given directory. This method essentially generates
	 * file names using monotonically increasing sequence of numbers 
	 * (starting with 1) and checks to see if a file with the given name
	 * exists. For example, given a directory <code>/home/raodm</code>
	 * and fileNamePrefix <code>test</code> and extension <code>.fasta</code>
	 * this method generates temporary file name in the form 
	 * <code>/home/raodm/test_1.fasta</code> and checks to see if the 
	 * file with this name exists. If it does it tries the file with
	 * the next higher number until a non-existent file is encountered.
	 *   
	 * @param directory The directory in which a unique file name is
	 * desired. This directory must exist. If the directory does not
	 * exist then the effectiveness of this method is unspecified and
	 * will be platform dependent.
	 * 
	 * @param fileNamePrefix The prefix to be used to generate the file name.
	 * This value cannot be null.
	 * 
	 * @param extension The file extension to be used. This value cannot be
	 * null.
	 * 
	 * @return The full path to the unique file within the given directory.
	 */
	public static String getUniqueFileName(final String directory, 
			final String fileNamePrefix, final String extension) {
		String uniqueFileName = null;
		int trial = 1;
		while (uniqueFileName == null) {
			final String trialFileName = directory + 
				File.separator + fileNamePrefix + "_" + trial + extension;
			final File trialFile = new File(trialFileName);
			if (!trialFile.exists()) {
				uniqueFileName = trialFile.getAbsolutePath();
			}
			// Try next file name
			trial++;
		}
		return uniqueFileName;
	}
	
	/**
	 * Convenience method to extract a word that occurs after a given phrase in
	 * a string.
	 * 
	 * This is a helper method can be used to extract a word that occurs
	 * after a given phrase in a string. For example given
	 * a string <code>"Generated 201035 Base Pairs\nInsertions: 6009"</code>
	 * and the phrase <code>"Insertions: "</code>, this method returns 
	 * the string <code>"6009"</code>. 
	 * 
	 * @param source The source string in which to search for the given
	 * phrase and return the word that occurs after the phrase.
	 * 
	 * @param phrase The phrase to search for. The word to be extracted is 
	 * assumed to start right after this phrase. Initial white spaces right
	 * after the phrase are skipped over. The word ends when the
	 * next white space (a space, tab, or newline) is encountered in the
	 * string or the end of string is hit.
	 * 
	 * @return If the word was found an extracted successfully, then this
	 * method returns the word. Otherwise this method returns null. If the
	 * phrase was at end-of-string then this method returns an empty string.
	 */
	public static String getWordAfter(final String source, final String phrase) {
		String retVal = null;
		final int phrasePos = source.indexOf(phrase); // find phrase.
		if (phrasePos != -1) {
			// Phrase found. Extract next word after the phrase.
			int wordPos = phrasePos + phrase.length();
			// Skip over white spaces until a non-white-space character
			while ((wordPos < source.length()) && Character.isWhitespace(source.charAt(wordPos))) {
				wordPos++;
			}
			// Gather characters after the phrase until a white space
			// is encountered.
			StringBuilder word = new StringBuilder(32);
			while ((wordPos < source.length()) && !Character.isWhitespace(source.charAt(wordPos))) {
				word.append(source.charAt(wordPos++));
			} 
			// Extract the word to be returned to the user
			retVal = word.toString();
		}
		return retVal;
	}
}
