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
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * <p>This class is a utility class that provides a general purpose "Find" dialog that
 * can be presented to the user to search for information. This dialog provides 
 * a couple of additional features that facilitates repeated use of the dialog to
 * perform searches. One of the objectives of this class is to enable reuse of the
 * find dialog across multiple tabs. Consequently, this class has been designed to be
 * used as a singleton class.</p>
 * 
 * <p>Typical use of this dialog class is performed in the following manner:
 * 
 * <ul>
 * 
 * <li>The reference to the singleton instance of this dialog must be obtained via
 * the getDialog() method in this class.</li>
 * 
 * <li>Next, the window that is interested in receiving and processing search events must
 * first set a FindListener that receives find requests via the setFindListener method.</li>
 *  
 * <li>Next, the find dialog is displayed by calling the static interface method 
 * showDialog()</li>
 * 
 * <li>Whenever the user clicks the "Find" button in the dialog, a FindEevent is dispatched
 * to the FindListener (set via call to setFindListener method) to perform the find. The FindEvent
 * contains additional information indicating options and parameters set by the user for the Find
 * operation.</li>
 * 
 * </ul>
 * 
 * </p>
 */
public class FindDialog extends JDialog implements ActionListener {
	/**
	 * <p>The primary API method to create create a singleton instance of the FindDialog.
	 * This method is called only once from the MainFrame class when it is instantiated.
	 * This method checks to see if a singleton instance already exists. If not, it 
	 * instantiates the process-wide unique instance of the find dialog for further use.</p>
	 *  
	 * @param parent The parent frame that logically owns this dialog and whose
	 * child views will utilize this dialog.
	 */
	public static synchronized void createFindDialog(JFrame parent) {
		if (singletonFindDialog == null) {
			singletonFindDialog = new FindDialog(parent);
		}
	}
	
	/**
	 * Obtain the reference to the process-wide unique (and shared) instance of
	 * the FindDialog.
	 * 
	 * @return A reference to the process-wide unique (and shared) instance of
	 * the FindDialog. If a dialog has not yet been created, then this method
	 * returns null.
	 */
	public static FindDialog getDialog() {
		return singletonFindDialog;
	}
	
	/**
	 * Set the listener to receive notifications when the "Find" button is clicked.
	 * 
	 * @param listener The listener to be notified when the user wishes to
	 * conduct another search.
	 */
	public void setFindListener(FindListener listener) {
		this.listener = listener;
	}
	
	/**
	 * <p>The <b>only</b> constructor that is used to create an instance of the FindDialog.
	 * This constructor is called only once (as this dialog is designed as a singleton)
	 * from the {@link #createFindDialog(JFrame)} API method (which is called from the
	 * MainFrame when it is instantiated).</p>
	 * 
	 * <p>The constructor essentially creates and lays out the various components constituting
	 * the find dialog.</p>
	 * 
	 * @param parent The parent frame that logically owns this dialog and whose
	 * child views will utilize this dialog.
	 */
	private FindDialog(JFrame parent) {
		super(parent, "Find Dialog", false);
		// Setup basic properties
		setLayout(new BorderLayout(5, 10));
		setIconImage(Utilities.getIcon("images/16x16/Find.png").getImage());
		
		JPanel rootPanel = new JPanel(new BorderLayout(0, 5));
		rootPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
		
		// Create and add the top-panel component(s)
		rootPanel.add(createTopPanel(), BorderLayout.NORTH);
		// Create and add the middle options panel.
		rootPanel.add(createOptionsPanel(), BorderLayout.CENTER);
		// Finally create the bottom panel with the two buttons.
		rootPanel.add(createButtonPanel(), BorderLayout.SOUTH);
		// Add root panel to the main window
		add(rootPanel, BorderLayout.CENTER);
		// Pack the layout preparing it for display
		pack();
	}

	/**
	 * Helper method to create the bottom button panel in the find dialog. This 
	 * method is called only once from the constructor. This method was introduced
	 * to streamline the constructor and keep the code clutter to a minimum.
	 * 
	 * @return This method returns a JPanel containing the "Find" and "Close" 
	 * buttons (suitably organized).
	 */
	JPanel createButtonPanel() {
		// Create buttons and setup their properties.
		JButton find = Utilities.createButton("images/16x16/Find.png", " Find ", "find", 
				this, "Do a find operation", true);
		find.setDefaultCapable(true);
		JButton close = Utilities.createButton("images/16x16/CloseWindow.png", " Close ", "close", 
				this, "Close this dialog box", true);
		close.setMnemonic(KeyEvent.VK_ESCAPE);
		// Wrap the buttons in a suitable panel.
		JPanel btnPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 20, 5));
		btnPanel.add(find);
		btnPanel.add(close);
		
		return btnPanel;
	}
	
	/**
	 * Helper method to create the middle options-panel in the find dialog. This 
	 * method is called only once from the constructor. This method was introduced
	 * to streamline the constructor and keep the code clutter to a minimum.
	 * 
	 * @return This method returns a JPanel containing the search options (like:
	 * forward/backwards, wrap search, case sensitive etc.)
	 */
	private JPanel createOptionsPanel() {
		// First create the directional radio button pair.
		JRadioButton forward  = new JRadioButton("Forward search", true);
		forward.setActionCommand("forward");
		JRadioButton backward = new JRadioButton("Backward search", false);
		backward.setActionCommand("backward");
		// Pack them into a ButtonGroup to make them mutually exclusive.
		forwardBackward = new ButtonGroup();
		forwardBackward.add(forward);
		forwardBackward.add(backward);
		// Pack the direction options into a nice bordered panel.
		JPanel dirPanel = new JPanel(new GridLayout(3, 1));
		dirPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Direction"));
		dirPanel.add(forward);
		dirPanel.add(backward);
		
		// Next create and layout the option check boxes.
		caseSensitive = new JCheckBox("Case sensitive", false);
		wrapSearch    = new JCheckBox("Wrap search", false);
		regExp        = new JCheckBox("Regular expression", false);
		// Add the three to a nice bordered panel.
		JPanel optPanel = new JPanel(new GridLayout(3, 1));
		optPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Options"));
		optPanel.add(caseSensitive);
		optPanel.add(wrapSearch);
		optPanel.add(regExp);
		
		// Create the warning label for display. By default the label does not have any content
		// but just occupies space for future use.
		warningMessage = new JLabel("", FindDialog.BlankIcon, JLabel.CENTER);
		
		// Create a bigger panel to wrap the dirPanel and optPanel.
		JPanel wrapper = new JPanel(new BorderLayout(5, 5));
		wrapper.add(dirPanel, BorderLayout.WEST);
		wrapper.add(optPanel, BorderLayout.EAST);
		wrapper.add(warningMessage, BorderLayout.SOUTH);
		
		return wrapper;
	}
	
	/**
	 * Helper method to create the top-panel in the find dialog. This method is called
	 * only once from the constructor. This method was introduced to streamline the
	 * constructor and keep the code clutter to a minimum.
	 * 
	 * @return This method returns a JPanel containing the search-term entry box along
	 * with some text and icons.
	 */
	private JPanel createTopPanel() {
		searchTerm = new JComboBox();
		searchTerm.setEditable(true);
		searchTerm.setActionCommand("searchTerm");
		searchTerm.addActionListener(this);
		
		JPanel findPanel = new JPanel(new BorderLayout(0, 2));
		findPanel.add(new JLabel("Find (string or regular expression):"), BorderLayout.NORTH);
		findPanel.add(searchTerm, BorderLayout.SOUTH);
		
		// Create a image to make things look nice
		JLabel img = new JLabel(Utilities.getIcon("images/32x32/Find.png"));
		JPanel topPanel = new JPanel(new BorderLayout(10, 5));
		topPanel.add(img, BorderLayout.WEST);
		topPanel.add(findPanel, BorderLayout.CENTER);
		return topPanel;
	}
	
	/**
	 * Helper method to validate if the inputs are valid. This method is invoked from
	 * the {@link #actionPerformed(ActionEvent)} method. This method performs the
	 * following tasks:
	 * 
	 * <ol>
	 * 
	 * <li>If the current search term is empty, then it sets up a warning and
	 * returns with false.</li>
	 * 
	 * <li>If the search term is not flagged to be used as a regular expression it
	 * returns true.</li>
	 * 
	 * <li>It validates to ensure that the search term is a valid regular expression.
	 * If not, it sets a warning message and returns with false. If the regular 
	 * expression is valid, then it returns true.
	 * 
	 * </ol>
	 * 
	 * @return This method returns true if the inputs are valid. Otherwise it returns
	 * false.
	 */
	private boolean validateInputs() {
		Object comboBoxInput = searchTerm.getEditor().getItem();
		// Get the search string entered by the user.
		String searchStr  = (comboBoxInput != null) ? comboBoxInput.toString() : null;
		String warningMsg = null;  // Set to a valid string upon warnings
		
		if ((searchStr == null) || (searchStr.length() < 1)) {
			warningMsg = "A search term has not been provided";
		} else if (regExp.isSelected()) {
			// Check to ensure that the specified search term is a valid
			// regular expression.
			try {
				// Compile regex to ensure it is valid.
				regExPattern = Pattern.compile(searchStr);
			} catch (PatternSyntaxException pse) {
				// The regular expression is invalid.
				warningMsg = "Invalid regular expression";
			}
		} else {
			// Simple string search.
			regExPattern = null;
		}
		
		// When control drops here, we have a valid search term if warning is null.
		if (warningMsg != null) {
			// Can't proceed with find operation due to error.
			warningMessage.setText(warningMsg);
			warningMessage.setIcon(FindDialog.ErrorIcon);
			Toolkit.getDefaultToolkit().beep();
			return false;
		}
		// The search term is good. First add it to choices in our combo box
		// for future reference (only if it is not already there).
		for(int index = 0; (index < searchTerm.getItemCount()); index++) {
			if (searchTerm.getItemAt(index).equals(searchStr)) {
				// Duplicate entry found. Remove it as a new entry will
				// be added at the top of the list.
				searchTerm.removeItemAt(index);
				break;
			}
		}
		searchTerm.insertItemAt(searchStr, 0);
		searchTerm.setSelectedIndex(0);
		// Clear out any search text (if any).
		warningMessage.setText("");
		warningMessage.setIcon(FindDialog.BlankIcon);
		return true;
	}
	
	/**
	 * Helper method to create and dispatch a find event when the "Find" button
	 * is clicked. This method is invoked from the {@link #actionPerformed(ActionEvent)}
	 * method after the inputs have been validated. This method creates and dispatches
	 * a FindEvent to the currently set find listener (if any).
	 * 
	 * @return This method returns true if a valid listener is present and it reported
	 * that a successful find operation was performed. Otherwise this method returns false.
	 */
	private boolean fireFindEvent() {
		if (listener != null) {
			// Extract the necessary flags and options for convenience.
			final String dirString   = forwardBackward.getSelection().getActionCommand();         
			final boolean searchDir  = dirString.equals("forward");
			final boolean wrapAround = wrapSearch.isSelected();
			final boolean caseSense  = caseSensitive.isSelected();
			// Create a suitable find event based on the flags.
			FindEvent fe;
			if (regExp.isSelected()) {
				fe = new FindEvent(this, regExPattern, searchDir, wrapAround, caseSense);
			} else {
				final String searchStr = searchTerm.getSelectedItem().toString();
				fe = new FindEvent(this, searchStr, searchDir, wrapAround, caseSense);
			}
			assert( fe != null );
			return listener.find(fe);
		}
		// No listener 
		return false;
	}
	
	/**
	 * This method is invoked from various buttons and options present in this dialog
	 * box. This method handles the various call backs generated when buttons are 
	 * clicked or other operations are performed by the user. This method uses the
	 * action commands associated with the various GUI elements in this dialog to
	 * perform the appropriate operations.
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		final String cmd = e.getActionCommand();
		if ("find".equals(cmd) || "searchTerm".equals(cmd)) {
			if (validateInputs()) {
				fireFindEvent();
			}
		} else if ("close".equals(cmd)) {
			this.setVisible(false);
		}
	}
	
	/**
	 * This text field is used to obtain the actual search string
	 * or regular expression from the user. This is created in the
	 * {@link #createTopPanel()} method and used in different spots in
	 * this class.
	 */
	JComboBox searchTerm;
	
	/**
	 * The button group that contains the pair of radio buttons that can be used by
	 * the user to indicate the search should proceed in the forward or backward 
	 * direction. The first entry is for the "Forward" search option while the 
	 * second entry (in the pair) is for the"Backward" search option. These radio
	 * buttons are created by the  {@link #createOptionsPanel()} method and placed 
	 * within this ButtonGroup to make them mutually exclusive.
	 */
	ButtonGroup forwardBackward;
	
	/**
	 * The check box to indicate if the search is case sensitive or case insensitive.
	 * This check box is created by the {@link #createOptionsPanel()} and is used
	 * by various methods in this class.
	 */
	JCheckBox caseSensitive;

	/**
	 * The check box to indicate if the search should logically wrap around.
	 * This check box is created by the {@link #createOptionsPanel()} and is used
	 * by various methods in this class.
	 */
	JCheckBox wrapSearch;

	/**
	 * The check box to indicate if the search should treat the search expression
	 * as a regular expression (or not).  This check box is created by the 
	 * {@link #createOptionsPanel()} and is used by various methods in this class.
	 */
	JCheckBox regExp;
	
	/**
	 * A simple message that is displayed (as needed) just above the buttons.
	 * For example, this label is used to report if a given regular expression
	 * is valid. 
	 */
	JLabel warningMessage;
	
	/**
	 * Reference to the process-wide unique (and shared) singleton instance of
	 * the find dialog that is shared between various views in the GUI subsystem.
	 */
	private static FindDialog singletonFindDialog = null;
	
	/**
	 * The listener that must receive notifications for searching for a specified
	 * string. The listener is set via the #setFindListener method.
	 */
	private FindListener listener;
	
	/**
	 * This instance variable is used to track the current regular expression
	 * pattern being used for performing the search. This instance variable
	 * is set in the {@link #validateInputs()} method and used in the
	 * {@link #fireFindEvent()} method.
	 */
	private Pattern regExPattern;
	
	/**
	 * A simple error icon that is used when displaying error messages in this dialog.
	 */
	private static final ImageIcon ErrorIcon = Utilities.getIcon("images/16x16/Error.png");
	
	/**
	 * A blank icon that is used instead of the error icon when regular messages (or
	 * no message) is displayed in this dialog.
	 */
	private static final ImageIcon BlankIcon = Utilities.getIcon("images/16x16/Blank.png");

	/**
	 * A generated serialization GUID (just to keep the compiler happy) 
	 */
	private static final long serialVersionUID = 1431317283690842814L;

	public static void main(String [] args) {
		JFrame frame = new JFrame();
		frame.setPreferredSize(new Dimension(600,600));
		FindDialog.createFindDialog(frame);
		FindDialog findDiag = FindDialog.getDialog();
		frame.setVisible(true);
		findDiag.setVisible(true);
	}
}