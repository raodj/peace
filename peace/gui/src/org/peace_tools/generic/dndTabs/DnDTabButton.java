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

package org.peace_tools.generic.dndTabs;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.peace_tools.generic.Utilities;

/**
 * This is a simple component that provides a custom close button on 
 * the currently active tab in a DnDTabbedPane.
 * 
 * This is a simple tab component that is set to render the title 
 * for each tab added to a DnDTabbedPane. This component (extends
 * a JPanel) uses a JLabel and JButton to draw the title. The JLabel
 * is used to draw the usual tab title (icon + text) while the button
 * serves as the close button. 
 */
public class DnDTabButton extends JPanel implements ActionListener, ChangeListener {
	/**
	 * A generated serialization UID to keep the complier happy
	 */
	private static final long serialVersionUID = -6226717869585754033L;

	/**
	 * The tab pane which which this tab button is logically
	 * associated. This value is set in the constructor when the
	 * button is created.
	 */
	private final DnDTabbedPane pane;

	/**
	 * The close button that is visible on the currently active
	 * tab. The close button simply has icons that are displayed
	 * at different times.
	 */
	private final JButton     button;

	/**
	 * The constructor.
	 * 
	 * The constructor lays out the panel with a jabel to the left
	 * and the close button to the right. The close button is made
	 * visible if the index (of the currently visible tab) matches
	 * the index where the tab component is going to be placed. 
	 * 
	 * @param pane The tab pane to which the tab is going to be
	 * added.
	 * 
	 * @param index The index where this tab component is going to be set.
	 */
	public DnDTabButton(final DnDTabbedPane pane, int index) {
		// Reset default FlowLayout' gaps
		super(new FlowLayout(FlowLayout.LEFT, 0, 0));
		assert (pane != null);
		this.pane = pane;
		// Get notification when tabs change to show/hide button
		pane.addChangeListener(this);
		setOpaque(false);

		// Make JLabel read titles from JTabbedPane. Here we get the
		// title and icon dynamically so that any changes to this
		// information in the tab pane is suitably reflected.
		JLabel label = new JLabel() {
			private static final long serialVersionUID = 6659971162564603892L;
			@Override
			public Icon getIcon() {
				int idx = pane.indexOfTabComponent(DnDTabButton.this);
				if (idx != -1) {
					return pane.getIconAt(idx);
				}
				return null;
			}
			@Override
			public String getText() {
				int idx = pane.indexOfTabComponent(DnDTabButton.this);
				if (idx != -1) {
					return pane.getTitleAt(idx);
				}
				return null;
			}
		};
		add(label);

		// -----------------------------------------------------------
		// Create close button.
		button = new JButton(Utilities.getIcon("images/16x16/CloseTabDefault.png"));
		button.setToolTipText("Close this tab");
		// Make it transparent
		button.setContentAreaFilled(false);
		// No need to be focusable
		button.setFocusable(false);
		// Set small border on top for this button.
		button.setBorder(new EmptyBorder(1, 3, 0, 0));
		button.setBorderPainted(false);
		// Set up roll over effect by changing icons.
		button.setRolloverEnabled(true);
		// Setup icons for regular and roll over state
		button.setRolloverIcon(Utilities.getIcon("images/16x16/CloseTabOver.png"));
		// Close the proper tab by clicking the button
		button.addActionListener(this);

		// Set initial visibility of this button.
		if (pane.getSelectedIndex() == -1) {
			button.setVisible(pane.getSelectedIndex() == index);
		} else {
			// No tab is selected right now.
			button.setVisible(index == 0);
		}
		// Make the button sensitive to roll overs etc.
		add(button);
	}

	/**
	 * Close tab when button is clicked.
	 * 
	 * This method implements the only method in the Action Listener
	 * interface. This method is called when the close button is
	 * clicked. This method closes the corresponding pane in the tab
	 * panel.
	 * 
	 * @param e The action event associated with this event. This
	 * parameter is currently unused.
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		int i = pane.indexOfTabComponent(this);
		if (i != -1) {
			pane.deleteTab(pane.getComponentAt(i));
		}
	}

	/**
	 * Intercept tab changes and show/hide the close button.
	 * 
	 * This method intercepts tab changes on the dnDTabPanel. It
	 * suitably shows or hides the close button.
	 */
	@Override
	public void stateChanged(ChangeEvent e) {
		// Show / hide kill button if this tab is not active.
		int i = pane.indexOfTabComponent(DnDTabButton.this);
		button.setVisible(i == pane.getSelectedIndex());
	}
}
