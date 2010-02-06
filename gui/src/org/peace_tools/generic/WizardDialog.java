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
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;

/**
 * A extended wizard dialog class.
 * 
 * This class extends the features of a standard JDialog class and provides some
 * custom graphics and capabilities to develop a Wizard. A wizard is essentially
 * a series of dialogs that guide the user through a sequence of steps that are
 * needed to accomplish a specific task. The specific features of this class
 * include:
 * 
 * <ul>
 * 
 * <li>Ability to add one or more WizardPanel objects in a specific sequence to
 * the pagePanel of this wizard dialog class. The page panel has a
 * CardLayout manager that displays only one WizardPanel at a time. The
 * WizardPanel objects are displayed in the order in which they are present in
 * the content pane.</li>
 * 
 * <li>A custom title image that is displayed as the header for this wizard
 * dialog.</li>
 * 
 * <li>A custom column image that is shared/reused by all the wizard pages
 * added to this dialog. This image is displayed on the left hand side of the
 * dialog (below the list of steps in the wizard)</li>
 * 
 * <li>A custom display of the sequence of steps in this wizard along with
 * indication on which sequence have already been completed. </li>
 * 
 * </ul>
 */
public class WizardDialog extends JDialog implements ActionListener {
	/**
	 * The constructor.
	 * 
	 * The constructor merely passes the parameters to the base class and
	 * initializes the shared (by all the WizardPage objects added to this
	 * wizard) static bitmap images in this class (if they have not yet been
	 * initialized).
	 * 
	 * @param parent The parent window to which this wizard dialog logically
	 * belongs.
	 */
	protected WizardDialog(Frame parent) {
		super(parent, "", true);
		// Set the preferred layout to be border layout.
		setLayout(new BorderLayout(0, 0));
		// Set up an overall border via a JPanel.
		JPanel contentPane = new JPanel(new BorderLayout(0, 0));
		CompoundBorder cb = new CompoundBorder(new CustomBorder("ddLL"),
				new CustomBorder("DDss"));
		setContentPane(contentPane);
		contentPane.setBorder(cb);
		// Initialize all the default panels.
		pagePanel = new CustomPanel(new CardLayout());
		// Create the title panel with a box layout.
		titlePanel = new CustomPanel(new BorderLayout(0, 0));
		titlePanel.setUseImageSize(false, true);
		// Next create the sequence panel to display sequence overview
		sequencePanel = new CustomPanel();
		sequencePanel.setLayout(new BoxLayout(sequencePanel, BoxLayout.Y_AXIS));
		// Create the four basic buttons to be displayed in the wizard.
		prevButton = Utilities.createButton("images/16x16/WizardPrev.png", 
				" Back ", "prev", this, 
				"Revert to previous page in this wizard", false);
		nextButton = Utilities.createButton("images/16x16/WizardNext.png", 
				" Next ", "next", this, 
				"Proceed to next page in this wizard", false);
		cancelButton = Utilities.createButton("images/16x16/WizardCancel.png", 
				"Cancel", "cancel", this, 
				"Exit from this wizard without completing it", false);
		finishButton = Utilities.createButton("images/16x16/WizardFinish.png", 
				"Finish", "finish", this, 
				"Complete all operations of this wizard", false);
		helpButton = Utilities.createButton("images/16x16/Help.png", 
				" Help ", "help", this, 
				"Show online help page for this wizard", false);
		// Setup the array of pages.
		pages = new ArrayList<WizardPage>();
		// Set up a handler for the close
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				if (cancel()) {
					setVisible(false);
					dispose();
					// Notify that the wizard is done.
					done(false);
				}
			}
		});
	}

	/**
	 * Set an image to be used for the title panel. This image is
	 * displayed in the background of the title panel. The title(s)
	 * are displayed on top of this image.
	 * 
	 * @param fileName The file name of the image to be set as the
	 * background for the title of this wizard.
	 * 
	 * @param bgColor The default background color for the title panel.
	 * If this parameter is null, then the default color is used.
	 */
	public void setTitleBackground(String fileName, Color bgColor) {
		titlePanel.setImage(fileName);
		if (bgColor != null) {
			titlePanel.setBackground(bgColor);
		}
	}
	
	/**
	 * Set an image to be used for the sequence panel. This image is
	 * displayed in the background of the sequence panel that is
	 * displayed to the left of the wizard. The sequence panel displays
	 * an overview of the various sequences involved in this wizard.
	 * 
	 * @param fileName The file name of the image to be set as the
	 * background for the sequence panel.
	 */
	public void setSequenceBackground(String fileName) {
		sequencePanel.setImage(fileName);
	}

	/**
	 * Adds a new wizard page to the list of wizard dialog pages.
	 * 
	 * This method adds the given wizard page to the set of pages
	 * to be displayed by this
	 *  
	 * @param page The wizard page to be added to the list of pages.
	 */
	public void addPage(WizardPage page) {
		pages.add(page);
	}

	/**
	 * Adds a given worker thread to the wizard.
	 * 
	 * This method can be used to add a worker thread used by a 
	 * wizard page. When the thread exits it must remove itself 
	 * via the remoteThread method. The threads are interrupted
	 * on the user decided to cancel out of the wizard.
	 * 
	 * @param t The thread to be added (and deleted when aborting)
	 * to the list of worker threads associated with this wizard.
	 */
	public synchronized void addThread(Thread t) {
		workerThreads.add(t);
	}
	
	/**
	 * Remove a given worker thread from the list of worker threads.
	 * 
	 * This method may be used to remove a worker thread from the
	 * list of worker threads associated with this wizard.
	 * 
	 * @param t
	 */
	public synchronized void removeThread(Thread t) {
		workerThreads.remove(t);
	}
	
	/**
	 * The main method that lays out the components in this wizard
	 * and displays the wizard. 
	 * 
	 * @param helpURL The help URL to be used for this wizard if any.
	 * This parameter can be null.
	 */
	public void showWizard(String helpURL) {
		// Set the help URL to be used
		this.helpURL = helpURL;
		// Create the layout for the wizard
		layoutWizard();
		// Pack the wizard forcing layout to happen
		pack();
		// Center it on the parent.
		setLocationRelativeTo(getOwner()); 
		// Set the current page being displayed
		currentPage = 0;
		// Now the indicate that the first page is to be displayed
		setButtonStatus();
		// Update the sequences
		updateSequence(0);
		// Inform the page about this change.
		pages.get(0).pageChanged(this, -1, 0);
		// Finally show the wizard.
		setVisible(true);
	}
	
	/**
	 * Helper method to update the current sequence information and
	 * highlight the current step in the wizard.
	 * 
	 * @param currPage The current page being displayed.
	 */
	private void updateSequence(int currPage) {
		// First set the labels for the completed sequences with
		// a little check mark.
		for(int pageID = 0; (pageID < pages.size()); pageID++) {
			JLabel label = seqTitles[pageID];
			// Reset the border for the sequence title.
			label.setBorder(new EmptyBorder(7, 7, 7, 25));
			// Make the label to be non-opaque and update color
			label.setOpaque(false);
			label.setForeground(Color.white);
			// Setup suitable icon.
			if (pageID < currPage) {
				label.setIcon(Utilities.getIcon("images/16x16/CheckedBoxWhite.png"));
			} else {
				label.setIcon(Utilities.getIcon("images/16x16/BoxWhite.png"));
			}
		}
		// Reset the information for the current page.
		JLabel label = seqTitles[currPage];
		// Ensure the empty border extends to the edge of the pane.
		int rightMargin = sequencePanel.getWidth() - label.getWidth() + 19;
		CustomBorder lineBorder = new CustomBorder(Color.black,
				Color.black, Color.black, null);
		label.setBorder(new CompoundBorder(lineBorder, new EmptyBorder(7, 7, 7, rightMargin)));
		// Make the label to be opaque so it is grey and make stuff
		// black to improve readability
		label.setOpaque(true);
		label.setForeground(Color.black);
		label.setIcon(Utilities.getIcon("images/16x16/Box.png"));
	}
	
	/** Method to suitably set status of buttons in the wizard.

		This method can be used to enable/disable the buttons
		in this wizard panel. 
	 */
	private void setButtonStatus() {
		prevButton.setEnabled(currentPage > 0);
		nextButton.setEnabled(currentPage < pages.size() - 1);
		nextButton.setVisible(currentPage < pages.size() - 1);
		
		finishButton.setEnabled(currentPage == pages.size() - 1);
		finishButton.setVisible(currentPage == pages.size() - 1);
		// Enable cancel button all the way to the end.
		cancelButton.setEnabled(true);
		// Set the titles for this page.
		WizardPage page = pages.get(currentPage);
		titles[0].setText(page.getTitle());
		titles[1].setText(page.getSubTitle());
		// Force a layout to ensure titles are layed out properly
		// titlePanel.doLayout();
	}

	/** Method to enable or disable back and/or next buttons.

    	This method can be used to enable or disable the back
    	and/or next buttons in a wizard dialog. This method accepts
    	integers with values -1, 0, and 1 that are used to represent
    	the following three states respectively: unchanged, disable,
    	or enable.

    	@param prev This value must be -1, 0, or 1 to ignore, 
    	disable, or enable the back button in the wizard dialog.

    	@param next This value must be -1, 0, or 1 to ignore, 
    	disable, or enable the next button in the wizard dialog.

    	@param cancel This value must be -1, 0, or 1 to ignore,
    	disable, or enable the cancel button in the wizard dialog.
	 */
	public void setButtonStatus(int prev, int next, int cancel) {
		 if (prev != -1) {
			 prevButton.setEnabled(prev == 1);
		 }
	     if (next != -1) {
	    	 nextButton.setEnabled(next == 1);
	    	 finishButton.setEnabled(next == 1);
	     }
	     if (cancel != -1) {
	    	 cancelButton.setEnabled(cancel == 1);
	     }
	}

	/**
	 * Obtain the index of the current page being displayed in this
	 * wizard.
	 * 
	 * @return The zero-based index of the current page being displayed
	 * in this wizard.
	 */
	public int getCurrentPage() { return this.currentPage; }

	/**
	 * Obtain wizard page associated with the given index.
	 * 
	 * 
	 * @param index The zero-based index of the page being requested.
	 * 
	 * @return The page associated with this index. 
	 */
	public WizardPage getPage(int index) { return pages.get(index); }

	/**
	 * This method is called from the showWizard method just before the
	 * wizard is made visible. This method organizes the various components
	 * in the wizard in the appropriate locations.
	 * 
	 */
	protected void layoutWizard() {
		// Create the titles and layout the title panel.
		layoutTitlePanel();
		// Create the sequence labels in the sequence panel.
		layoutSequencePanel();
		// Create a JPanel to contain the pages and buttons
		JPanel centerPanel = new JPanel(new BorderLayout(0, 0));
		// Layout the button panel.
		centerPanel.add(layoutButtons(), BorderLayout.SOUTH);
		// Layout all the pages in the page panel.
		for(WizardPage page: pages) {
			String pageName = page.getTitle();
			pagePanel.add(page.getPage(), pageName);
		}
		// add the page panel.
		centerPanel.add(pagePanel, BorderLayout.CENTER);
		// Finally add the centerPanel to the center of this wizard
		add(centerPanel, BorderLayout.CENTER);
	}
	
	/**
	 * Helper method to layout the button panel to the bottom of the wizard.
	 * This method is invoked only once from the layoutWizard() method just 
	 * before the wizard is about to be displayed on the screen. This method
	 * lays out the buttons in a panel to the bottom of the wizard.
	 * 
	 * @return The button panel containing the buttons.
	 */
	protected Component layoutButtons() {
		// Create a button panel to place all the buttons.
		JPanel btnPanel = new JPanel();
		btnPanel.setLayout(new BoxLayout(btnPanel, BoxLayout.X_AXIS));
		// Add the buttons such that they are right justified.
		btnPanel.add(Box.createHorizontalGlue());
		btnPanel.add(prevButton); btnPanel.add(Box.createHorizontalStrut(5));
		btnPanel.add(nextButton); 
		btnPanel.add(finishButton); btnPanel.add(Box.createHorizontalStrut(5));
		if (helpURL != null) {
			helpButton.setEnabled(true);
			btnPanel.add(helpButton); 
			btnPanel.add(Box.createHorizontalStrut(5));
		}
		btnPanel.add(cancelButton); btnPanel.add(Box.createHorizontalStrut(5));
		// A a empty border to make the buttons panel look good
		CompoundBorder cb = new CompoundBorder(new CustomBorder("dnnn"),
				new CustomBorder("Lnnn"));
		btnPanel.setBorder(new CompoundBorder(cb, new EmptyBorder(15, 5, 15, 5)));
		// Return the created button panel
		return btnPanel;
	}
	
	/**
	 * Helper method to layout the sequence panel in the wizard. This method is
	 * invoked only once from the layoutWizard() method just before the wizard
	 * is about to be displayed on the screen. This method creates one label
	 * for each wizard page and organizes them vertically in the sequence
	 * panel.
	 */
	protected void layoutSequencePanel() {
		// Add some blank space at the top to make things pretty.
		sequencePanel.setBorder(new EmptyBorder(25, 5, 5, 0));
		sequencePanel.setUseImageSize(false, false);
		// Get the empty black box icon for each step.
		ImageIcon box = Utilities.getIcon("images/16x16/BoxWhite.png");
		// Create a sequence of labels for each wizard page.
		seqTitles = new JLabel[pages.size()];
		for(int pageID = 0; (pageID < pages.size()); pageID++) {
			WizardPage page = pages.get(pageID);
			// Create label with title.
			JLabel label = new JLabel(page.getTitle(), box, JLabel.LEFT);
			label.setForeground(Color.white);
			label.setBorder(new EmptyBorder(7, 7, 7, 25));
			// Add label to sequence pane.
			sequencePanel.add(label);
			// Save the label for future reference.
			seqTitles[pageID] = label;
		}
		// Add sequence panel to the left of the wizard
		add(sequencePanel, BorderLayout.WEST);
	}
	
	/**
	 * Helper method to layout the title panel in the wizard. This method
	 * is invoked only once from the layoutWizard() method just before
	 * the wizard is about to be displayed on the screen. This method
	 * creates and lays out the title and the sub-title components for
	 * this wizard.
	 */
	protected void layoutTitlePanel() {
		titles = new JLabel[2];
		for(int i = 0; (i < 2); i++) {
			titles[i] = new JLabel("", JLabel.LEFT);
			titles[i].setAlignmentX(0);
		}
		// Set the main title to have a larger font.
		Utilities.adjustFont(titles[0], 2, 12, 1);
		// Add titles to the title pane. Recollect that title pane has
		// a vertical box layout. Add vertical glue around titles to
		// center them in the title panel.
		Box titleBox = new Box(BoxLayout.Y_AXIS);
		titleBox.setBorder(new EmptyBorder(1, 10, 1, 1));
		titleBox.setBackground(Color.white);
		// Layout the titles and vertically center it.
		titleBox.add(Box.createVerticalGlue());
		titleBox.add(titles[0]);
		titleBox.add(Box.createVerticalStrut(5));
		titleBox.add(titles[1]);
		titleBox.add(Box.createVerticalGlue());
		// Add titles in the title panel.
		titlePanel.add(titleBox, BorderLayout.CENTER);
		titlePanel.setBorder(new CompoundBorder(new CustomBorder("nnLn"),
				new CustomBorder("nnDn")));
		add(titlePanel, BorderLayout.NORTH);
	}
	
	/**
	 * Helper method invoked when user clicks cancel button.
	 * 
	 * This is a helper method that derived classes can override. This method is
	 * invoked when the user clicks the cancel button in the wizard. This method
	 * is used to display a confirmation dialog to ensure that the user really
	 * wants to exit out of the wizard.
	 * 
	 * <p>
	 * <b>Note:</b>Derived classes that override this method must call the base
	 * class if the user wishes to quit so that any worker threads associated
	 * with this wizard are interrupted and deleted (and they don't simply hang
	 * around).
	 * </p>
	 * 
	 * @return This method always returns true to indicate the user wants to
	 *         quit.
	 */
	protected synchronized boolean cancel() {
		// Interrupt any worker threads currently working.
		for(Thread t: workerThreads) {
			// Interrupt the thread.
			t.interrupt();
		}
		return true;
	}
	
	/**
	 * Helper method invoked once the wizard is done.
	 * 
	 * This is a feedback method that derived classes can override.
	 * This method is invoked just after the wizard has completed
	 * its operations. This method can be used to trigger other
	 * operations after the wizard is done. 
	 * 
	 * <p><b>Note:</b>  The default method does absolutely nothing.</p>
	 *
	 * @param success This flag is true if the wizard completed
	 * successfully and the user clicks the finish button. Otherwise
	 * this flag is set ot false, indicating an premature exit.
	 */
	protected void done(boolean success) {
	}
	
	/**
	 * Helper method to perform the tasks related with page switching.
	 * 
	 * This is a helper method that is invoked whenever the pages
	 * in the wizard need to be switched. This method operates as
	 * follows:
	 * 
	 * <ol>
	 * 
	 * <li>It checks with the current page to ensure that the change
	 * in pages is OK.</li>
	 * 
	 * <li>If the page can be changed, it calls the changePage()
	 * method to actually change the page.</li>
	 * 
	 * </ol>
	 * 
	 * @param currPage Index of the current page being displayed.
	 * @param nextPage The next page that needs to be displayed.
	 */
	protected void checkChangePage(int currPage, int nextPage) {
		WizardPage page = pages.get(currPage);
		// Check with current page to see if we can change.
		if (!page.pageChanging(this, currPage, nextPage)) {
			// The page does not want to permit changing.
			return;
		}
		// Use helper method to actually do the change
		changePage(currPage, nextPage);
	}

	/**
	 * Helper method to perform the tasks related with page switching.
	 * 
	 * This method can be used whenever the pages
	 * in the wizard need to be switched. This method operates as
	 * follows:
	 * 
	 * <ol>
	 * 
	 * <li>It switches pages and notifies the new page that
	 * the change has occurred.</li>
	 * 
	 * <li>It updates the wizard display.</li>
	 * 
	 * </ol>
	 * 
	 * @param currPage Index of the current page being displayed.
	 * @param nextPage The next page that needs to be displayed.
	 */
	public void changePage(int currPage, int nextPage) {
		WizardPage page = pages.get(currPage);
		// If this was the last page we have nothing further to do.
		if (nextPage == pages.size()) {
			setVisible(false);
			dispose();
			// Notify that the wizard is done.
			done(true);
			return;
		}
		// Update the current page and let  the new page know
		// the change has occurred.
		int prevPage = currPage;
		currentPage  = nextPage;
		// Now update the status of the buttons.
		setButtonStatus();
		// Let the new page know about the change.
		page = pages.get(currentPage);
		page.pageChanged(this, currentPage, prevPage);
		// Update the sequence display
		updateSequence(currentPage);
		// Show the previous or next wizard page.
		CardLayout layout = (CardLayout) pagePanel.getLayout();
		layout.show(pagePanel, page.getTitle());
	}
	
	@Override
	public void actionPerformed(ActionEvent event) {
		String cmd = event.getActionCommand();
		if ("cancel".equals(cmd)) {
			if (cancel()) {
				// The user really want's to quit. So bail out.
				setVisible(false);
				dispose();
				// Notify that the wizard is done.
				done(false);
			}
		} else if ("help".equals(cmd)) {
			// The user clicked the help button. Launch help.
			HelpHandler.showHelp(getParent(), helpURL);
		} else if ("finish".equals(cmd)) {
			// The user clicked finish and the wizard is done.
			checkChangePage(currentPage, currentPage + 1);
		} else if ("prevnext".indexOf(cmd) > -1) {
			// When control drops here either it is the next or prev
			// buttons that were clicked. Add +1 or -1 to current page
			// (depending on button click) to determine next page.
			int nextPage = currentPage + ("prev".equals(cmd) ? -1 : 1);
			checkChangePage(currentPage, nextPage);
		}
	}
	
	/** Obtain the page pane to which wizard pages must be added.
	 * 
	 *  <p><b>Note:</b>The number of pages in this component and the currently
	 *  visible page determine the behavior or the buttons displayed
	 *  by this wizard.</p>
	 *  
	 * @return The page pane which contains all the pages of this wizard. 
	 */
	public CustomPanel getPagePane() { return pagePanel; }
		
	/**
	 * A custom panel that is used to hold all the pages of this wizard. This
	 * panel is added to the center of the wizard and uses a CardLayout
	 * to display one wizard page at a time.
	 */
	private CustomPanel pagePanel;

	/**
	 * A custom panel that is used as the title panel for this wizard. The
	 * title panel is displayed on top of the wizard.
	 */
	private CustomPanel titlePanel;
	
	/**
	 * The array of two titles that are displayed at the top of this
	 * wizard dialog. These values are automatically set whenever a new page
	 * is displayed by this wizard. The labels are created when the
	 * showWizard() method is invoked.
	 */
	private JLabel titles[];

	/**
	 * The array of sequence titles that are displayed to the left of the
	 * wizard dialog. These values are automatically updated whenever a 
	 * new page is displayed by this wizard. The labels are created when
	 * the showWizard() method is invoked.
	 */
	private JLabel seqTitles[];
	
	/**
	 * The "Prev" button that is displayed at the bottom of the wizard
	 * in a separate panel. This button is created in the constructor 
	 * but are organized into a panel only when the showWizard() method
	 * is called.
	 */
	private JButton prevButton;
	
	/**
	 * The "Next" button that is displayed at the bottom of the wizard
	 * in a separate panel. This button is created in the constructor 
	 * but are organized into a panel only when the showWizard() method
	 * is called.
	 */
	private JButton nextButton;

	/**
	 * The "Cancel" button that is displayed at the bottom of the wizard
	 * in a separate panel. This button is created in the constructor 
	 * but is organized into a panel only when the showWizard() method
	 * is called.
	 */
	private JButton cancelButton;

	/**
	 * The "Finish" button that is displayed at the bottom of the wizard
	 * in a separate panel. This button is created in the constructor 
	 * but are organized into a panel only when the showWizard() method
	 * is called. By default this button is hidden and is displayed
	 * only when the last panel is shown (instead of the next button).
	 */
	private JButton finishButton;

	/**
	 * The "Help" button that is displayed at the bottom of the wizard
	 * in a separate panel. This button is created in the constructor 
	 * but are organized into a panel only when the showWizard() method
	 * is called. By default this button is hidden and is displayed
	 * only if a valid help URL is specified.
	 */
	private JButton helpButton;
	
	/**
	 * A custom column panel that is displayed to the left of the wizard.
	 * This panel displays the sequence of configuration steps associated
	 * with this wizard.
	 */
	private CustomPanel sequencePanel;

	/**
	 * The array list containing the set of pages to be displayed by
	 * this wizard in the same sequence as they are to be displayed.
	 * Entries are appended to this list via the addPage() method.
	 */
	private ArrayList<WizardPage> pages;
	
	/**
	 * The current page that is being displayed by this wizard. This
	 * value is set to zero and works it way up to the last page 
	 * in this wizard.
	 */
	private int currentPage;

	/**
	 * The URL to the online help page associated with this wizard.
	 * This value is set via the showWizard() method. If the Help
	 * URL is set then the "Help" button is displayed.
	 */
	private String helpURL;

	/**
	 * The list of worker threads that various pages may have spawned.
	 * This list is maintained so that the threads can be killed when
	 * this wizard exits. Various pages may add and remove threads
	 * as they are created and destroyed.
	 */
	private ArrayList<Thread> workerThreads = new ArrayList<Thread>();
	
	/**
	 * Generated serial version ID for serializing this class. 
	 */
	private static final long serialVersionUID = -3960338752571957543L;
}
