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

package org.peace_tools.core;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.DefaultListCellRenderer;
import javax.swing.Icon;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;

import org.peace_tools.generic.CustomPanel;
import org.peace_tools.generic.Utilities;

public class AboutDialog extends JDialog {
 	/** The constructor.
     * 
     * The constructor essentially sets up the default properties of
     * this dialog along with the top image for display. The core setup
     * of various tabs is done using helper methods.
     */
    public AboutDialog(JFrame owner) throws Exception {
    	super(owner, "About PEACE", true);
    	setLayout(new BorderLayout(0, 0));
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        this.getContentPane().setBackground(Color.white);
        setResizable(false);
        // Setup icon in title bar.
        this.setIconImage(Utilities.getIcon("images/16x16/PEACE.png").getImage());
        // Create a CustomPanel with the peace logo in it. This panel
        // will contain other components.
        cp = new CustomPanel(new BorderLayout(0, 0));
        cp.setBackground(Color.white);
        cp.setImage("images/peace_blue_header.png");
        // Use the image width to setup empty border for logo spacer
        JLabel imgHolder = new JLabel(Utilities.getIcon("images/peace_blue_header.png"));
        // Use the logo's width as the preferred width for this dialog.
        Dimension size = imgHolder.getPreferredSize();
        size.height = 475;
        setPreferredSize(size);
        setSize(size);
        // Add logo and tabs to the main dialog
        cp.add(Box.createVerticalStrut(100), BorderLayout.NORTH);
        // Create the various tabs in this about box.
        JTabbedPane tabs = new JTabbedPane();
        // Create and add a tab to the tabbed pane.
        createAuthorList(tabs);
        createLicenceTab(tabs);
        cp.add(tabs, BorderLayout.CENTER);
        // Add the custom pane to the dialog box
        add(cp, BorderLayout.CENTER);
    }
    
    /**
     * Helper method to load author information from various files.
     * 
     * This method is a helper method that was primarily introduced
     * to streamline the code better. This method uses the brief names
     * of the authors (statically fixed in this class) to load the
     * data and images for each one into arrays. The data is used to
     * create the list and the images are suitably rendered.
     * 
     * @note This method adds images for authors in the AuthorImage
     * array list.
     * 
     * @return The list of strings containing the author info.
     * 
     * @throws Exception This method throws exceptions if errors occur
     * when loading the information.
     */
    private String[] loadAuthorInfo() throws Exception {
    	final String Path = "installFiles/authors/";
    	ArrayList<String> authorInfo = new ArrayList<String>();
    	for(int i = 0; (i < AuthorList.length); i++) {
    		String infoFileName = Path + AuthorList[i] + ".txt";
    		String info = Utilities.readSmallTextFile(infoFileName);
    		info = "<html>" + info + "</html>";
    		authorInfo.add(info);
    		// Load the image for the author as well.
    		String imgFile = Path + AuthorList[i] + ".png";
    		AuthorImages.add(Utilities.getIcon(imgFile));
    	}
    	return authorInfo.toArray(new String[authorInfo.size()]);
    }
    
    /**
     * Helper method to create the "Authors" tab in the dialog.
     * 
     * This method is invoked from the constructor only once when the
     * dialog is created. This method was introduced to streamline the
     * code better and ease creation of the various tabs in the dialog.
     * 
     * @param tabs The tabbed pane to which the author tab is to be added.
     * 
     * @throws Exception This method propagates exception that may occur
     * when loading author information from various files.
     */
    private void createAuthorList(JTabbedPane tabs) throws Exception {
    	JList authorList = new JList(loadAuthorInfo());
    	authorList.setBackground(Color.white);
    	// Create a custom renderer for the list to display images
    	authorList.setCellRenderer(new DefaultListCellRenderer() {
    		private static final long serialVersionUID = 1304243355755229299L;
			/**
    		 * Override default implementation to setup the image
    		 * of the author adjacent to the description.
    		 */
    		@Override
    		public Component getListCellRendererComponent(JList list,
                    Object value,
                    int index,
                    boolean isSelected,
                    boolean cellHasFocus) {
    			// Let base class do the default setup
    			super.getListCellRendererComponent(list, value,
                        index, isSelected, cellHasFocus);
    			// Setup image to be displayed.
    			this.setIcon(AuthorImages.get(index));
    			this.setBorder(new EmptyBorder(3, 3, 3, 3));
    			// Return the label to be rendered
    			return this;
    		}
    	});
    	// Wrap the list in a scroll pane to let it scroll
    	JScrollPane jsp = new JScrollPane(authorList);
        // Create label to draw the title with a bit of border to
        // make the tab a bit larger than usual (for fancy)
        JLabel tabTitle = new JLabel("Authors", 
        		Utilities.getIcon("images/16x16/ProgLog.png"),
                        JLabel.LEFT);
        tabTitle.setBorder(new EmptyBorder(4, 0, 1, 0));
        tabs.add(jsp);
        tabs.setTabComponentAt(0, tabTitle);
    }
    
	/**
	 * This is a helper method that is used to create the license
	 * tab in the about dialog. This method was introduced to streamline
	 * the code. It is a relatively straightforward method.
	 * 
	 * @param tabs The tabbed pane to which the license tab is to be added.
     * 
     * @throws Exception This method propagates exception that may occur
     * when loading author information from various files.
	 */
	private void createLicenceTab(JTabbedPane tabs) throws Exception {
		JTextArea space = new JTextArea();
		space.setEditable(false);
		// Load a set the license in this space.
		space.setText(Utilities.readSmallTextFile("installFiles/License.txt"));
		space.setCaretPosition(0);
		// Wrap the license space in a scroll pane to let users
		// scroll through the license information.
		JScrollPane jsp = new JScrollPane(space);
		jsp.setBorder(new EtchedBorder(EtchedBorder.RAISED));
		
		tabs.addTab("License", 
				Utilities.getIcon("images/16x16/GPL.png"), jsp);
	}
	
    /**
     * A custom panel that actually has the PEACE logo as the 
     * background image. This panel is used to hold all the 
     * components just to make this dialog look nice.
     */
    private CustomPanel cp;

    /**
     * The list of authors and their information to be loaded from
     * text files in PEACE install directory.
     */
    private final String AuthorList[] = {
    	"raodm", "ozden", "karroje", "liang", "moler", "zhang"
    };
    
    /**
     * The icons to be displayed for each author to make the about
     * dialog a bit more fancy.
     */
    private final ArrayList<Icon> AuthorImages = 
    	new ArrayList<Icon>(10);
    
    /**
	 * A generated serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 7642494335728910205L;
}
