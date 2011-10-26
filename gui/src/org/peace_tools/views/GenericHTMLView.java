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

package org.peace_tools.views;

import java.awt.BorderLayout;
import java.awt.Component;
import java.net.URL;

import javax.swing.JButton;
import javax.swing.JEditorPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToolBar;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.text.AbstractDocument;
import javax.swing.text.Document;
import javax.swing.text.EditorKit;
import javax.swing.text.html.HTMLEditorKit;

import org.peace_tools.core.MainFrame;
import org.peace_tools.generic.HelpHandler;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;

/**
 * A HTML-based view to display simple HTML data.
 * 
 * This is a simple view to display any HTML data. For instance,
 * this view is used to display the welcome.html file in 
 * the installFiles folder. This view displays the HTML file
 * in a JEditorPane. Note that the HTML rendering in the editor
 * pane is rather simple. Consequently, don't use complex tags
 * and elaborate formatting in HTML.
 */
public class GenericHTMLView extends JPanel implements HyperlinkListener {
	/**
	 * The default and only constructor.
	 * 
	 * The constructor performs all the tasks of creating and laying
	 * out the welcome screen. This includes loading the welcome
	 * HTML file.
	 * 
	 * @param fileOrURL The file name (which is assumed to be in PEACE
	 * source hierarchy/jar) or the HTTP URL to be displayed. Note that it
	 * is assumed that URLs are prefixed with http:// while files are not.
	 * 
	 * @param mainFrame A reference to the main frame that owns this view.
	 * This reference is used to launch default browser when links are clicked.
	 * 
	 * @exception This method throws an exception when the HTML view
	 * could not be successfully created.
	 */
	public GenericHTMLView(String fileOrURL, MainFrame mainFrame) throws Exception {
		super(new BorderLayout(0, 0));
		// Save reference to the main frame to launch urls
		this.mainFrame = mainFrame;
		// Check to see if the specified fileOrURL and handle it 
		URL url = null;
		if (fileOrURL.startsWith("http")) {
			url = new URL(fileOrURL);
		} else {
			// We need to check two different paths to handle difference between
			// running via a Jar file and running in development environment (Eclipse).
			// Try the Jar file choice first so that regular runs go faster.
			url = Utilities.class.getResource("/" + fileOrURL);
			if (url == null) {
				url = Utilities.class.getResource(Utilities.PATH_PREFIX + fileOrURL);
			}
			ProgrammerLog.log("HTML editor is loading: " + url);
		}
		// Force loading HTML page synchronously to work around Java bug
		JEditorPane htmlPane = createHTMLPane();
		htmlPane.setPage(url);
		htmlPane.addHyperlinkListener(this);
		// Wrap the editor pane in a scroll pane.
		JScrollPane jsp = new JScrollPane(htmlPane);
		// Finally add scroll pane to ourselves
		add(jsp, BorderLayout.CENTER);
	}

	/**
	 * Helper method to create a read-only panel to display HTML data.
	 * 
	 * This method is a helper method that has been exposed for use by
	 * other classes to create a suitable GUI component to display HTML
	 * data. This method creates a read-only panel in which HTML data
	 * can be displayed by calling setText() method
	 * or by specifying a URL via the setPage() method. In addition
	 * the component may also be configured to react to clicking of
	 * links via the addHyperlinkListener() method.
	 *  
	 * @return A suitably configured JEditorPane object in which HTML
	 * content can be displayed.
	 */
	public static JEditorPane createHTMLPane() {
		JEditorPane htmlPane = new JEditorPane();
		htmlPane.setEditorKit(getEditorKit());
		// htmlPane.setPage(url);
		htmlPane.setDragEnabled(false);
		htmlPane.setDropTarget(null);
		htmlPane.setEditable(false);
		// htmlPane.addHyperlinkListener(this);
		return htmlPane;
	}
	
	/**
	 * Create custom editor kit to work around Java bug.
	 * 
	 * This method is a helper method that is used to create a custom
	 * HTML editor kit. The custom editor kit forces the loading of
	 * the page to proceed synchronously to avoid some race condition
	 * in Java. This a known bug and this a work around provided
	 * by Sun.
	 * 
	 * @return An HTML editor kit that forces loading to proceed
	 * synchronously.
	 */
	private static EditorKit getEditorKit() {
		HTMLEditorKit kit = new HTMLEditorKit() {
			private static final long serialVersionUID = 1L;
			@Override
            public Document createDefaultDocument() {
                Document doc = super.createDefaultDocument();
                if (doc instanceof AbstractDocument) {
                    // force non-threaded loading, fixes CR 952223
                    // there are JEditorPane threading bugs, see  http://forums.sun.com/thread.jspa?threadID=5391050
                    ((AbstractDocument) doc).setAsynchronousLoadPriority(-1);
                }
                return doc;
            }
        };
        return kit;
	}
	
	/**
	 * Method to intercept hyper link clicks in HTML.
	 * 
	 * This method is invoked whenever the user clicks on a hyper link
	 * in an HTML document opened in this pane. This method intercepts
	 * the events and performs appropriate operations to handle special
	 * hyper links that translate to main menu options in PEACE.
	 * 
	 * @param hlEvent The hyper link event to be processed by this method.
	 */
	@Override
	public void hyperlinkUpdate(HyperlinkEvent hlEvent) {
		if (hlEvent.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
			String targetURL = hlEvent.getDescription();
			if (targetURL.startsWith("##")) {
				// This is a URL for a tool bar. Handle it separately.
				triggerTool(targetURL.substring(2));
			} else {
				// This is a standard URL. Load it in default browser.
				HelpHandler.showHelp(mainFrame, targetURL);
			}
		}
	}
	
	/**
	 * Method to trigger a tool bar button action.
	 * 
	 * This method is invoked whenever the user clicks on a hyper link
	 * that has two "##" in front it. The remainder of the string is
	 * assumed to be the action command associated with a tool bar 
	 * button. This method searches the list of tools in the tool bar
	 * for the given action and if it was found, it triggers the action.
	 * 
	 * @param actionCmd The action command to be processed by this
	 * method.
	 */
	private void triggerTool(String actionCmd) {
		JToolBar toolbar = mainFrame.getToolBar();
		for(int i = 0; (i < toolbar.getComponentCount()); i++) {
			Component tool = toolbar.getComponentAtIndex(i);
			if (tool instanceof JButton) {
				JButton button = (JButton) tool;
				if (actionCmd.equals(button.getActionCommand())) {
					// Found a valid action. Trigger it.
					button.doClick();
					break;
				}
			}
		}
	}
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * A generated serialization UID to keep the compiler happy. 
	 */
	private static final long serialVersionUID = -7553639166871732314L;
}
