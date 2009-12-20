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
			url = GenericHTMLView.class.getResource(fileOrURL);
		}
		// Force loading HTML page synchronously to work around Java bug
		JEditorPane htmlPane = new JEditorPane();
		htmlPane.setEditorKit(getEditorKit());
		htmlPane.setPage(url);
		htmlPane.setDragEnabled(false);
		htmlPane.setDropTarget(null);
		htmlPane.setEditable(false);
		htmlPane.addHyperlinkListener(this);
		// Wrap the editor pane in a scroll pane.
		JScrollPane jsp = new JScrollPane(htmlPane);
		// Finally add scroll pane to ourselves
		add(jsp, BorderLayout.CENTER);
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
	private EditorKit getEditorKit() {
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
	 * @param event The hyper link event to be processed by this method.
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
				mainFrame.showHelp(targetURL);
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
