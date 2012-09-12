package org.peace_tools.generic;

import java.awt.Color;
import java.awt.Component;
import java.awt.Rectangle;

import javax.swing.JTextPane;
import javax.swing.SwingUtilities;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.plaf.ComponentUI;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultStyledDocument;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;

import org.peace_tools.core.PEACEInstaller;

/**
 * Helper class to display outputs (standard output & standard error)
 * from a process.
 * 
 * This class is a helper class that is used to display output streams
 * (may it be standard output or standard error) from a running process.
 * This class provides a styled document that provides different styles 
 * to display different forms of output. The styles are essentially 
 * extensions of Java's standard StyledDocument class. The current styles 
 * supported by this class include:
 * 
 * <table>
 * <tr><td>Style Kind</td><td>Description</td></tr>
 * <tr><td><tt>stdout</tt></td><td>Preferred style for standard output</td></tr>
 * <tr><td><tt>stderr</tt></td><td>Preferred style for standard error</td></tr>
 * <tr><td><tt>warning</tt></td><td>Preferred style for warnings (if any)</td></tr>
 * <tr><td><tt>info</tt></td><td>Preferred style for additional information</td></tr>
 * <tr><td><tt>heading</tt></td><td>Preferred style for line of heading</td></tr>
 * </table>
 * 
 * This class also creates a text pane that automatically tracks the
 * end of the document as content is appended to the display pane. This
 * is a convenience feature that permits the user to view the current output
 * (as in a standard terminal).
 * 
 * @see PEACEInstaller
 */
public class ProcessOutputDisplay implements DocumentListener {
	/**
	 * The styled document that is used to display the output from
	 * the process is different colors and styles
	 */
	private final DefaultStyledDocument outputDocument;
	
	/**
	 * The text pane where the output from the installation 
	 * process is displayed.
	 */
	private final JTextPane textPane;
	
	/**
	 * The various pre-defined styles that can be used to create
	 * entries in this class.
	 */
	public enum StyleKind {
		STDOUT, STDERR, WARNING, INFO, HEADING, BOLD_INFO
	}
	
	/**
	 * The constructor.
	 * 
	 * The constructor creates a styled document with the default
	 * set of styles and adds the styled document to a text pane
	 * created by the constructor.
	 */
	public ProcessOutputDisplay() {
		// Create the styled document to hold various formatted outputs
		outputDocument = createStyledDocument();
		// Create the text pane to display output logs
		textPane = new JTextPane(outputDocument) {
			private static final long serialVersionUID = 3247209363068354903L;
			@Override
			public boolean getScrollableTracksViewportWidth() {
			    Component parent = getParent();
			    ComponentUI ui = getUI();
			    return parent != null ? (ui.getPreferredSize(this).width <= parent
			        .getSize().width) : true;
			  }
		};
		textPane.setEditable(false);
	}
	
	/**
	 * Obtain the styled document encapsulated by this class.
	 * 
	 * @return The styled document created when this object
	 * was instantiated.
	 */
	public DefaultStyledDocument getDoc() {
		return outputDocument;
	}
	
	/**
	 * Obtain the text pane GUI element that contains the 
	 * styled document.
	 * 
	 * @return The text pane that actually displays the styled 
	 * document created when this object was instantiated.
	 */
	public JTextPane getTextPane() {
		return textPane;
	}
	
	/**
	 * Helper method to add string to output and scroll it.
	 * 
	 * This is a helper method that is called from various methods to
	 * show more detailed progress information to the user.
	 * 
	 * @param entry The log entry to be appended.
	 * 
	 * @param style The style (that determines color & font) for the entry.
	 */
	public void addLog(String entry, StyleKind style) {
		try {
			// Add the output in given style
			outputDocument.insertString(outputDocument.getLength(), entry, 
					outputDocument.getStyle(style.toString()));
			// Scroll the output to the bottom is automatically
			// done by the ProcessOutputDisplay class.
		} catch (BadLocationException e) {
			ProgrammerLog.log(e);
		}
	}
	
	/**
	 * Internal helper method to create a styled document.
	 * 
	 * This helper method was introduced to streamline the constructor.
	 * This method creates a styled document and sets the various
	 * default styles supported by the document.
	 * 
	 * @return The newly created styled document for use in the
	 * constructor.
	 */
	private DefaultStyledDocument createStyledDocument() {
		// Create our document with a suitable style context.
		StyleContext sc = new StyleContext();
		DefaultStyledDocument styleDoc = new DefaultStyledDocument(sc);
		styleDoc.addDocumentListener(this);
		// Setup some standard styles we are going to use.
		Style style = sc.addStyle(StyleKind.STDOUT.toString(), null);
		style.addAttribute(StyleConstants.Foreground, Color.green.darker());
		style = sc.addStyle(StyleKind.STDERR.toString(), null);
		style.addAttribute(StyleConstants.Foreground, Color.cyan);
		style = sc.addStyle(StyleKind.WARNING.toString(), null);
		style.addAttribute(StyleConstants.Foreground, Color.red);
		style = sc.addStyle(StyleKind.INFO.toString(), null);
		style.addAttribute(StyleConstants.Foreground, Color.blue);
		// Create heading style
		style = sc.addStyle(StyleKind.HEADING.toString(), null);
		style.addAttribute(StyleConstants.Foreground, Color.blue);
		style.addAttribute(StyleConstants.FontSize, new Integer(14));
		style.addAttribute(StyleConstants.FontConstants.Bold, Boolean.TRUE);
		// Create the bold info style
		style = sc.addStyle(StyleKind.BOLD_INFO.toString(), null);
		style.addAttribute(StyleConstants.Foreground, Color.black);
		style.addAttribute(StyleConstants.FontConstants.Bold, Boolean.TRUE);
		
		// Return the newly created style document.
		return styleDoc;
	}
	
	/**
	 * Implementation for method in DocumentListener to scroll to
	 * end of the document. This is not the best of implementation
	 * for this method but is sufficient for this class.
	 * 
	 * @param arg0 The document associated with the call. Currently,
	 * this event is not used.
	 */
	@Override
	public void changedUpdate(DocumentEvent arg0) {
		insertUpdate(arg0);
	}

	/**
	 * Implementation for method in DocumentListener to scroll to
	 * end of the document. This is not the best of implementation
	 * for this method but is sufficient for this class.
	 * 
	 * @param arg0 The document associated with the call. Currently,
	 * this event is not used.
	 */
	@Override
	public void insertUpdate(DocumentEvent arg0) {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				textPane.scrollRectToVisible(new Rectangle(0,
						textPane.getHeight()-2, 1, 1));	
			}
		});
	}

	/**
	 * Implementation for method in DocumentListener to scroll to
	 * end of the document. This is not the best of implementation
	 * for this method but is sufficient for this class.
	 * 
	 * @param arg0 The document associated with the call. Currently,
	 * this event is not used.
	 */
	@Override
	public void removeUpdate(DocumentEvent arg0) {
		insertUpdate(arg0);
	}
}
