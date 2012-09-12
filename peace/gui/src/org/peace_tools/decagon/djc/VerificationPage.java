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

package org.peace_tools.decagon.djc;

import java.awt.BorderLayout;
import java.util.List;

import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.border.EmptyBorder;
import javax.swing.text.DefaultStyledDocument;

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.decagon.DecagonVariables;
import org.peace_tools.decagon.jaxb.Job;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.views.GenericHTMLView;
import org.peace_tools.workspace.Workspace;

/**
 * This wizard page displays read-only information about the job(s) 
 * to be created for the pipeline. This page provides the user with
 * a brief summary of the jobs that are going to be created for
 * the various steps in the assemblers pipeline. This is the last 
 * page that the user can back track or cancel the job. In the next
 * page the jobs are created on the target server and started.
 * 
 */
public class VerificationPage extends GenericWizardPage implements SummaryWriter {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: a couple of labels
	 * and a text area summary information about the jobs.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public VerificationPage(DecagonJobCreator wizard) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Job(s) Summary", 
				"Verify information and submit job(s)");
		setBorder(new EmptyBorder(5, 5, 5, 5));

		// Create the informational labels.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/32x32/Information.png"), 
				JLabel.LEFT);
		JLabel caution = new JLabel(CAUTION_MSG, 
				Utilities.getIcon("images/24x24/Warning.png"), 
				JLabel.LEFT);
		
		// Create the summary tabs with styled documents in which
		// summary information is to be displayed (later on)
		JTabbedPane summaryTab = new JTabbedPane();
		jobSummary = createHTMLDocument(summaryTab, 
					"Job(s) Summary", "JobSet.png");
		// Display list of DADX parameters for user's reference.
		paramSummary = createHTMLDocument(summaryTab, 
				"DADX Parameters", "EAST.png");
		// Display list of variables for user's reference.
		varSummary = createHTMLDocument(summaryTab, 
				"DECAGON Variables", "EAST.png");

		// Pack the display fields into a suitable panel
		JPanel subPanel = new JPanel(new BorderLayout(5, 5));
		subPanel.add(info, BorderLayout.NORTH);
		subPanel.add(summaryTab, BorderLayout.CENTER);
		subPanel.add(caution, BorderLayout.SOUTH);
		// Finally add the sub panel to this panel.
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * Helper method to create a styled document to log summary
	 * information. This helper method has been introduced to
	 * streamline the constructor and keep the code clutter to
	 * a minimum. This method performs the following tasks:
	 * 
	 * <ol>
	 * <li> It creates a styled document and configures the 
	 * document's custom styles.</li>
	 * <li>It wrap the styled document into a suitably configured
	 * scroll pane.</li>
	 * <li>It adds the scroll pane to the provided tab panel
	 * with supplied (via parameters) information. 
	 * 
	 * @param parentPanel The tab panel to which the scroll pane
	 * containing the styled document created by this method should
	 * be added. This parameter cannot be null.
	 * 
	 * @param tabIconName The name of the image file that contains
	 * the icon to be used for the tab to be created and added to 
	 * the parentPanel. The file name is prefixed with the path 
	 * prefix "images/16x16/" by this method.
	 * 
	 * @return The styled document to which summary information
	 * is to be added.
	 */
	private JEditorPane createHTMLDocument(JTabbedPane parentPanel,
			String tabName, String tabIconName) {
		// Create suitable HTML display panel
		JEditorPane htmlDoc = GenericHTMLView.createHTMLPane();
		htmlDoc.putClientProperty(JEditorPane.HONOR_DISPLAY_PROPERTIES, Boolean.TRUE);
		htmlDoc.setFont(getFont());

		// Add the text pane to a scroll pane to handle large outputs
		JScrollPane jsp = new JScrollPane(htmlDoc);
		// Add the scroll pane to the supplied tab pane as a tab
		parentPanel.addTab(tabName, Utilities.getIcon("images/16x16/" + tabIconName), jsp);
		return htmlDoc;
	}

	/**
	 * Helper method to start a new section in summary.
	 * 
	 * This method is used to add a new summary section when creating summary
	 * information for a given job.
	 * 
	 * @param currentSummary The string builder (that is typically passed in as a parameter
	 * from this class to the EASTJobWizard) to which an HTML-row of section
	 * title is to be written.
	 * 
	 * @param sectionTitle The title string for the section. A new line
	 * character is automatically added to the section title string to ensure
	 * consistent formatting. Leading and trailing white spaces are trimmed.
	 */
	public void addSection(String sectionTitle) {
		currentSummary.append("<tr><td colspan=\"2\">");
		currentSummary.append("<font size=\"+1\"><b>");
		currentSummary.append(sectionTitle);
		currentSummary.append("</b></font></td></tr>");
		// Reset flag to indicate first line of section
		oddLine = true;
	}
	
	/**
	 * Helper method to add a new summary line.
	 * 
	 * This method is used to add a new line of summary information to a
	 * given summary content. 
	 * 
	 * @param sb The string builder (that is typically passed in as a parameter
	 * from this class to the EASTJobWizard) to which an HTML-row of summary
	 * information is to be written.
	 * 
	 * @param name The name of the parameter or property for which a summary
	 * line is being added to the given dsd. This parameter can be null.
	 * 
	 * @param value The value for the given parameter to be added to the dsd.
	 * This parameter can be null.
	 * 
	 * @param desc An optional description to be associated with the summary
	 * line. This parameter can be null. 
	 */
	@Override
	public void addSummary(String name,	String value, String desc) {
		addSummary(name, value, desc, "black", false, false);
	}
	
	/**
	 * Helper method to add a new sub-section summary line.
	 * 
	 * This method is used to add a new line of summary information to a
	 * given summary content. Except that, this method permits the 
	 * implementation to suitably highlight (if possible) the sub-section
	 * entry. 
	 * 
	 * <p>Note that some or all of the parameters 
	 * can be null. It is up to the implementation to suitably handle,
	 * format, and display the summary information.</p> 
	 * 
	 * @param name The name of the parameter or property for which a summary
	 * line is being added. This parameter can be null.
	 * 
	 * @param value The value for the given parameter to be added to the dsd.
	 * This parameter can be null.
	 * 
	 * @param desc An optional description to be associated with the summary
	 * line. This parameter can be null as it is ignored by this method as
	 * this style does not use it. 
	 */
	public void addSubSection(String name, String value, String desc) {
		addSummary(name, null, null, "black", false, true);
	}
	
	/**
	 * Helper method to add a new sub-summary line.
	 * 
	 * This method is used to add a new line of summary information to 
	 * logically appear under a sub-section. The need for a sub-summary
	 * line to appear under a summary is not enforced. It is up to the
	 * caller to ensure a suitable sub-section summary is created 
	 * first. 
	 * 
	 * <p>Note that some or all of the parameters  can be null. It is up
	 * to the implementation to suitably handle, format, and display the
	 * summary information.</p> 
	 * 
	 * @param name The name of the parameter or property for which a summary
	 * line is being added. This parameter can be null.
	 * 
	 * @param value The value for the given parameter to be added to the dsd.
	 * This parameter can be null.
	 * 
	 * @param desc An optional description to be associated with the summary
	 * line. This parameter can be null. 
	 */
	public void addSubSummary(String name, String value, String desc) {
		addSummary(name, value, desc, "black", true, false);
	}
	
	/**
	 * Internal helper method to add a new summary line.
	 * 
	 * This method is used by the various summary generation methods in
	 * this class (with different parameters) generate suitable summary
	 * lines in HTML format. 
	 * 
	 * @param sb The string builder (that is typically passed in as a parameter
	 * from this class to the EASTJobWizard) to which an HTML-row of summary
	 * information is to be written.
	 * 
	 * @param name The name of the parameter or property for which a summary
	 * line is being added to the given dsd. This parameter can be null.
	 * 
	 * @param value The value for the given parameter to be added to the dsd.
	 * This parameter can be null.
	 * 
	 * @param desc An optional description to be associated with the summary
	 * line. This parameter can be null.
	 * 
	 *  @param bgColor An optional coloring scheme for the entire row of the
	 *  summary. This parameter cannot be null.
	 *  
	 *  @param indent If this flag is true then the summary information is
	 *  placed in the second column of the summary-table (see {@link #HTML_BEGIN}).
	 *  Otherwise the summary information spans the two columns of the table.
	 *  
	 *  @param bold Flag to indicate if the characters in the line must be in
	 *  bold/underlined (to highlight a sub-section)
	 */
	private void addSummary(String name, String value, String desc, String bgColor, 
			boolean indent, boolean bold) {
		// Create the summary information to be placed in a column(s)
		StringBuilder sb = new StringBuilder(1024);
		if (bold) {
			sb.append("<b><u>");
		}
		if (name != null) {
			sb.append(name);
		}
		if (value != null) {
			if (name != null) {
				sb.append(':');
			}
			sb.append("<i>"); sb.append(value); sb.append("</i>");
		}
		if (desc != null) {
			if (sb.length() > 0) {
				// There was either name or value added earlier. Place
				// description on the next line.
				sb.append("<br/>");
			}
			sb.append("<font size=\"-2\">");
			sb.append(desc);
			sb.append("</font>");
		}
		if (bold) {
			sb.append("</u></b>");
		}

		// Based on indent flag place the summary information (in sb)
		// in two columns or just the second column.
		currentSummary.append("<tr>");
		if (indent) {
			currentSummary.append("<td/><td>");
			currentSummary.append(sb);
			currentSummary.append("</td>");
		} else {
			currentSummary.append("<td colspan=\"2\">");
			currentSummary.append(sb);
			currentSummary.append("</td>");
		}
		currentSummary.append("</tr>");
		// Toggle odd line flag
		oddLine = !oddLine;
	}
	
	/**
	 * Method to setup this page just before it is displayed.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before this page is displayed. This method populates the
	 * summary information using the data supplied by the user
	 * in the various wizard pages (along with default values).
	 * 
	 * @param dialog The wizard dialog to which this page logically
	 * belongs. Currently this parameter is not used and the
	 * {@link #wizard} object is used directly.
	 * 
	 * @param currPage The zero-based logical index of this page
	 * in the sequence of pages on this wizard.
	 * 
	 * @param prevPage The zero-based logical index of the page 
	 * from where the user navigated to this page.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Create the overall job summary first.
		currentSummary = new StringBuilder(1024);
		currentSummary.append(HTML_BEGIN);
		wizard.getDataSet().summarize(this);
		addSummary(HORIZ_LINE, null, null);
		// Create summary information about each job in the pipeline. For this
		// we obtain the job information and job-server-configuration details
		final JobConfigPage[] jobConfigs = wizard.getJobConfigs();
		final List<Job>       jobList    = wizard.getDADXHelper().getAssemblerDetails().getJob();
		assert ( jobConfigs.length == jobList.size());
		long jobID = Workspace.get().getJobList().getSeqCounter();
		for(int i = 0; (i < jobList.size()); i++, jobID++) {
			// Setup a tentative job ID for the job(s) for now.
			DecagonVariables decVars = wizard.getDADXHelper().getVariables();
			decVars.addVariable(DecagonVariables.JOB_ID, "job" + jobID);
			wizard.getDADXHelper().summarize(this, jobList.get(i), 
					jobConfigs[i].getServerInfoPanel());
			addSummary(HORIZ_LINE, null, null);
		}		
		currentSummary.append(HTML_END);
		jobSummary.setText(currentSummary.toString());
		jobSummary.setCaretPosition(0);
		
		// Create the DADX parameter summary information next.
		currentSummary = new StringBuilder(1024);
		currentSummary.append(HTML_BEGIN);
		wizard.getDADXHelper().summarize(this, false);
		currentSummary.append(HTML_END);
		paramSummary.setText(currentSummary.toString());
		paramSummary.setCaretPosition(0);
		
		// Create the DECAGON variables summary information next.
		currentSummary = new StringBuilder(1024);
		currentSummary.append(HTML_BEGIN);
		wizard.getDADXHelper().getVariables().summarize(this, false);
		currentSummary.append(HTML_END);
		varSummary.setText(currentSummary.toString());
		varSummary.setCaretPosition(0);
		
		// Let the super class do any additional work.
		super.pageChanged(dialog, currPage, prevPage);
	}

	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final DecagonJobCreator wizard;

	/**
	 * The HTML document that is used to display the summary information 
	 * about various jobs with a pretty style.
	 */
	private JEditorPane jobSummary;
	
	/**
	 * The HTML document that is used to display the summary information 
	 * about various DECAGON variables with a pretty style.
	 */
	private JEditorPane varSummary;
	
	/**
	 * The HTML document that is used to display the summary information 
	 * about various DADX parameter values with a pretty style. 
	 */
	private JEditorPane paramSummary;
	
	/**
	 * This is an intermediate/convenience object that is used to 
	 * gather HTML summary information from various methods that
	 * implement the  SummaryWriter interface. The data gathered
	 * in this object is used to populate the {@link #jobSummary},
	 * {@link #varSummary}, and {@link #paramSummary} objects.  
	 */
	private transient StringBuilder currentSummary;
	
	/**
	 * Flag to indicate if the line being added to a summary document
	 * is an odd or an event line. This flag works based on the 
	 * assumption that a single complete section of information is
	 * added to a given document. This flag is reset to true by the
	 * {@link #addSection(DefaultStyledDocument, String)} method and
	 * toggled by the {@link #addSummaryLine(DefaultStyledDocument, String)}
	 * method.
	 */
	private boolean oddLine;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some contextual information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"<html>You have provided the following information regarding<br/>" +
		"the job(s) to be scheduled. Please verify the information. Next<br/>" +
		"the wizard will perform the initial job setup.</html>";

	/**
	 * A simple caution message that is displayed at the bottom of the
	 * wizard page.
	 */
	private static final String CAUTION_MSG = 
		"<html><b>On clicking the 'Next' button the job(s) will be submitted!</b><br>" +
		"You will not be able to backtrack from the next page.</html>";

	/**
	 * Static text used to display a separator line in the summary information.
	 * This text intentionally does not use <hr/> tag as it did not look
	 * very good in the Java HTML display.
	 */
	private static final String HORIZ_LINE =  
		". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ." +
		". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .";
	

	/**
	 * A static string that is used to start a summary HTML document.
	 * Summary information is organized as a single HTML-table. The first
	 * column is merely used to indent sub-summaries.  Regular summaries
	 * are displayed to span two columns. Each summary entry occupies a
	 * row in the table. This string is used in the 
	 * {@link #pageChanged(WizardDialog, int, int)} method.
	 */
	private static final String HTML_BEGIN = 
		"<html><table cellspacing=\"2\" cellpadding=\"0\" width=\"600\">" +
		"<tr><td width=\"5%\"/><td/></tr>";
	
	/**
	 * A static string that is used to end a summary HTML document.
	 * This string is used in the {@link #pageChanged(WizardDialog, int, int)}
	 * method.
	 */
	private static final String HTML_END = "</table></html>";

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8538523942750752146L;
}
