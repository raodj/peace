package org.peace_tools.core;

/**
 * An interface to streamline creation of summary information.
 *
 * This interface has been introduced to streamline the process of creating
 * summary information to be displayed to the user at the end of wizard.
 * The summarizer provides a convenient interface that can be used by the
 * various components to create summary information.
 */
public interface SummaryWriter {
	/**
	 * Helper method to start a new section in summary.
	 * 
	 * This method is used to add a new summary section when creating summary
	 * information for a given job.
	 * 
	 * @param sectionTitle The title string for the section. A new line
	 * character is automatically added to the section title string to ensure
	 * consistent formatting. Leading and trailing white spaces are trimmed.
	 */
	public void addSection(String sectionTitle);
	
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
	 * line. This parameter can be null. 
	 */
	public void addSubSection(String name, String value, String desc);
	
	/**
	 * Helper method to add a new summary line.
	 * 
	 * This method is used to add a new line of summary information to a
	 * given summary content. Note that some or all of the parameters 
	 * can be null. It is up to the implementation to suitably handle,
	 * format, and display the summary information. 
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
	public void addSummary(String name, String value, String desc);
	
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
	public void addSubSummary(String name, String value, String desc);
}
