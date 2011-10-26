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
