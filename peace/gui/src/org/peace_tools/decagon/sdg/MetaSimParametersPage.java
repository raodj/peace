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

package org.peace_tools.decagon.sdg;

import java.awt.BorderLayout;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;

import org.peace_tools.decagon.helpers.ParameterSetGUIHelper;
import org.peace_tools.decagon.jaxb.Parameter;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.SeqTechType;

/**
 * Generic page to obtain a set of parameters used to run MetaSim.
 * 
 * This class serves as a generic page to obtain inputs from the
 * user based on a parameter definition file. The parameter definitions
 * are loaded by the {@link SyntheticDataSetGenerator#create(org.peace_tools.core.MainFrame)
 * method. This class facilitates displaying the parameters within
 * the wizard and performing standard operations. Several instances of
 * this class is created by the SyntehticDataSetGenerator to obtain
 * different sets of parameters.
 */
public class MetaSimParametersPage extends GenericWizardPage {

	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components are created by the
	 * {@link ParameterSetGUIHelper#createGUI(JPanel)} method.
	 * 
	 * @param wizard The wizard that logically owns this page. This
	 * parameter cannot be null.
	 * 
	 * @param psgh The set of parameters to be managed by this
	 * page. This parameter cannot be null.
	 * 
	 * @param techType The type of sequencing-technology for which this class 
	 * is managing parameters. This value must be one of <code>Sanger</code>,
	 * <code>454</code>, or <code>Illumina</code>. This value is also reflected
	 * as the prefixes for various parameter-names in the corresponding parameter
	 * configuration XML file.
	 * 
	 * @param subTitle A sub-title to be displayed for this wizard page.
	 */
	public MetaSimParametersPage(SyntheticDataSetGenerator wizard, 
			ParameterSetGUIHelper psgh, final SeqTechType techType,
			final String subTitle) {
		this.wizard   = wizard;
		this.psgh     = psgh;
		this.techType = techType;
		
		assert(this.wizard   != null);
		assert(this.psgh     != null);
		assert(this.techType != null);
		
		if (psgh.getParameterValueAsBoolean("Enable" + techType) == null) {
			throw new IllegalArgumentException("The parameter-set and the supplied " +
					"techType don't seem to match. Check to ensure that the parameter- " + 
					"names in the XML file all start with " + techType + " prefix.");
		}

		// Setup the title(s) for this page and border
		final String mainTitle = techType + " Settings";
		setTitle(mainTitle, subTitle);
		setBorder(new EmptyBorder(5, 5, 5, 5));
		
		// Create the GUI fields for parameters.
		JPanel subPanel = psgh.createGUI(null);
		// subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Pack the parameters into a scroll pane in case they are too big
		JScrollPane jsp = new JScrollPane(subPanel);
		jsp.setBorder(null);
		// Set the contents to this page
		add(jsp, BorderLayout.CENTER);
	}
	
	/**
	 * Determine if synthetic reads are to be generated.
	 * 
	 * This is a convenience method that can be used to determine if
	 * the user has requested generation of synthetic reads for
	 * the given technology type.
	 * 
	 * @return This method returns true if the reads are to be
	 * generated. Otherwise this method returns false.
	 */
	protected boolean isTechEnabled() {
		return psgh.getParameterValueAsBoolean("Enable" + techType);
	}

	/**
	 * Obtain the sequencing-technology associated with this page.
	 * 
	 * @return The type of sequencing-technology for which this class 
	 * is managing parameters. 
	 */
	protected SeqTechType getTechType() {
		return techType;
	}
	
	/**
	 * Obtain the parameter GUI helper set for this page.
	 * 
	 * @return The parameter GUI helper set for this page when
	 * it was created. The returned value is never null.
	 */
	protected ParameterSetGUIHelper getParameterSetGUIHelper() {
		return psgh;
	}
	
	/**
	 * Helper method to create read-count summary information.
	 * 
	 * This is a helper method that is invoked from the 
	 * {@link SyntheticDataSetGenerator#getSummary()} method to generate
	 * read-count summary information. This method assumes that 
	 * the supplied string builder object already has necessary
	 * HTML formatting in it and adds a table-row entry to the
	 * count summary. This method is invoked three times, once for
	 * each of the following technologies: Sanger, 454, and Illumina.
	 * If the given type of technology is not enabled, then this
	 * method does not create a summary line for it.
	 * 
	 * @param useMatePairForCount If this is true, then the mate-pair
	 * probability value is used to appropriately increase the net number
	 * of reads.
	 *  
	 * @param sb The string builder object to which a HTML-table-row entry
	 * is to be added if the technology is enabled. If this parameter
	 * is null then entries are not added to this string builder.
	 * 
	 * @return This method returns the number of reads for the given
	 * technology type that will be present in the final generated file.
	 * This value is used to report total number of reads in the generated
	 * file. If the technology is not enabled, then this method returns
	 * -1. 
	 */
	protected long getCountSummary(final boolean useMatePairForCount, StringBuilder sb) {
		if (!psgh.getParameterValueAsBoolean("Enable" + techType)) {
			// Generation of synthetic reads is disabled. Count summary
			// is not required.
			return -1;
		}
		// The technology is enabled. Get basic number of reads.
		long readCount = psgh.getParameterValueAsNumber(techType + "Count").longValue();
		if (useMatePairForCount) {
			// Use mate-pair probability to scale count of number of reads.
			// Typically this is done only for sanger technology.
			final double matePairProb = psgh.getParameterValueAsNumber(techType + "MatePairProb").doubleValue();
			readCount += Math.round(readCount * matePairProb);
		}
		// Get the average length of each read and variance in length.
		String avgReadLenStr = psgh.getParameter(techType + "AvgLen").getValue();
		if (psgh.getParameter(techType + "LenSD") != null) {
			// We have variance in length for this technology
			avgReadLenStr += " (&plusmn; " + 
				psgh.getParameter(techType + "LenSD").getValue() + ")"; 
		}
		// Form the row to be displayed in the supplied string builder
		if (sb != null) {
			sb.append("<tr><td>" + techType + "-type reads</td><td>~" +
					readCount + "</td><td>" + avgReadLenStr + " nt</td></tr>");
		}
		// Return number of reads added.
		return readCount;
	}
	
	/**
	 * Convenience method used to display information about parameters.
	 * 
	 * This method is invoked from the {@link SyntheticDataSetGenerator#getSummary()}
	 * method to generate summary information about the parameters for
	 * the given technology type.
	 * 
	 * @return A fragment of HTML containing summary information about the 
	 * parameters in this class. If sequences for this technology are not
	 * enabled, then this method returns an empty string.
	 */
	protected String getSummary() {
		if (!psgh.getParameterValueAsBoolean("Enable" + techType)) {
			// Generation of synthetic reads is disabled for this technology.
			return "";
		}
		// Create summary information.
		StringBuilder sb = new StringBuilder();
		final String heading = "Parameters for " + techType + 
			"-type reads";
		sb.append("<table><tr><td colspan=2 bgColor='lightGrey'><b>" + heading + "</b></td></tr>");
		for(Parameter param: psgh.getParameters()) {
			if (param.isController() || (param.getSummary() == null) ||
				param.getSummary().isEmpty() || (param.isHidden())) {
				// Empty string or this parameter is not a controller.
				continue;
			}
			sb.append("<tr><td>" + param.getSummary() + "</td>");
			sb.append("<td>" + param.getValue() + "</td></tr>\n");
		}
		sb.append("</table>");
		return sb.toString();
	}
	
	/**
	 * Method to setup this page just before it is displayed.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before this page is displayed. 
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		// Nothing to be done for now.
	}
	
	/**
	 * Method to ensure a valid data set is selected before proceeding to
	 * the next page.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before the user navigates to a different page. 
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int prevPage) {
		// Proceed to next page.
		return true;
	}

	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final SyntheticDataSetGenerator wizard;
	
	/**
	 * The set of parameters to be managed by this page. This
	 * value is set in the constructor.
	 */
	private final ParameterSetGUIHelper psgh;

	/**
	 * The type of sequencing-technology for which this class is managing
	 * parameters. This value must be one of <code>Sanger</code>, <code>454</code>, 
	 * or <code>Illumina</code>. This value is also reflected as the
	 * prefixes for various parameter-names in the corresponding parameter
	 * configuration XML file.
	 */
	private final SeqTechType techType;
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 2834710857108474822L;
}
