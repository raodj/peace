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

package org.peace_tools.workspace;

import java.io.PrintWriter;
import java.util.regex.Pattern;

import org.w3c.dom.Element;

/**
 * Database (DB) Classifier to distinguish ESTs from different databases. 
 *
 * This class corresponds to a "DBClassifier" element in a PEACE work 
 * space XML data. This class encapsulates the core information 
 * associated with a classifier.  Users create these entries by specifiying
 * a suitable regular expression to identify ESTs from a data base. The data
 * is persisted in the XML workspace for future reference so that users can
 * repeatedly reuse a classifier. 
 */
public class DBClassifier {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * DBClassifier entry. This method is typically  used to create a suitable
	 * entry when loading a Work space into the GUI.
	 * 
	 * @param dbNode The DOM element to be used for creating the classifier
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created classifier entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static DBClassifier create(Element jobNode) throws Exception {
		// First extract the necessary information from the DOM tree.
		String regExp   = DOMHelper.getStringValue(jobNode, "RegExp", true);
		regExp          = (regExp != null) ? regExp : "";
		String desc     = DOMHelper.getStringValue(jobNode, "Description", true);
		desc            = (desc != null) ? desc : "";
		String colorStr = DOMHelper.getStringValue(jobNode, "Color");
		String enableStr= DOMHelper.getStringValue(jobNode, "Enabled");
		// Convert color to an integer.
		int color       = Integer.parseInt(colorStr, 10);
		// Convert enabled string to a boolean
		boolean enabled = "true".equals(enableStr);
		// Now that we have sufficient information create the core
		// job object.
		DBClassifier classifier = new DBClassifier(regExp, desc, color, enabled);
		// Return the newly created object.
		return classifier;
	}
	
	/**
	 * Constructor to create a classifier object with the fixed value fields
	 * initialized to specific values.
	 * 
	 * @param regExp The regular expression to be associated with the given
	 * classifier.
	 * 
	 * @param description A user-supplied description for this entry. This
	 * maybe an empty string (but cannot be null).
	 * 
	 * @param color The RGB color code to be associated with this entry.
	 * 
	 * @param enabled Flag to indicate if this classifier is enabled. Only
	 * enabled classifiers are used for classifying ESTs.
	 */
	public DBClassifier(String regExp, String description, int color, boolean enabled) {
		this.description = description;
		this.regExp      = regExp;
		this.pattern     = null;
		this.color       = color;
		this.enabled     = enabled;
	}

	/**
	 * Obtain the user supplied description for this this classifier.
	 *  
	 * This value is set when new entries are added. The description is
	 * persisted in the work space configuration file and loaded when a 
	 * work space is opened in the GUI.
	 * 
	 * @return This method returns the user supplied description
	 * set for this job.
	 */
	public String getDescription() { return description; }

	/**
	 * Set a description for this entry 
	 * 
	 * The description is persisted in the work space configuration
	 * file and loaded when a work space is opened in the GUI.
	 * @param desc
	 */
	public void setDescription(String desc) {
		description = desc;
	}
	
	/**
	 * Set the regular expression for this DB Classifier entry.
	 * 
	 * This method must be used to set the regular expression to be
	 * associated with this classifier entry.
	 * 
	 * @note This method resets the pattern associated with this
	 * classifier.
	 * 
	 * @param regExp The regular expression to be set for this
	 * entry. It is assumed that the regular expression is valid
	 * and this method does not perform any specific validation.
	 */
	public void setRegExp(String regExp) {
		this.regExp  = regExp;
		this.pattern = null;
	}
	
	/**
	 * Obtain the raw regular expression set for this classifier.
	 * 
	 * This method returns the raw string expression set by the user
	 * for this classifier. To utilize the compiled regular expression
	 * use the getPattern() method in this class.
	 * 
	 * @return The raw regular expression set for this classifier.
	 */
	public String getRegExp() {
		return regExp;
	}

	/**
	 * Obtain the color code for shading database entry.
	 *  
	 * This value is set when new entries are added. The color is
	 * persisted in the work space configuration file and loaded when a 
	 * work space is opened in the GUI.
	 * 
	 * @return This method returns the user supplied color to be used
	 * to shade the class of ESTs identified by this classifier.
	 */
	public int getColor() { return color; }

	/**
	 * Set a color to shade ESTs identified by this classifier. 
	 * 
	 * The color is persisted in the work space configuration
	 * file and loaded when a work space is opened in the GUI.
	 * 
	 * @param color The RGB color code for the color to shade this
	 * data base entry.
	 */
	public void setColor(int color) {
		this.color = color;
	}

	/**
	 * Determine if this classifier is enabled.
	 * 
	 * Only classifiers that are enabled are to be used for classifying
	 * ESTs.
	 * 
	 * @return This method returns true to indicate that this classifier
	 * is enabled.
	 */
	public boolean isEnabled() {
		return enabled;
	}
	
	/**
	 * Indicate if this classifier must be used to group/classify ESTs.
	 * 
	 * @param enabled If the parameter is true, then this classifier entry
	 * will be used for grouping/classifying ESTs.
	 */
	public void setEnabled(boolean enabled) {
		this.enabled = enabled;
	}
	
	/**
	 * This method must be used to obtain the pre-compiled pattern 
	 * associated with this classifier.
	 * 
	 * This method creates and returns patterns in a "lazy" manner. A
	 * pattern is created the first time this method is called (after
	 * the regular expression is changed) and subsequent calls simply
	 * return the same pattern.
	 * 
	 * @return The pre-compiled pattern to be returned back to the caller.
	 */
	public Pattern getPattern() {
		if (pattern == null) {
			pattern = Pattern.compile(regExp);
		}
		return pattern;
	}
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent ClassifierList node in the 
	 * DOM tree.
	 * 
	 * @param dbClassList The DOM element corresponding to the "DBClassifier"
	 * node that contains this entry.
	 */
	public final void marshall(Element dbClassList) {
		// Create a top-level entry for this "DBClassifier"
		Element dbClas = DOMHelper.addElement(dbClassList, "DBClassifier", null);
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(dbClas, "RegExp", regExp);
		DOMHelper.addElement(dbClas, "Description", description);
		DOMHelper.addElement(dbClas, "Color", Integer.toString(color));
		DOMHelper.addElement(dbClas, "Enabled", enabled ? "true" : "false");
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t";
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		// Create a top-level server entry for this server
		out.printf("%s<DBClassifier>\n", Indent); 
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "Description", description);
		out.printf(STR_ELEMENT, "RegExp", regExp);
		out.printf(STR_ELEMENT, "Color", Integer.toString(color));
		out.printf(STR_ELEMENT, "Enabled", enabled ? "true" : "false");
		// Close the job tag
		out.printf("%s</DBClassifier>\n", Indent);
	}
	
	/**
	 * A user defined description for this DB classifier enttry. This 
	 * description is set when a new classifier entry is created by the
	 * user via the GUI. The information is persisted in the work space
	 * configuration for future references.
	 */
	private String description;
	
	/**
	 * The regular expression that is used to identify a given database
	 * name. This information is set when a new classifier entry is 
	 * created by the user via the GUI. The information is persisted in
	 * the work space configuration for future references.
	 */
	private String regExp;
	
	/**
	 * The RGB color code associated with this entry. The color code is
	 * marshalled and unmarshalled into hexadecimal codes corresponding
	 * to RGB values.
	 */
	private int color;
	
	/**
	 * Flag to indicate if this classifier is enabled and must be used to
	 * classify ESTs. If a classifier is not enabled, then it is not used
	 * for classifying ESTs.
	 */
	private boolean enabled;
	
	/**
	 * The pre-compiled pattern corresponding to this regular expression.
	 *
	 * The pre-compiled pattern is created in a "lazy" manner the first
	 * time the getPattern() method is called. The pattern is reset each
	 * time the setRegExp() method is invoked.
	 */
	transient private Pattern pattern;
}
