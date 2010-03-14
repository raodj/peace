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
import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * This class is meant to encapsulate the information and parameters for various
 * filters used to filter out ESTs that may negatively impact clustering. Similar 
 * to many of the other classes, this class provides a convenient interface to 
 * marshall, unmarshall, and use work space configuration information. This 
 * class is a generic Filter class that is used to store information regarding 
 * all the filters used within a Job element.
 * 
 */
public class Filter {
	/**
	 * Different enumerations defining the different types of Filters that are
	 * currently supported. Note that the enumeration of filter types is identical
	 * to the names in the schema file to ease marshalling to XML.
	 */
	public enum FilterType {
		/**
		 * This filter is used to filter out ESTs that are shorter than a 
		 * given length.
		 */
		LengthFilter,
		/**
		 * This filter is used to filter out ESTs that have a Low Complexity
		 * section in them.
		 */
		LCFilter,
	};
	
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * Filter entry. This method is typically  used to create a suitable
	 * entry when loading a Work space into the GUI.
	 * 
	 * @param filter The DOM element to be used for creating the
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created filter entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static Filter create(Element filter) throws Exception {
		// First extract the necessary information from the DOM tree.
		String name = DOMHelper.getStringValue(filter, "Name");
		// Convert the filter name to suitable enumeration type.
		FilterType type = FilterType.valueOf(FilterType.class, name);
		
		// Create the core filter class and then add parameters directly
		Filter node = new Filter(type);
		
		// Process list of parameters for the filter.
		NodeList params = filter.getElementsByTagName("Param");
		// Process each parameter in the list
		for(int idx = 0; (idx < params.getLength()); idx++) {
			Element param    = (Element) params.item(idx);
			String paramName = DOMHelper.getStringValue(param, "Name");
			String paramValue= DOMHelper.getStringValue(param, "Value");
			// Add the parameter information to the parameters list.
			node.parameters.add(new Param(paramName, paramValue));
		}
		// Return the newly created object.
		return node;
	}
	
	/**
	 * The only constructor for this class.
	 * 
	 * @param name The name of the heuristic with which this class is
	 * associated.
	 */
	public Filter(FilterType filterType) {
		this.filterType = filterType;
		this.parameters = new ArrayList<Param> ();
	}
	
	/**
	 * Obtain the type of filter represented by this class.
	 * 
	 * This method provides the name in a form that can be readily 
	 * dispatched to the PEACE clustering engine.
	 * 
	 * @return This method returns the type of filter represented by this 
	 * class.
	 */
	public String getName() {
		final String FilterNames[] = {"lengthFilter", "lcFilter"};
		return FilterNames[filterType.ordinal()]; 
	}
	
	/**
	 * Obtain the full list of parameters passed to this filter.
	 * 
	 * @return This method returns the full list of parameters supplied to
	 * this heuristic to fine tune its runtime operations.
	 */
	public ArrayList<Param> getParameters() { return parameters; }
	
	/**
	 * Add a parameter to the filter.
	 * 
	 * This method may be used to add a parameter to the filter.
	 * The name of the parameter must be one of the command line 
	 * parameters in PEACE. The value must be suitably set to be
	 * compatible with the parameter. This method does not perform
	 * any special checks (at least not yet, maybe it should be)
	 * to verify that the parameter is valid. 
	 * 
	 * @param p The parameter to be added to the list. This value
	 * must not be null.
	 */
	public void addParameter(Param p) {
		parameters.add(p);
	}
	
	/**
	 * Provides a one line information about this heuristics.
	 * Parameters are included on the line as comma separated list
	 * of values.
	 */
	@Override
	public String toString() {
		String paramInfo = "";
		for(Param p: parameters) {
			// Add a comma to separate parameters as needed.
			paramInfo += (paramInfo.length() > 0) ? ", " : "";
			// Add parameter name ":" value information.
			paramInfo += p.getName() + ":" + p.getValue();
		}
		return filterType.toString() + "[Parameters: " + paramInfo + "]";
	}
	
	/**
	 * Provides a multi-line information about this filter.
	 * Parameters are displayed on separate lines with each line indented
	 * by a single tab. This is usually used to display summary information
	 * to the user.
	 */
	public String getSummary() {
		String summary = "Filter: " + filterType + "\n";
		for(Param p: parameters) {
			// Add parameter name ":" value information.
			summary += "\t" + p.getName() + ":" + p.getValue() + "\n"; 
		}
		return summary;
	}
	
	/**
	 * Return the information in the form of a partial PEACE command line 
	 * parameter.
	 * 
	 * This method can be used to obtain the information needed to
	 * generate MST data while applying this heuristic. 
	 * 
	 * @return Return the information as a command line parameter.
	 */
	public String toCmdLine() {
		String cmdLine = "";
		for(Param p: parameters) {
			cmdLine += " --" + p.getName() + " " + p.getValue();
		}
		return cmdLine;
	}
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent FilterList node in the DOM tree.
	 * 
	 * @param filterList The DOM element corresponding to the "FilterList"
	 * node that contains this entry.
	 */
	public final void marshall(Element filterList) {
		// Create a top-level entry for this "Job"
		Element filter = DOMHelper.addElement(filterList, "Filter", null);
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(filter, "Name", filterType.toString());
		for(Param param : parameters) {
			Element paramNode = DOMHelper.addElement(filter, "Param", null);
			DOMHelper.addElement(paramNode, "Name",  param.getName());
			DOMHelper.addElement(paramNode, "Value", param.getValue());
		}
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t\t\t";
		final String STR_ELEMENT = Indent + "\t\t" + "<%1$s>%2$s</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<Filter>\n", Indent);
		out.printf("%s\t<Name>%s</Name>\n", Indent, filterType.toString());
		// Add new sub-elements for each parameter
		for(Param param : parameters) {
			out.printf("%s\t<Param>\n", Indent);
			out.printf(STR_ELEMENT, "Name",  param.getName());
			out.printf(STR_ELEMENT, "Value", DOMHelper.xmlEncode(param.getValue()));
			out.printf("%s\t</Param>\n", Indent);
		}
		// Close the heuristic tag
		out.printf("%s</Filter>\n", Indent);
	}
	
	/**
	 * The name/type of this filter. This value is set when this class
	 * is instantiated and is never changed.
	 */
	private final FilterType filterType;
	
	/**
	 * The list of parameters that are simply managed as a name-value
	 * pair within this class.
	 */
	private ArrayList<Param> parameters;
}
