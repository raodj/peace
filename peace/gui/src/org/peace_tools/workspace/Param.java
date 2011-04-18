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

import org.w3c.dom.Element;

/** 
 * A simple class to encapsulate a <Name, Value> pair that is used to
 * describe parameters. These parameters are are supplied to heuristics
 * and filters to fine tune their operations. This class is primarily
 * designed to be used in the Heuristic and Job classes.
 */
public class Param {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * parameter entry. This method is typically  used to create a suitable
	 * Job entry when loading a Work space into the GUI.
	 * 
	 * @param paramNode The DOM element to be used for creating the 
	 * parameter entry and populating with the needed data.
	 * 
	 * @return The newly created parameter entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static Param create(Element paramNode) throws Exception {
		// First extract the name information from the DOM tree.
		String name  = DOMHelper.getStringValue(paramNode, "Name");
		// Extract value if we have one.
		String value = null;
		if (DOMHelper.hasElement(paramNode, "Value")) {
			value = DOMHelper.getStringValue(paramNode, "Value");
		}
		return new Param(name, value);
	}
	
	/**
	 * The constructor merely initializes the <Name, Value> pair 
	 * encapsulated by this object.
	 * 
	 * @param name The name of the parameter.
	 * @param value The value set for the parameter.
	 */
	public Param(String name, String value) {
		this.name  = name;
		this.value = value;
	}

	/**
	 * The name set for this parameter.
	 * 
	 * @return This method returns the name set for this parameter.
	 */
	public String getName() { return name; }

	/**
	 * The value associated with this parameter.
	 * 
	 * @return The value associated with this parameter.
	 */
	public String getValue() { return value; }
	
	/**
	 * Method to marshal the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to an appropriate node in the DOM tree.
	 * 
	 * @param parent The DOM element corresponding to the parent node
	 * that contains this (and other) parameter nodes.
	 */
	public final void marshall(Element parent) {
		// Create a top-level entry for this "Job"
		Element param = DOMHelper.addElement(parent, "Param", null);
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(param, "Name", name);
		if (value != null) {
			DOMHelper.addElement(param, "Value", value);
		}
	}

	/**
	 * Method to marshal the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t\t";
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		// Create a top-level server entry for this server
		out.printf("%s<Param>\n", Indent); 
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "Name", name);
		if (value != null) {
			out.printf(STR_ELEMENT, "Value", value);
		}
		out.printf("%s</Param>\n", Indent);
	}
	
	@Override
	public String toString() {
		return name + (value != null ? " " + value : "");
	}
	
	/**
	 * The name set for this parameter. This value is set when this parameter
	 * object is instantiated.
	 */
	private final String name;

	/**
	 * The value set for this parameter. This value is set when this parameter
	 * object is instantiated.
	 */
	private final String value;
}
