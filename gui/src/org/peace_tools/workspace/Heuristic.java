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
 * heuristics used to accelerate clustering. Similar to many of the other
 * classes, this class provides a convenient interface to marshall, unmarshall,
 * and use work space configuration information. This class is a generic
 * Heuristic class that is used to store information regarding all the
 * heuristics used within a MSTData element.
 * 
 */
public class Heuristic {
	/** 
	 * A simple class to encapsulate a <Name, Value> pair that is used to
	 * describe the various parameters that are supplied to heuristics
	 * to fine tune their operations.
	 */
	public class Param {
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
	
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * Heuristic entry. This method is typically  used to create a suitable
	 * entry when loading a Work space into the GUI.
	 * 
	 * @param heuristic The DOM element to be used for creating the
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created Heuristic entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static Heuristic create(Element heuristic) throws Exception {
		// First extract the necessary information from the DOM tree.
		String name = DOMHelper.getStringValue(heuristic, "Name");
		// Create the core heuristic class and then add parameters directly
		Heuristic node = new Heuristic(name);
		// Process list of parameters for the heuristic.
		NodeList params = heuristic.getElementsByTagName("Param");
		// Process each parameter in the list
		for(int idx = 0; (idx < params.getLength()); idx++) {
			Element param    = (Element) params.item(idx);
			String paramName = DOMHelper.getStringValue(param, "Name");
			String paramValue= DOMHelper.getStringValue(param, "Value");
			// Add the parameter information to the parameters list.
			node.parameters.add(node.new Param(paramName, paramValue));
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
	public Heuristic(String name) {
		this.name = name;
		this.parameters = new ArrayList<Param> ();
	}
	
	/**
	 * Obtain the name of the heuristic associated with this class.
	 * 
	 * @return This method returns the name of the heuristic associated
	 * with this class.
	 */
	public String getName() { return name; }
	
	/**
	 * Obtain the full list of parameters passed to this heuristic.
	 * 
	 * @return This method returns the full list of parmeters supplied to
	 * this heuristic to fine tune its runtime operations.
	 */
	ArrayList<Param> getParameters() { return parameters; }
	
	/**
	 * Add a parameter to the heuristic.
	 * 
	 * This method may be used to add a parameter to the heuristic.
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
		return name + "[Parameters: " + paramInfo + "]";
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
	 * passed in corresponds to the parent HeuristicList node in the DOM tree.
	 * 
	 * @param heurList The DOM element corresponding to the "HeuristicList"
	 * node that contains this entry.
	 */
	public final void marshall(Element heurList) {
		// Create a top-level entry for this "Job"
		Element heur = DOMHelper.addElement(heurList, "Heuristic", null);
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(heur, "Name", name);
		for(Param param : parameters) {
			Element paramNode = DOMHelper.addElement(heur, "Param", null);
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
		out.printf("%s<Heuristic>\n", Indent);
		out.printf("%s\t<Name>%s</Name>\n", Indent, name);
		// Add new sub-elements for each parameter
		for(Param param : parameters) {
			out.printf("%s\t<Param>\n", Indent);
			out.printf(STR_ELEMENT, "Name",  param.getName());
			out.printf(STR_ELEMENT, "Value", DOMHelper.xmlEncode(param.getValue()));
			out.printf("%s\t</Param>\n", Indent);
		}
		// Close the heuristic tag
		out.printf("%s</Heuristic>\n", Indent);
	}
	
	/**
	 * The name of the heuristic. This value is set when this class
	 * is instantiated and is never changed.
	 */
	private final String name;
	
	/**
	 * The list of parameters that are simply managed as a name-value
	 * pair within this class.
	 */
	private ArrayList<Param> parameters;
}