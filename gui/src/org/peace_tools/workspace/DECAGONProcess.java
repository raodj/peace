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
 * <p>This class corresponds to a "DECAGONProcess" element in a 
 * PEACE work space XML data. This class encapsulates the core 
 * information associated with a DECAGON process. A non-null set of
 * DECAGON-processes constitute a DECAGON job. This class serves 
 * merely as a read-only type class that is created when a new 
 * job is run. The data is persisted in the XML work space for 
 * future reference so that users can determine jobs that were 
 * scheduled and check on status of long-running jobs. Currently,
 * this class is rather simple and contains just the full comamnd-line
 * used to run the process. However,  
 * 
 * <p>In addition to encapsulating the data this class also provides
 * convenient interfaces for marshaling and un-marshaling XML data compatible
 * with the PEACE GUI configuration XML file format. </p>
 */
public class DECAGONProcess {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * process entry. This method is typically  used to create a suitable
	 * process entry when loading a Work space into the GUI.
	 * 
	 * @param procNode The DOM element to be used for creating the process
	 * entry and populating with the needed data. This node must
	 * correspond to "DECAGONProcessType" element.
	 * 
	 * @return The newly created process entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static DECAGONProcess create(Element procNode) throws Exception {
		// First extract the necessary information to create a 
		// process node.
		final String cmdLine = DOMHelper.getStringValue(procNode, "CmdLine");
		
		// Now we have the information to create a minimal decagon
		// job class. So do it and use it for further operations.
		DECAGONProcess process = new DECAGONProcess(cmdLine);
		
		// Return the fully populated entry back to the caller
		return process;
	}
		
	/**
	 * Constructor to create a DECAGON process object.
	 *  
	 * This constructor is typically used to create a DECAGON job
	 * object that is part of a DECAGON job. This constructor
	 * is used in the {@link #create(Element)} method.
	 * 
	 * @param cmdLine The command-line that was used to run an 
	 * instance of the process.
	 */
	public DECAGONProcess(String cmdLine) {
		// Save information about heuristics and filters
		this.cmdLine = cmdLine;
	}
		
	/**
	 * Obtain the command line used for running this process.
	 * 
	 * This method can be used to obtain the command-line (stored in this
	 * object) for running this information needed to run the
	 * job based on the supplied information in the form of a command line. 
	 * 
	 * @return Return the information as a command line that was used to
	 * run this job.
	 */
	public String getCmdLine() {
		return cmdLine;
	}
	
	@Override
	public String toString() {
		return "process (" + cmdLine + ")";
	}
	
	/**
	 * Method to marshal the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent DECAGON job node in the DOM tree.
	 * 
	 * @param decJob The DECAGON job element corresponding to the 
	 * "DECAGONJob" node that contains this entry.
	 */
	public final void marshall(Element decJob) {
		// Create a top-level entry for this "Job"
		Element proc = DOMHelper.addElement(decJob, "Process", null);
		DOMHelper.addElement(proc, "CmdLine", cmdLine);
	}
	
	/**
	 * Method to marshal the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 * 
	 * @param indentPrefix The extra indentation to be done to make the
	 * output look nice. If no additional indentation is needed then
	 * an empty string ("") must be passed in.
	 */
	public final void marshall(PrintWriter out, final String indentPrefix) {
		final String Indent = indentPrefix + "\t";
		// Create a top-level server entry for this server
		out.printf("%s<Process>\n", indentPrefix);
		// Now marshal the information regarding the process.
		out.printf("%s<CmdLine>%s</CmdLine>\n", Indent, 
				DOMHelper.xmlEncode(cmdLine));
		// Close the process tag
		out.printf("%s</Process>\n", indentPrefix);
	}
		
	/**
	 * This element contains the full, final command-line that was used 
	 * to run this process. This command-line is loaded from the PEACE
	 * workspace configuration file. The list of parameters is originally
	 * constructed by DECAGON from the DADX file description for this 
	 * process.
	 */
	private String cmdLine;
}
