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
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * A class to encapsulate information about a list of DBClassifier entries
 * that have already been configured in this work space. This class is 
 * instantiated from the Work space class. This class is relatively 
 * straightforward in that it merely contains a list of DBClassifier entires. 
 * In addition, it facilitates marshalling and unmarshalling of data and
 * handling generation of notification events on change.
 */
public class ClassifierList {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * ClassifierList entry. This method is typically  used to create a suitable
	 * ClassifierList entry using data from a given DOM element. For each 
	 * DBClassifier entry this method uses the DBClassifier.create method to
	 * create suitable entries.
	 * 
	 * @param classifierListNode The DOM element to be used for creating the 
	 * list and to be used for creating the DBClassifier entries.
	 * 
	 * @return The newly created ClassifierList entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static ClassifierList create(Element dbClassListNode) throws Exception {
		// Obtain all the classifier elements into an array list.
		ArrayList<DBClassifier> list = new ArrayList<DBClassifier>();
		NodeList classifiers = dbClassListNode.getElementsByTagName("DBClassifier");
		for(int idx = 0; (idx < classifiers.getLength()); idx++) {
			Node tmpNode = classifiers.item(idx);
			// Ensure that this node is actually a classifier node
			if ((tmpNode.getNodeType() == Node.ELEMENT_NODE) &&
				("DBClassifier".equals(tmpNode.getNodeName())) && 
				(tmpNode instanceof Element)) {
				// OK, safe to type cast and try to parse..
				Element dbNode = (Element) tmpNode;
				DBClassifier classifier = DBClassifier.create(dbNode);
				// Add the valid classifier to the temporary list.
				list.add(classifier);
			}
		}
		// Create and return a suitable classifier list object.
		ClassifierList dbClassList = new ClassifierList();
		dbClassList.classifiers    = list;
		return dbClassList;
	}
	
	/**
	 * The default constructor. It merely initializes all the instance
	 * variables to their default initial value. This creates an empty
	 * DBClassifier list.
	 */
	public ClassifierList() {
		classifiers = new ArrayList<DBClassifier>();
	}

	/**
	 * Method to set a new set of classifier entries.
	 * 
	 * This method sets a new set of classifiers to be used for classifying
	 * ESTs based on their data base source information. If the supplied
	 * list is empty then all classifiers are cleared.
	 * 
	 * @param classifierList The new list of classifier entries to be set for
	 * classifiying EST entries based on their origin data base names.
	 */
	public void set(final ArrayList<DBClassifier> classifierList) {
		// First make a copy for ourselves.
		classifiers = new ArrayList<DBClassifier>(classifierList);
		// Fire notification to listeners to update GUIs
		WorkspaceEvent we = new WorkspaceEvent(this);
		Workspace.get().fireWorkspaceChanged(we);
	}
	
	/**
	 * Obtain reference to a <b>clone</b> of list of classifiers. 
	 * 
	 * This method must be used to obtain the list of classifiers currently
	 * configured for this work space. This method makes a copy of the 
	 * classifiers and returns the copy. Consequently, the actual classifier
	 * list cannot be directly modified. Instead, the set() method must be
	 * used to actually udpate the list of entries.
	 *
	 * @note The list of classifiers returned by this method is a clone and
	 * modifying them does not modify the classifiers in this class.
	 * 
	 * @param classifierID The generated work space wide unique classifier ID. 
	 * Note that checks are case sensitive. 
	 * 
	 * @return A copy of the list of classifiers associated with this work space.
	 */
	public ArrayList<DBClassifier> getClassifiers() {
		return new ArrayList<DBClassifier>(classifiers);
	}
	
	/**
	 * Convenience method to determine the number of classifiers in this list.
	 * 
	 * This method returns the number of classifiers currently configured in
	 * this list.
	 * 
	 * @return The number of classifiers currently in the classifier list.
	 */
	public int getSize() {
		return classifiers.size();
	}
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent Workspace node in the DOM tree.
	 * 
	 * @param workspace The DOM element corresponding to the Workspace
	 * node that contains this entry.
	 */
	public final void marshall(Element workspace) {
		// Create a top-level classifier list entry for this class
		Element classifierList = DOMHelper.addElement(workspace, "ClassifierList", null);
		// Add sub-elements for each classifier in our classifier list
		for (DBClassifier classifier : classifiers) {
			classifier.marshall(classifierList);
		}
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragement. The XML fragement is guaranteed to be compatible
	 * with the PEACE work space configuration data. This method provides
	 * better control on the XML formatting to generate a more readable
	 * XML output when compared to the DOM tree.
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t";
		// Create a top-level classifier list entry for this class
		out.printf("%s<ClassifierList>\n", Indent);
		// Add sub-elements for each classifier in our classifier list
		for (DBClassifier srvr : classifiers) {
			srvr.marshall(out);
		}
		// Close the classifier list element.
		out.printf("%s</ClassifierList>\n\n", Indent);
	}
	
	/**
	 * The list of DB classifier objects that have been configured and added
	 * to this list. This list is populated either via the create 
	 * method or new entries are set via the setClassifiers() method
	 * in this class. 
	 */
	private ArrayList<DBClassifier> classifiers;
}
