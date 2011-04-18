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

import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * This class contains a collection of static utility method that are used to
 * retrieve or create elements in a given DOM tree. All the methods in this 
 * class are static methods and can be directly invoked. This class is not
 * meant to be (and cannot be) instantiated. 
 */
public class DOMHelper {
	/**
	 * This string defines the name space for the XML elements created by the
	 * utility methods.
	 */
	public static String PEACE_NS = "http://www.peace-tools.org/";
	
	/**
	 * Utility method to obtain the value of a given element as a string.
	 * 
	 * This method can be used to search for a given child-element within the
	 * root document and if the element is found, it returns the value of the
	 * child element as a String. 
	 * 
	 * <p><b>Note:</b> This method does not recursively search child elements.
	 * Only the immediate underlying child elements are searched.</p>
	 * 
	 * @param domDoc The root of the DOM document tree whose immediate underlying
	 * child elements are to be searched.
	 * @param elementName The name of the element whose value is to be converted
	 * to an integer and returned.
	 * 
	 * @return The value of the element. If the element was not found, then this
	 * method will raise an exception. So prior to calling this method
	 * check to ensure the element exists via a call to DOMHelper.hasElement() 
	 */
	public static String getStringValue(Document domDoc, String elementName) {
		NodeList nodes = domDoc.getElementsByTagName(elementName);
		return getStringValue(nodes, elementName, false);
	}

	/**
	 * Utility method to obtain the value of a given element as a string.
	 * 
	 * This method can be used to search for a given child-element within another
	 * element and if the element is found, it returns the value of the
	 * child element as a string. 
	 * 
	 * <p><b>Note:</b>  This method does not recursively search its 
	 * child elements. Only the immediate underlying child elements
	 * are searched.</p>
	 * 
	 * @param parent The parent element whose immediate child elements are to 
	 * be searched.
	 * @param elementName The name of the element whose value is to be converted
	 * to an integer and returned.
	 * 
	 * @return The value of the element. If the element was not found, then this
	 * method will raise an exception. So prior to calling this method
	 * check to ensure the element exists via a call to DOMHelper.hasElement() 
	 */
	public static String getStringValue(Element parent, String elementName) {
		NodeList nodes = parent.getElementsByTagName(elementName);
		return getStringValue(nodes, elementName, false);
	}
	
	/**
	 * Utility method to obtain the value of a given element as a string.
	 * 
	 * This method can be used to search for a given child-element within another
	 * element and if the element is found, it returns the value of the
	 * child element as a string. 
	 * 
	 * <p><b>Note:</b> This method does not recursively search its child elements. Only the
	 * immediate underlying child elements are searched.</p>
	 * 
	 * @param parent The parent element whose immediate child elements are to 
	 * be searched.
	 * 
	 * @param elementName The name of the element whose value is to be converted
	 * to an integer and returned.
	 * 
	 * @param canBeNull This flag indicates that the data can be null.
	 * 
	 * @return The value of the element. If the element was not found, then this
	 * method will raise an exception. So prior to calling this method
	 * check to ensure the element exists via a call to DOMHelper.hasElement() 
	 */
	public static String getStringValue(Element parent, String elementName, 
			boolean canBeNull) {
		NodeList nodes = parent.getElementsByTagName(elementName);
		return getStringValue(nodes, elementName, canBeNull);
	}
	
	/**
	 * Utility method to obtain the value of a given element as a string or
	 * use a default value.
	 * 
	 * This method can be used to search for a given child-element within another
	 * element and if the element is found, it returns the value of the
	 * child element as a string. Otherwise this method returns the supplied
	 * default value. 
	 * 
	 * <p><b>Note:</b> This method does not recursively search its child elements.
	 * Only the immediate underlying child elements are searched.</p>
	 * 
	 * @param parent The parent element whose immediate child elements are to 
	 * be searched.
	 * 
	 * @param elementName The name of the element whose value is to be converted
	 * to an integer and returned.
	 * 
	 * @param valueIfNull This value is returned in case the element does not
	 * have a value specified.
	 * 
	 * @return The value of the element if it is not null. If the element has
	 * a null value then the supplied valueIfNull object is returned. If the 
	 * element was not found, then this method will raise an exception. So prior
	 * to calling this method  check to ensure the element exists via a call to 
	 * {@link #hasElement(Element, String)}. 
	 */
	public static String getStringValue(Element parent, String elementName, 
			String valueIfNull) {
		NodeList nodes = parent.getElementsByTagName(elementName);
		String retVal = getStringValue(nodes, elementName, true);
		return (retVal != null ? retVal : valueIfNull);
	}
	
	/**
	 * Utility method to obtain the value of a given element as a String.
	 * 
	 * This method can be used to search for a given child-element within a
	 * list of nodes and if the element is found, it returns the value of the
	 * child element as a String. No interpretations are performed on the String
	 * value returned by this method.
	 * 
	 * <p><b>Note:</b>  This method does not recursively search the list of nodes. It
	 * only checks the given list of nodes.</p>
	 * 
	 * @param nodes The list of nodes to be searched.
	 * @param elementName The name of the element whose value is to be returned
	 * as a String.
	 * @param canBeNull This flag indicates that the data can be null.
	 * 
	 * @return The value of the element. If the element was not found, then this
	 * method will raise an exception. So prior to calling this method
	 * check to ensure the element exists via a call to DOMHelper.hasElement() 
	 */
	public static String getStringValue(NodeList nodes, String elementName, 
			boolean canBeNull) {
		if ((nodes == null) || (nodes.getLength() == 0)) {
			throw new DOMException(DOMException.NOT_FOUND_ERR,
					"Node with name '" + elementName + "' not found.");
		}
		// Get the first node off the list.
		Node node = nodes.item(0);
		// Check to ensure that it has a first child and it is a text node.
		if ((node.getFirstChild() == null)
				|| (node.getFirstChild().getNodeType() != Node.TEXT_NODE)) {
			if (canBeNull) {
				return null;
			}
			throw new DOMException(DOMException.TYPE_MISMATCH_ERR,
					"Node with name '" + elementName +
					"' did not have data associated with it.");
		}
		// The first child is a text node as expected. Return its value.
		return node.getFirstChild().getNodeValue();
	}

	/**
	 * Utility method to obtain the value of a given element as an integer.
	 * 
	 * This method can be used to search for a given child-element within another
	 * element and if the element is found, it returns the value of the
	 * child element as an integer. This method assumes that the value of the
	 * element is expected to be an integer and it has passed validation.
	 * 
	 * <p><b>Note:</b>  This method does not recursively search its child elements. Only the
	 * immediate underlying child elements are searched.</p>
	 * 
	 * @param parent The parent element whose immediate child elements are to 
	 * be searched.
	 * @param elementName The name of the element whose value is to be converted
	 * to an integer and returned.
	 * @return The value of the element. If the element was not found, then this
	 * method will raise an exception. So prior to calling this method
	 * check to ensure the element exists via a call to DOMHelper.hasElement() 
	 */
	public static int getIntValue(Element parent, String elementName) {
		NodeList nodes = parent.getElementsByTagName(elementName);
		String str = getStringValue(nodes, elementName, false);
		return Integer.parseInt(str);
	}


	/**
	 * Utility method to obtain the value of a given element as a long value.
	 * 
	 * This method can be used to search for a given child-element within another
	 * element and if the element is found, it returns the value of the
	 * child element as an integer. This method assumes that the value of the
	 * element is expected to be an integer and it has passed validation.
	 * 
	 * <p><b>Note:</b>  This method does not recursively search its child elements. Only the
	 * immediate underlying child elements are searched.</p>
	 * 
	 * @param parent The parent element whose immediate child elements are to 
	 * be searched.
	 * @param elementName The name of the element whose value is to be converted
	 * to an integer and returned.
	 * @return The value of the element. If the element was not found, then this
	 * method will raise an exception. So prior to calling this method
	 * check to ensure the element exists via a call to DOMHelper.hasElement() 
	 */
	public static long getLongValue(Element parent, String elementName) {
		NodeList nodes = parent.getElementsByTagName(elementName);
		String str = getStringValue(nodes, elementName, false);
		return Long.parseLong(str);
	}

	/**
	 * Utility method to obtain the value of a given element as a double.
	 * 
	 * This method can be used to search for a given child-element within another
	 * element and if the element is found, it returns the value of the
	 * child element as an integer. This method assumes that the value of the
	 * element is expected to be an integer and it has passed validation.
	 * 
	 * <p><b>Note:</b>  This method does not recursively search its child elements. Only the
	 * immediate underlying child elements are searched.</p>
	 * 
	 * @param parent The parent element whose immediate child elements are to 
	 * be searched.
	 * @param elementName The name of the element whose value is to be converted
	 * to an integer and returned.
	 * @return The value of the element. If the element was not found, then this
	 * method will raise an exception. So prior to calling this method
	 * check to ensure the element exists via a call to DOMHelper.hasElement() 
	 */
	public static double getDoubleValue(Element parent, String elementName) {
		NodeList nodes = parent.getElementsByTagName(elementName);
		String str = getStringValue(nodes, elementName, false);
		return Double.parseDouble(str);
	}
	
	/**
	 * Utility method to determine if a given DOM element has an optional sub-element.
	 * 
	 * This is a utility method that can be used to determine if  an immediate
	 * underlying child element is present within a DOM document.
	 * 
	 * @param parent The DOM element within which to search.
	 * 
	 * @param elementName The name of element to be located and returned.
	 * 
	 * @return This method returns true if the DOM element exists. Otherwise
	 * this method returns false.
	 */
	public static boolean hasElement(Element parent, String elementName) {
		NodeList nodes = parent.getElementsByTagName(elementName);
		if ((nodes == null) || (nodes.getLength() < 1)) {
			return false;
		}
		// Get the first node off the list.
		Node node = nodes.item(0);
		// Check to ensure that it has a first child and it is an element node.
		if ((node == null) || (node.getNodeName() == null) ||
			(!elementName.equals(node.getNodeName())) || 
			(node.getNodeType() != Node.ELEMENT_NODE) ||
			(!(node instanceof Element))) {
			return false;
		}
		// This is the element we are looking for. We found it.
		return true;
	}
	
	/**
	 * Utility method to search for a given immediate child element.
	 * 
	 * This is a utility method that can be used to search for an immediate
	 * underlying child element within a DOM document.
	 * 
	 * @param parent The root DOM document within which to search.
	 * @param childName The name of element to be located and returned.
	 * @return This method returns a valid DOM element object if an element with
	 * the given (childName) was found. If the element was not found then this
	 * method returns null.
	 */
	public static Element getElement(Document parent, String childName) {
		return getElement(parent.getElementsByTagName(childName), childName);
	}

	/**
	 * Utility method to search for a given immediate child element.
	 * 
	 * This is a utility method that can be used to search for an immediate
	 * underlying child element within a DOM element.
	 * 
	 * @param parent The parent DOM element within which to search.
	 * @param childName The name of element to be located and returned.
	 * @return This method returns a valid DOM element object if an element with
	 * the given (childName) was found. If the element was not found then this
	 * method returns null.
	 */
	public static Element getElement(Element parent, String childName) {
		return getElement(parent.getElementsByTagName(childName), childName);
	}
	
	/**
	 * Helper method to search within a given list of nodes to locate a given
	 * element.
	 * 
	 * This method is an internal helper method that is called from other
	 * overloaded getElement() methods to locate, test, and return the element
	 * corresponding to a given element name.
	 * 
	 * <p>
	 * <b>Note:</b> This method does not search recursively within the list of
	 * nodes. Only the immediately underlying child elements are searched for a
	 * match.
	 * </p>
	 * 
	 * @param nodes
	 *            The list of DOM nodes in which to search for a given element.
	 * @param elementName
	 *            The name of element to search for.
	 * @return The element node within the list of DOM nodes. If the element is
	 *         not found, this method returns null.
	 */
	private static Element getElement(NodeList nodes, String elementName) {
		if ((nodes == null) || (nodes.getLength() < 1)) {
			throw new DOMException(DOMException.NOT_FOUND_ERR,
					"Node with name '" + elementName + "' not found.");
		}
		// Get the first node off the list.
		Node node = nodes.item(0);
		// Check to ensure that it has a first child and it is a text node.
		if ((node == null) || (node.getNodeName() == null) ||
			(!elementName.equals(node.getNodeName())) || 
			(node.getNodeType() != Node.ELEMENT_NODE) ||
			(!(node instanceof Element))) {
			throw new DOMException(DOMException.NOT_FOUND_ERR,
					"Node with name '" + elementName + "' not found.");
		}
		// This is the element we are looking for. Return it.
		return (Element) node;
	}
	
	/**
	 * Helper method to add a new element to a given DOM node.
	 * This method is a helper method that is used to add a simple
	 * element (under the PEACE name space) with a simple text value
	 * to a given DOM element.
	 * 
	 * @param parent The DOM element to which a new element is to be added.
	 * @param name The name of the new element to be added to the DOM tree.
	 * @param value The value to be associated with the new element.
	 * @return If the element is added successfully, then this method returns
	 * the newly created element. On errors it returns null. 
	 */
	public static Element addElement(Element parent, String name,
			String value) {
		Document document = parent.getOwnerDocument();
		if (document == null) {
			return null;
		}
		// Create a new element for the child node.
		Element child = document.createElementNS(PEACE_NS, name);
		// Set its value via a suitable text node 
		if (value != null) {
			Node valueNode = document.createTextNode(value);
			child.appendChild(valueNode);
		}
		// Now add the child node to the parent.
		parent.appendChild(child);
		// Return the newly created child node.
		return child;
	}
	
	/**
	 * Helper method to add a new element to a given DOM node.
	 * This method is a helper method that is used to add a simple
	 * element (under the PEACE name space) with a simple integer value
	 * to a given DOM element.
	 * 
	 * @param parent The DOM element to which a new element is to be added.
	 * @param name The name of the new element to be added to the DOM tree.
	 * @param value The value to be associated with the new element.
	 * @return If the element is added successfully, then this method returns
	 * the newly created element. On errors it returns null. 
	 */
	public static Element addElement(Element parent, String name,
			long value) {
		return addElement(parent, name, "" + value);
	}
	
	/**
	 * Helper method to add a new element to a given DOM node.
	 * This method is a helper method that is used to add a simple
	 * element (under the PEACE name space) with a simple float value
	 * to a given DOM element.
	 * 
	 * @param parent The DOM element to which a new element is to be added.
	 * @param name The name of the new element to be added to the DOM tree.
	 * @param value The value to be associated with the new element.
	 * @return If the element is added successfully, then this method returns
	 * the newly created element. On errors it returns null. 
	 */
	public static Element addElement(Element parent, String name,
			float value) {
		return addElement(parent, name, "" + value);
	}
	
    /** 
     * Encodes any text as PCDATA.
     *
     * This is a helper method that is used to encode a given string
     * into XML safe character set. Specifically, this method translates
     * the following characters: &amp;, &lt;, &gt;, &quot;, and &apos;. 
     * 
     * @param text The text to be encoded to XML. If this parameter is
     * null, this method immediately exits with null.
     * 
     * @return The string character data suitably encoded into a XML
     * safe string.
     */
	public final static String xmlEncode(String text) {
		if (text == null) {
			return null;
		}
		// Use helper strings for encoding.
		String StringChars = "&<>\"'";
		String XMLChars[] = {"&amp;", "&lt;", "&gt;", "&quot;", "&apos;"};
		// Encode characters into a string buffer.
		StringBuilder codedStr = new StringBuilder(text.length() * 2);
		for (int i = 0; (i < text.length()); i++) {
			char c    = text.charAt(i);
			int index;
			if ((index = StringChars.indexOf(c)) != -1) {
				codedStr.append(XMLChars[index]);
			} else {
				codedStr.append(c);
			}
		}
		// Return the encoded string.
		return codedStr.toString();
	}

	/**
	 * This class is not meant to be instantiated. Therefore the
	 * constructor has been made private to enforce this policy.
	 */
	private DOMHelper() {
		// Intentionally supressed constructor.
	}
}
