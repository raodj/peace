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
 * A simple class to encapsulate information about the frame/word
 * analyzer that was used to generate Minimum Spanning Tree (MST) data.
 * This class provides the necessary helper methods to marhsall and 
 * unmarshall the needed data from-and-to XML.
 */
public class FWAnalyzer {
	/**
	 * The different kinds of frame/word analyzers that are currently
	 * supported in PEACE. This enumeration provides a mechanism to
	 * consistently refer to one of the algorithms used for computing
	 * distance/similarity metric.
	 */
	public enum FWAnalyzerType {
		TWOPASSD2, D2, D2ZIM, CLU
	};

	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * FWAnalyzer entry. This method is typically  used to create a suitable
	 * analyzer entry when loading a Work space into the GUI.
	 * 
	 * @param analyzer The DOM element to be used for creating the Job
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created analyzer entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static FWAnalyzer create(Element analyzer) throws Exception {
		// First extract the necessary information from the DOM tree.
		String metric   = analyzer.getAttribute("metric");
		String typeStr  = DOMHelper.getStringValue(analyzer, "Type");
		int    frame    = DOMHelper.getIntValue(analyzer, "WindowSize");
		int    word     = DOMHelper.getIntValue(analyzer, "WordSize");
		String cache    = DOMHelper.getStringValue(analyzer, "CacheType");
		int   cacheSize= DOMHelper.getIntValue(analyzer, "CacheSize");
		// Translate the typeStr to an enumeration value.
		typeStr = typeStr.toUpperCase();
		FWAnalyzerType type = FWAnalyzerType.valueOf(FWAnalyzerType.class, typeStr);
		boolean isDistance  = "distance".equals(metric);
		// Now that we have sufficient information create the core
		// FWAnalyzer object.
		return new FWAnalyzer(type, isDistance, frame, word, cache, cacheSize); 
	}
	
	/**
	 * The constructor merely initializes the immutable instance variables to
	 * their values using appropriate parameters.
	 * 
	 * @param type
	 *            The type of f/w analyzer that this object contains information
	 *            about.
	 * @param isDistance If this f/w analyzer provides a distance metric value 
	 * then this flag must be true. Otherwise it is assumed that this f/w analyzer
	 * provides a similarity metric.
	 * @param windowSize
	 *            The frame/window size that was used to analyze EST sequences.
	 * @param wordSize
	 *            The word size that was used to analyze EST sequences. The word
	 *            size represents subdivisions with a frame/window.
	 * @param cacheType
	 *            The type of cache that was used to cache distance/similarity
	 *            metrics generated during MST construction.
	 * @param cacheSize
	 *            The size of each cache entry (per EST) in number of
	 *            distance/similarity metrics that are to be cached to avoid
	 *            having to recompute metrics.
	 */
	public FWAnalyzer(FWAnalyzerType type, boolean isDistance, 
			int windowSize, int wordSize,
			String cacheType, int cacheSize) {
		this.type       = type;
		this.isDistance = isDistance;
		this.windowSize = windowSize;
		this.wordSize   = wordSize;
		this.cacheType  = cacheType;
		this.cacheSize  = cacheSize;
	}
	
	/**
	 * Determine the type of f/w analyzer that this object contains information
	 * about.
	 * 
	 * @return The type of f/w analyzer that was set when this object was 
	 * created (or loaded from an XML configuration file).
	 */
	FWAnalyzerType getType() { return type; }
	
	/**
	 * Determine if this f/w analyzer uses a distance metric or a similarity
	 * metric.
	 * 
	 * @return This method returns true if the f/w analyzer uses a distance
	 * metric.
	 */
	boolean isDistance() { return isDistance; }
	
	/**
	 * The frame/window size that was used to analyze EST sequences.
	 * 
	 * @return This method returns the frame or window size that was
	 * set when this object was created (or loaded from an XML configuration
	 * file).
	 */
	int getFrameSize() { return windowSize; }

	/**
	 * The word size that was used to analyze EST sequences. The word size
	 * represents subdivisions with a frame/window.
	 * 
	 * @return This method returns the word size that was set when this object
	 *         was created (or loaded from an XML configuration file).
	 */
	int getWordSize() { return wordSize; }

	/**
	 * The type of cache that was used to cache distance/similarity metrics
	 * generated during MST construction.
	 * 
	 * @return This method returns the type of cache that was set when this object
	 *         was created (or loaded from an XML configuration file).
	 */
	String getCacheType() { return cacheType; }

	/**
	 * The size of each cache entry (per EST) in number of distance/similarity
	 * metrics that are to be cached to avoid having to recompute metrics.
	 * 
	 * @return This method returns the cache size that was set when this object
	 *         was created (or loaded from an XML configuration file).
	 */
	int getCacheSize() { return cacheSize; }

	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent MSTData node in the DOM tree.
	 * 
	 * @param mstData The DOM element corresponding to the "MSTData"
	 * node that contains this entry.
	 */
	public final void marshall(Element mstData) {
		// Create a top-level server entry for this server
		Element analyzer = DOMHelper.addElement(mstData, "FWAnalyzer", null);
		analyzer.setAttribute("metric", (this.isDistance ? "distance" : "similarity"));
		// Add new sub-elements for each value.
		DOMHelper.addElement(analyzer, "Type", type.toString().toLowerCase());
		DOMHelper.addElement(analyzer, "WindowSize", "" + windowSize);
		DOMHelper.addElement(analyzer, "WordSize", "" + wordSize);
		DOMHelper.addElement(analyzer, "CacheType", cacheType);
		DOMHelper.addElement(analyzer, "CacheSize", "" + cacheSize);
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t\t";
		// Printf format strings for string & numerical elements
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		final String NUM_ELEMENT = Indent + "\t" + "<%1$s>%2$d</%1$s>\n";
		
		// Create a top-level entry
		out.printf("%s<FWAnalyzer metric=\"%s\">\n", Indent,
				(isDistance() ? "distance" : "similarity")); 
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "Type", type.toString().toLowerCase());
		out.printf(NUM_ELEMENT, "WindowSize", windowSize);
		out.printf(NUM_ELEMENT, "WordSize", wordSize);
		out.printf(STR_ELEMENT, "CacheType", cacheType);
		out.printf(NUM_ELEMENT, "CacheSize", cacheSize);
		// Close the element
		out.printf("%s</FWAnalyzer>\n", Indent);
	}
	
	/**
	 * Return the information in the form of a partial PEACE command 
	 * line.
	 * 
	 * This method can be used to obtain the information associated
	 * with this analyzer as a partial PEACE command line.
	 * 
	 * @return Return the information as a command line parameter.
	 */
	public String toCmdLine() {
		final String[] AnalyzerName = {"twopassD2", "d2", "d2zim", "clu"};
		String cmdLine = "--analyzer " + AnalyzerName[type.ordinal()];
		cmdLine += " --frame " + windowSize;
		cmdLine += " --word "  + wordSize;
		cmdLine += " --cacheType " + cacheType;
		cmdLine += " --cache " + cacheSize;
		return cmdLine;
	}
	
	/**
	 * The size of the window or the frame size that was used for analyzing
	 * EST sequences.
	 */
	private final int windowSize;

	/**
	 * The word size that was used for analyzing EST sequences.
	 */
	private final int wordSize;
	
	/**
	 * The type of cache that was used to hold the data distance/similarity
	 * metric when the analysis was conducted.
	 * 
	 */
	private final String cacheType;
	
	/**
	 * The size of each cache entry (per EST) in number of distance/similarity
	 * metrics that are to be cached to avoid having to recompute metrics.
	 */
	private final int cacheSize;
	
	/**
	 * The type of f/w analyzer that this object contains information about.
	 */
	private final FWAnalyzerType type;
	
	/**
	 * Flag to indicate if this f/w analyzer provides distance metrics
	 * (rather than similarity) for identifying similar ESTs.
	 */
	private final boolean isDistance;
}
