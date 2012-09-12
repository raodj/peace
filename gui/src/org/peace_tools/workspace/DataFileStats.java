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

package org.peace_tools.workspace;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import org.peace_tools.core.SummaryWriter;
import org.w3c.dom.Element;

public class DataFileStats {
	/**
	 * Constructor to create a fully populated data file statistics object.
	 * 
	 * @param filePath The full path to the cDNA file about which this
	 * statistics object is being created.
	 * 
	 * @param lastUpdate The last update time stamp indicating when this
	 * file was last updated. This value is typically obtained from the
	 * {@link java.io.File#lastModified()} method. If this value is -1,
	 * then the constructor uses the specified {@link #filePath} to determine
	 * and set the last update time. 
	 *  
	 * @param count The total number of cDNA fragments in the file.
	 * @param avgLength The average length of the cDNA fragments in the data
	 * file.
	 * 
	 * @param lengthSD The standard deviation in the length of the cDNA 
	 * fragments in the data file.
	 * 
	 * @param minLength The minimum length of the cDNA fragments in the
	 * data file.
	 * 
	 * @param maxLength The maximum length of the cDNA fragments in the
	 * data file.
	 */
	public DataFileStats(final String filePath, long lastUpdate, int count, 
			float avgLength, float lengthSD, int minLength, int maxLength) {
		this.filePath   = filePath;
		if (lastUpdate == -1) {
			File temp = new File(filePath);
			lastUpdate = temp.lastModified();
		}
		this.lastUpdate = lastUpdate;
		this.count      = count;
		this.avgLength  = avgLength;
		this.lengthSD   = lengthSD;
		this.minLength  = minLength;
		this.maxLength  = maxLength;
	}
	
	public DataFileStats(final String filePath, long lastUpdate, 
			final List<DataFileStats> statsList) {
		this.filePath   = filePath;
		if (lastUpdate == -1) {
			File temp = new File(filePath);
			lastUpdate = temp.lastModified();
		}
		this.lastUpdate = lastUpdate;
		// Build aggregate statistics from the list of values supplied.
		int count       = 0;
		float avgLen    = 0, avgLenSD   = 0;
		int minLen      = 0, maxLen     = 0;
		float insErrSum = 0, delErrSum  = 0, subErrSum = 0;
		if ((statsList != null) && (statsList.size() > 0)) {
			minLen = statsList.get(0).getMinLength();
			maxLen = statsList.get(0).getMaxLength();
			for(DataFileStats dfs: statsList) {
				count      += dfs.getCount();
				avgLen     += dfs.getAvgLength();
				avgLenSD   += dfs.getLengthSD();
				minLen      = Math.min(minLen, dfs.getMinLength());
				maxLen      = Math.max(maxLen, dfs.getMaxLength());
				insErrSum  += dfs.getInsErrorRate();
				delErrSum  += dfs.getDelErrorRate();
				subErrSum  += dfs.getSubErrorRate();
			}
		}
		// Compute aggregate averages while handling case when no stats 
		// were passed in.
		final int statsCount = ((statsList == null) || (statsList.size() == 0) ? 1 : statsList.size());
		this.count        = count;
		this.avgLength    = avgLen   / statsCount;
		this.lengthSD     = avgLenSD / statsCount;
		this.minLength    = minLen;
		this.maxLength    = maxLen;
		this.techType     = SeqTechType.UNKNOWN;
		this.insErrorRate = insErrSum / statsCount;
		this.delErrorRate = delErrSum / statsCount;
		this.subErrorRate = subErrSum / statsCount;
	}
	
	/**
	 * Set sequencing information for this statistics object.
	 * 
	 * This method is used to setup additional technology and associated
	 * information in this object.
	 * 
	 * @param techType The technology type associated with this object.
	 * 
	 * @param insErrorRate The insertion rate associated with this
	 * statistics object. If the rate is not known then this value 
	 * must be -1. The valid range of values is 0 to 1.0 (inclusive).
	 * 
	 * @param delErrorRate The deletion rate associated with this
	 * statistics object. If the rate is not known then this value 
	 * must be -1. The valid range of values is 0 to 1.0 (inclusive).
	 * 
	 * @param subErrorRate The substitution error rates associated with 
	 * the given technology type. If the rate is not known then this 
	 * value must be -1. The valid range of values is 0 to 1.0 (inclusive).
	 */
	public void setSeqInfo(final SeqTechType techType, 
			final float insErrorRate, final float delErrorRate,
			final float subErrorRate) {
		this.techType     = techType;
		this.insErrorRate = insErrorRate;
		this.delErrorRate = delErrorRate;
		this.subErrorRate = subErrorRate;
	}
	
	/**
	 * Convenience constructor to load data from a DOM Tree
	 * 
	 * This constructor is typically used when loading data from an XML
	 * file. This constructor extracts the necessary information from a
	 * given root DOM element. 
	 * 
	 * @param fileName The path to the cDNA data file about which this
	 * object contains statistics information. 
	 * 
	 * @param statsData The DOM element corresponding to the "Stats"
	 * object from where the data is to be extracted by this method.
	 */
	public static DataFileStats create(final String fileName, final Element statsData) {
		assert( statsData != null );
		
		int count         = DOMHelper.getIntValue(statsData, "Count");
		float avgLength   = (float) DOMHelper.getDoubleValue(statsData, "AvgLength");
		float lengthSD    = (float) DOMHelper.getDoubleValue(statsData, "LengthSD");
		int minLength     = DOMHelper.getIntValue(statsData, "MinLength");
		int maxLength     = DOMHelper.getIntValue(statsData, "MaxLength");
		long lastModified = DOMHelper.getLongValue(statsData,"LastUpdate");
		// Create a new object to be returned to the caller.
		DataFileStats dfs = 
			new DataFileStats(fileName, lastModified, count, avgLength, 
				lengthSD, minLength, maxLength);
		// Set additional sequencing technology information and
		// error rates (if available).
		String techType   = statsData.getAttribute("seqTech");
		techType          = ((techType == null) || (techType.length() == 0) ? null : techType);
		SeqTechType seqTechType = (techType != null ? SeqTechType.valueOf(techType) : 
			SeqTechType.UNKNOWN);
		dfs.setSeqInfo(seqTechType, getErrorRate(statsData, "InsErrRate"), 
				getErrorRate(statsData, "DelErrRate"), 
				getErrorRate(statsData, "SubErrRate"));
		// Return the newly created object for further use
		return dfs;
	}

	/**
	 * Helper method to get a specific error rate.
	 * 
	 * This is a a helper method that is called from the 
	 * {@link #create(String, Element)} method to obtain the 
	 * error rate value associated with the statistics object.
	 * 
	 * @param statsData The DOM element from where the error rate
	 * is to be extracted.
	 * 
	 * @param elementName The name of the child element from where
	 * the error rate is to be extracted.
	 * 
	 * @return The error rate associated with the given child element.
	 * If the child element was not found in the statsData then this
	 * method returns -1.
	 */
	private static float getErrorRate(final Element statsData, 
			final String elementName) {
		float retVal = -1;
		if (DOMHelper.hasElement(statsData, elementName)) {
			retVal = (float) DOMHelper.getDoubleValue(statsData, elementName);
		}
		return retVal;
	}
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent Workspace node in the DOM tree.
	 * 
	 * @param workspace The DOM element corresponding to the "DataSet"
	 * node that contains this entry.
	 */
	public final void marshall(Element dataSet) {
		// Create a top-level entry for this "DataSet"
		Element stats = DOMHelper.addElement(dataSet, "Stats", null);
		// Setup technology attribute value
		stats.setAttribute("seqTech", techType.name());
		// Add new sub-element for the stats node
		DOMHelper.addElement(stats, "Count",      count);
		DOMHelper.addElement(stats, "AvgLength",  avgLength);
		DOMHelper.addElement(stats, "LengthSD",   lengthSD);
		DOMHelper.addElement(stats, "MinLength",  minLength);
		DOMHelper.addElement(stats, "MaxLength",  maxLength);
		DOMHelper.addElement(stats, "LastUpdate", lastUpdate);
		// Add error rate information.
		DOMHelper.addElement(stats, "InsErrRate", insErrorRate);
		DOMHelper.addElement(stats, "DelErrRate", delErrorRate);
		DOMHelper.addElement(stats, "SubErrRate", subErrorRate);
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a XML
	 * fragment. The XML fragment is guaranteed to be compatible with the PEACE
	 * work space configuration data.
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent        = "\t\t\t";
		final String INT_ELEMENT   = Indent + "\t<%1$s>%2$d</%1$s>\n";
		final String FLOAT_ELEMENT = Indent + "\t<%1$s>%2$f</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<Stats seqTech=\"%s\">\n", Indent, techType.name()); 
		// Add new sub-elements for the ESTData element
		out.printf(INT_ELEMENT,   "Count",      count);
		out.printf(FLOAT_ELEMENT, "AvgLength",  avgLength);
		out.printf(FLOAT_ELEMENT, "LengthSD",   lengthSD);
		out.printf(INT_ELEMENT,   "MinLength",  minLength);
		out.printf(INT_ELEMENT,   "MaxLength",  maxLength);
		out.printf(INT_ELEMENT,   "LastUpdate", lastUpdate);
		// Add error rates to the XML file
		out.printf(FLOAT_ELEMENT, "InsErrRate", insErrorRate);
		out.printf(FLOAT_ELEMENT, "DelErrRate", delErrorRate);
		out.printf(FLOAT_ELEMENT, "SubErrRate", subErrorRate);
		// End the stats object
		out.printf("%s</Stats>\n", Indent); 
	}
	
	/**
	 * Obtain the path to the cDNA data file about this this object
	 * is currently holding meta data about. 
	 * 
	 * @return The path to the cDNA data file about which this object
	 * is currently holding statistics about.
	 */
	public String getPath() {
		return filePath;
	}
	
	/**
	 * Obtain the number of cDNA fragments stored in a given data file. 
	 * 
	 * This method merely returns the number of cDNA fragment entries
	 * stored in the data file associated with this entry.
	 * 
	 * @return The number of cDNA fragments in the data file associated
	 * with this object.
	 */
	public int getCount() {
		return count;
	}
	
	/**
	 * Obtain the average length of cDNA fragments stored in a 
	 * given data file. 
	 * 
	 * This method merely returns the average length of the cDNA 
	 * fragment entries stored in the data file associated with this entry.
	 * 
	 * @return The average length of cDNA fragments in the data file associated
	 * with this object.
	 */
	public float getAvgLength() {
		return avgLength;
	}

	/**
	 * Obtain the standard deviation (SD) in average length of cDNA 
	 * fragments stored in a given data file. 
	 * 
	 * This method merely returns the SD of average length of the cDNA 
	 * fragment entries stored in the data file associated with this entry.
	 * 
	 * @return The SD in average length of cDNA fragments in the data file 
	 * associated with this object.
	 */
	public float getLengthSD() {
		return lengthSD;
	}

	/**
	 * Obtain the length of the shortest cDNA fragment stored in a 
	 * given data file. 
	 * 
	 * This method merely returns the length of the shortest cDNA 
	 * fragment stored in the data file associated with this entry.
	 * 
	 * @return The length of the shortest cDNA fragment in the data file 
	 * associated with this object.
	 */
	public int getMinLength() {
		return minLength;
	}

	/**
	 * Obtain the length of the longest cDNA fragment stored in a 
	 * given data file. 
	 * 
	 * This method merely returns the length of the longest cDNA 
	 * fragment stored in the data file associated with this entry.
	 * 
	 * @return The length of the longest cDNA fragment in the data file 
	 * associated with this object.
	 */
	public int getMaxLength() {
		return maxLength;
	}
	
	/**
	 * The sequencing technology associated with this statistics object.
	 * 
	 * This method returns the technology type associated with this
	 * statistics object. 
	 *  
	 * @return The genomic sequencing technology type associated with
	 * this object (if any).
	 */
	public SeqTechType getTechType() {
		return techType;
	}

	/**
	 * The insertion error rate associated with this statistics object.
	 *  
	 * @return This method returns the insertion rate (in the range 0.0
	 * to 1.0) set for this statistics object. If a value is not know
	 * then this method returns -1. Zero indicates no insertion error.
	 */
	public float getInsErrorRate() {
		return insErrorRate;
	}

	/**
	 * The deletion error rate associated with this statistics object.
	 *  
	 * @return This method returns the deletion rate (in the range 0.0
	 * to 1.0) set for this statistics object. If a value is not know
	 * then this method returns -1. Zero indicates no deletion error.
	 */
	public float getDelErrorRate() {
		return delErrorRate;
	}

	/**
	 * The substitution error rates associated with this 
	 * statistics object.
	 *  
	 * @return This method returns the substitution error rate 
	 * (in the range 0.0 to 1.0) set for this statistics object. 
	 * If a value is not know then this method returns -1. Zero 
	 * indicates no substitution error.
	 */
	public float getSubErrorRate() {
		return subErrorRate;
	}

	/**
	 * Obtain the total error rates associated with this statistics object.
	 * 
	 * @return This method returns the total error rate which is the
	 * sum of insertion, deletion, and substitution error rates.
	 */
	public float getTotErrorRate() {
		return insErrorRate + delErrorRate + subErrorRate;
	}
	
	/**
	 * Method to write summary information using the data in this
	 * object.
	 * 
	 * This method is a convenience method that meant to be called
	 * from {@link DataSet#summarize(SummaryWriter)} method. This
	 * method adds summary information regarding the aggregate
	 * statistics about the cDNA fragments in the data file. 
	 * 
	 * @param sw The summary writer to which the data is to be written.
	 */
	public void summarize(SummaryWriter sw) {
		sw.addSummary("#cDNA fragments", "" + count, null);
		sw.addSummary("Avg. length of cDNAs", 
				String.format("%.2f (%.2f) nt", avgLength, lengthSD), 
				"cDNA lengths/size are in nucleotides (nt) or base pairs. " +
				"Value in brackets is standard deviation in cDNA length");
		sw.addSummary("Shortest cDNA", "" + minLength + " nt", null);
		sw.addSummary("Longest  cDNA", "" + maxLength + " nt", null);
	}

	/**
	 * The path to the data file about which this object currently holds
	 * the aggregate statistics/meta-data.
	 */
	private final String filePath;
	
	/**
	 * Number of cDNA fragments stored in a given data file referred to by
	 * {@link #filePath}. This value is set by the constructor and is never
	 * changed during the life time of this object.
	 */
	private final int count;

	/**
	 * Average length of cDNA fragments stored in a given data file referred to by
	 * {@link #filePath}. This value is set by the constructor and is never
	 * changed during the life time of this object.
	 */
	public final float avgLength;
	
	/**
	 * Standard Deviation in average length of cDNA fragments stored in a 
	 * given data file referred to by {@link #filePath}. This value is set
	 * by the constructor and is never changed during the life time of 
	 * this object.
	 */
	public final float lengthSD;

	/**
	 * The shortest cDNA fragments stored in a given data file referred
	 * to by {@link #filePath}. This value is set by the constructor 
	 * and is never changed during the life time of this object.
	 */
	public final int minLength;
	
	/**
	 * The longest cDNA fragments stored in a given data file referred
	 * to by {@link #filePath}. This value is set by the constructor 
	 * and is never changed during the life time of this object.
	 */	
	public final int maxLength;
	
	/**
	 * The last updated time stamp of the source cDNA file from where
	 * the aggregate statistics were calculated. This time stamp is used
	 * to determine if the statistics are up to date with the file.
	 */
	public final long lastUpdate;
	
	/**
	 * The sequencing technology value associated with this entry.
	 * This value is optional and by default it is initialized
	 * to {@link SeqTechType#UNKNOWN}.
	 */
	public SeqTechType techType = SeqTechType.UNKNOWN;
	
	/**
	 * The insertion error rate (if known) associated with this
	 * data file. If the error rate is not known then this value
	 * is set to -1.
	 */
	public float insErrorRate = -1;
	
	/**
	 * The deletion error rate (if known) associated with this
	 * data file. If the error rate is not known then this value
	 * is set to -1.
	 */
	public float delErrorRate = -1;
	
	/**
	 * The substitution type error rate (if known) associated with this
	 * data file. If the error rate is not known then this value
	 * is set to -1.
	 */
	public float subErrorRate = -1;
}
