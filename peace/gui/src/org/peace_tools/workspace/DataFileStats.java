package org.peace_tools.workspace;

import java.io.File;
import java.io.PrintWriter;

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
	
	/**
	 * Convenience constructor to load data from a DOM Tree
	 * 
	 * This constructor is typically used when loading data from an XML
	 * file. This constructor extracts the necessary information from a
	 * given root DOM element. This element is typically generated
	 * via a call to the 
	 * @param fileName
	 * @param statsData
	 */
	public static DataFileStats create(final String fileName, final Element statsData) {
		assert( statsData != null );
		
		int count         = DOMHelper.getIntValue(statsData, "Count");
		float avgLength   = (float) DOMHelper.getDoubleValue(statsData, "AvgLength");
		float lengthSD    = (float) DOMHelper.getDoubleValue(statsData, "LengthSD");
		int minLength     = DOMHelper.getIntValue(statsData, "MinLength");
		int maxLength     = DOMHelper.getIntValue(statsData, "MaxLength");
		long lastModified = DOMHelper.getLongValue(statsData,"LastUpdate");
		// Create a new object and return it back to the caller
		return new DataFileStats(fileName, lastModified, count, avgLength, 
				lengthSD, minLength, maxLength);
	}

	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent Workspace node in the DOM tree.
	 * 
	 * @param workspace The DOM element corresponding to the "Workspace"
	 * node that contains this entry.
	 */
	public final void marshall(Element dataSet) {
		// Create a top-level entry for this "DataSet"
		Element stats = DOMHelper.addElement(dataSet, "Stats", null);
		// Add new sub-element for the stats node
		DOMHelper.addElement(stats, "Count",      count);
		DOMHelper.addElement(stats, "AvgLength",  avgLength);
		DOMHelper.addElement(stats, "LengthSD",   lengthSD);
		DOMHelper.addElement(stats, "MinLength",  minLength);
		DOMHelper.addElement(stats, "MaxLength",  maxLength);
		DOMHelper.addElement(stats, "LastUpdate", lastUpdate);
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a XML
	 * fragment. The XML fragment is guaranteed to be compatible with the PEACE
	 * work space configuration data.
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent        = "\t\t";
		final String INT_ELEMENT   = Indent + "<%1$s>%2$d</%1$s>\n";
		final String FLOAT_ELEMENT = Indent + "<%1$s>%2$f</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<Stats>\n", Indent); 
		// Add new sub-elements for the ESTData element
		out.printf(INT_ELEMENT, "Count",        count);
		out.printf(FLOAT_ELEMENT, "AvgLength",  avgLength);
		out.printf(FLOAT_ELEMENT, "LengthSD",   lengthSD);
		out.printf(INT_ELEMENT,   "MinLength",  minLength);
		out.printf(INT_ELEMENT,   "MaxLength",  maxLength);
		out.printf(INT_ELEMENT,   "LastUpdate", lastUpdate);
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
}
