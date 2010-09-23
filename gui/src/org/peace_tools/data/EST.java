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

package org.peace_tools.data;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PushbackInputStream;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.peace_tools.workspace.DBClassifier;

public class EST {
    /** EST Constructor.
     * 
     * This constructor is used to instantiate an EST method.
     * 
     * @param id The unqiue ID value to be set for this EST.
     * 
     * @param info The name and other information associated with the EST.
     * This information is typically the first header line read from a FASTA
     * file.  This information can be null.
     * 
     * @param sequence The actual sequence of base pairs associated with this
     * EST.  The sequence information can be null.
     */
	EST(int id, String info, String sequence) {
		this.id           = id;
		this.info         = info;
		this.sequence     = sequence;
		this.dbClassifier = -1;
	}

	/** EST copy Constructor.
     * 
     * This constructor is used to instantiate an EST object.
     * 
     * @param src The source EST from where the data is to be copied. The unqiue ID value to be set for this EST.
     */
	EST(final EST src) {
		this.id           = src.id;
		this.info         = src.info;
		this.sequence     = src.sequence;
		this.dbClassifier = src.dbClassifier;
	}
	
	/** Obtain the ID of this EST.
	 * 
	 * @return The ID of the EST that was set when this EST was created.
	 */
	public int getID() { return id; }

    /** Obtain the information associated with this EST.
     * 
     *  The name and other information associated with the EST.  This
     *  information is typically the first header line read from a FASTA
     *  file.
     *  
     *  @return Any information available for this EST.
     */
	public String getInfo() { return info; }

    /** Obtain the actual sequence of base pairs for this EST.
     * 
     *  Note that sequence inforamtion for an EST can be null if it does
     *  not have a valid sequence information associated with it.
     *  
     *  @return The actual sequence of base paris for this EST.
     */
	public String getSequence() { return sequence; }

	/**
	 * Obtain the classifier index with which this EST is associated.
	 * 
	 * The classifier which which this EST is associated. The return
	 * value from this method is meanigful only after the 
	 * ESTList.classify() method has been invoked to classify ESTs 
	 * using the current set of classifiers.
	 * 
	 * @return The classifier list with which this EST is associated.
	 * If this EST is not associated with a classifier then this 
	 * method returns -1.
	 */
	public int getDBClassifier() { return dbClassifier; }
	
	/**
	 * Method to dump the EST data out in FASTA file format. 
	 * 
	 * This method can be used to dump out the EST information in a 
	 * FASTA compatible format. 
	 * 
	 * @param os The output stream to which the EST must be written in
	 * a FASTA format.
	 */
	public void write(PrintStream os) {
		os.println(">" + info);
		// Write sequences out with no more than 80 chars per line.
		String seqRemaining = sequence;
		while (seqRemaining.length() > 0) {
			int bp2Print  = Math.min(80, seqRemaining.length());
			String subSeq = seqRemaining.substring(0, bp2Print);
			seqRemaining  = seqRemaining.substring(bp2Print);
			os.println(subSeq);
		}
	}
	
	/**
	 * Method to read an EST from a given input stream.
	 * 
	 * This method can be used to read an EST from a given input stream. This
	 * method requires an push-back stream so that at least one character
	 * can be undone. This one character push-back is needed to detect the
	 * end of an EST sequence and the beginning of the other by checking to
	 * see if the next character to be read is a ">" sign.
	 * 
	 * <p><b>Note:</b>  This method updates only the info and sequence 
	 * information associated with this EST and not its id value.</p>
	 * 
	 * @param is The input stream from where the EST is to be read.
	 */
	public void read(PushbackInputStream is) throws IOException {
		// First check to ensure that the first character is indeed a ">"
		int greater = is.read();
		if (greater != '>') {
			throw new IOException("EST entry in FASTA file did not start with a '>'");
		}
		// Read rest of the line as the info
		info = readline(is);
		if (info == null) {
			throw new IOException("EST entry in FASTA file did not have identifier after '>'");
		}
		info = info.trim();
		// Keep reading fasta sequences until we hit EOF or the next sequence
		StringBuffer bases = new StringBuffer(1000);
		do {
			// Peek at the next character to see if it is a '>'
			greater = is.read();
			if (greater != -1) {
				is.unread(greater);
			}
			// if the character was a ">" or -1 (eof) then bail out.
			if ((greater == '>') || (greater == -1)) {
				break;
			}
			// Accumulate base pairs into a buffer
			bases.append(readline(is).trim());
		} while (is.available() > 0);
		// Drop any blanks off
		sequence = bases.toString().trim();
	}
	
	/**
	 * Helper method to read a line from a given input stream. This method
	 * repeatedly reads characters from the input stream until a newline
	 * character is read. 
	 *  
	 * @param is The input stream from where a line is to be read.
	 * 
	 * @return A string containing the line read. If the no data was available
	 * then this method returns null. Otherwise it reaturns the next line read
	 * without the trailing new line character. 
	 * 
	 * @throws IOException This method throws an exception one errors.
	 */
	protected String readline(InputStream is) throws IOException {
		if (is.available() < 1) {
			return null;
		}
		// Use a string buffer that can hold at least 256 characters
		StringBuffer retVal = new StringBuffer(256);
		int nextChar = 0;
		while ((nextChar = is.read()) != -1) {
			if (nextChar == '\n') {
				// End of line detected.
				break;
			}
			retVal.append((char) nextChar);
		}
		// return the buffer back to the caller.
		return retVal.toString();
	}
	
	@Override
	public String toString() {
		return "" + id + ": " + info;
	}
	
	/**
	 * Utility methods to facilitate searching for information in this EST.
	 * This method searches for the given searchStr within all the information
	 * including: EST Index (ID), the FASTA header, and nucleotide sequence.
	 * 
	 * @param searchStr The string to search for.
	 * @return This method returns true if the EST matches or contains the
	 * given searchStr.
	 * 
	 * @param caseSensitive If this flag is false, then this method converts string
	 * to lower case to make a case insensitive search.
	 */
	public boolean contains(String searchStr, boolean caseSensitive) {
		String searchInfo[] = {"" + id, info, sequence};
		for(String str: searchInfo) {
			if (!caseSensitive) {
				str = str.toLowerCase();
			}
			if (str.indexOf(searchStr) != -1) {
				return true;
			}
		}
		// Search string not found.
		return false;
	}
	
	/**
	 * Utility methods to facilitate searching for information in this EST.
	 * This method searches for the given searchStr within all the information
	 * including: EST Index (ID), the FASTA header, and nucleotide sequence.
	 * 
	 * @param regExp The regular expression that is being used for the search.
	 * 
	 * @return This method returns true if the EST's information matches the
	 * given regular expression's grammar.
	 */
	public boolean contains(Pattern regExp) {
		final String info = toString() + ": " + sequence;
		Matcher m = regExp.matcher(info);
		return m.matches();
	}
	
	/**
	 * Method to determine the classifier that identifies this EST.
	 * 
	 * This method is typically invoked from the ESTList.classify() method
	 * to determine the first, enabled classifier that identifies or matches
	 * the information (from the FASTA header) associated with this EST.
	 * This method iterates over the list of classifiers and uses the 
	 * regular expression in each classifier to identify the matching
	 * classification for this EST. 
	 * 
	 * @param classifiers The list of classifiers to be used to determine
	 * the matching classifier. Typically, this list is the same as the 
	 * list of classifiers currently configured in the work space.
	 */
	protected void classify(final ArrayList<DBClassifier> classifiers) {
		dbClassifier = -1;
		for(int i = 0; (i < classifiers.size()); i++) {
			if (!classifiers.get(i).isEnabled()) {
				// This classifier is not enabled and should not be used.
				continue;
			}
			Pattern pattern = classifiers.get(i).getPattern();
			Matcher matcher = pattern.matcher(info);
			if (matcher.matches()) {
				// Found a matching classifier.
				dbClassifier = i;
				break;
			}
		}
	}
	
	/**
	 * The unique ID for this EST. This member holds the unique ID for this EST.
	 * The ID is set when the EST is instantiated and is never changed during
	 * the life time of this EST. The id is used to access and extract EST
	 * information.
	 */
	private int id;

	/**
	 * The name and other information associated with the EST. This information
	 * is typically the first header line read from a FASTA file.
	 */
	private String info;

    /** The actual sequence of base pairs associated with this EST.
     * 
     * This information is typically read from a FASTA file.  
     * Currently it is assumed that an EST sequence can have only the
     * four standard base pairs (A, T, G, and G).
     */
	private String sequence;

	/**
	 * Index of the classifier in the Work space's classifier list with
	 * which this EST is associated. This information changes whenever a
	 * new classifier list is set. If this EST does not belong to a
	 * specific classifier then this value is set to -1.
	 */
	private int dbClassifier;
}
