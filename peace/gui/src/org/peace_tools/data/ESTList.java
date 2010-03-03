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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PushbackInputStream;
import java.util.ArrayList;

import javax.swing.ProgressMonitor;

import org.peace_tools.workspace.ClassifierList;
import org.peace_tools.workspace.DBClassifier;
import org.peace_tools.workspace.Workspace;

/**
 * Class to encapsulate a list of ESTs.
 * 
 * This class has been designed to contain a list of ESTs that are typically
 * loaded from a FASTA file. The ESTs are maintained as an ordered list of
 * values so that the index of the ESTs are used as IDs to look up ESTs.
 * 
 */
public class ESTList {	
	/**
	 * Create an empty list of ESTs tagged with a given name.
	 * 
	 * @param name An identifier to be associated with this ESTList. 
	 * Typically the absolute path to the FASTA file from where the data 
	 * was loaded is used as the identifier. This information helps to 
	 * locate ESTs that were loaded from a given file. 
	 */
	public ESTList(String name) {
		this.name = name;
		this.ests = new ArrayList<EST>();
	}
	
	/**
	 * The name associated to identify this ESTList.
	 * 
	 * @return The name set for this ESTList when it was created. The
	 * name can be null if a name was not set.
	 */
	public String getName() { return name; }
	
	/**
	 * Obtain the list of ESTs encapsulated by this ESTList.
	 * 
	 * @return The list of ESTs encapsulated by this EST list. The list
	 * may be empty.
	 */
	public ArrayList<EST> getESTs() { return ests; }
	
	/**
	 * Load ESTs from a given FASTA file and build a new ESTList.
	 * 
	 * This method can be used to load ESTs from a given FASTA file into
	 * a new ESTList. The absolute path to the FASTA file is used as the
	 * name for the newly created ESTList.
	 *  
	 * @param fastaFile The FASTA file from where the ESTs are to be loaded.
	 * 
	 * @return The newly cerated EST list containing all the ESTs in the
	 * FASTA file.
	 * 
	 * @exception IOException This method throws an IO exception on errors.
	 */
	public static ESTList loadESTs(File fastaFile) throws IOException {
		// Create a new input stream to read ESTs
		InputStream is = new FileInputStream(fastaFile);
		is = new BufferedInputStream(is);
		// Use the helper method to load the data.
		return loadESTs(fastaFile.getAbsolutePath(), is);
	}
	
	/**
	 * Load ESTs from a stream that provides FASTA data.
	 * 
	 * This method can be used to load ESTs from a given input stream into
	 * a new ESTList. The absolute path to the FASTA file to be used to 
	 * name/tag the newly created ESTList must be provided.
	 *  
	 * @param fileName The absolute path to the FASTA file to be used to
	 * name/tag the newly created ESTList.
	 * 
	 * @param is The input stream from where the EST data is to be read.
	 * 
	 * @return The newly created EST list containing all the ESTs in the
	 * FASTA file.
	 * 
	 * @exception IOException This method throws an IO exception on errors.
	 */
	public static ESTList loadESTs(String fileName, InputStream is) throws IOException {
		// First create an empty list.
		ESTList list = new ESTList(fileName);
		// Get the ests from the list for handy reference
		ArrayList<EST> ests = list.getESTs();
		// Wrap input stream into push back stream to load FASTA data
		PushbackInputStream fasta = new PushbackInputStream(is);
		// Load ESTS...
		while (fasta.available() > 0) {
			// Create dummy est.
			EST est = new EST(ests.size(), null, null);
			// Get EST to load itself.
			est.read(fasta);
			// Add est to the list of ESTs
			ests.add(est);
		}
		// Close the stream.
		is.close();
		// Let the caller deal with this list of ESTs.
		return list;
	}
	
	/**
	 * Writes ESTs in FASTA format to a given writer.
	 * 
	 * This method can be used to write the ESTs encapsulated in this list
	 * to a given stream in FASTA format.
	 *  
	 * @param writer The output stream to which the ESTs are to be written.
	 * 
	 * @exception IOException This method throws an IO exception on errors.
	 */
	public void writeESTs(PrintStream writer) throws IOException {
		for(int index = 0; (index < ests.size()); index++) {
			ests.get(index).write(writer);
		}
	}
	
	/**
	 * Writes ESTs in FASTA format to a given file.
	 * 
	 * This method can be used to write the ESTs encapsulated in this list
	 * to a given stream in FASTA format.
	 *  
	 * @param outFile The file to which the ESTs are to be written.
	 * 
	 * @exception IOException This method throws an IO exception on errors.
	 */
	public void writeESTs(File outFile) throws IOException {
		// Create an print writer that wraps an output file stream.
		PrintStream writer = new PrintStream(outFile);
		for(int index = 0; (index < ests.size()); index++) {
			ests.get(index).write(writer);
		}
	}

	/**
	 * Method to compute and return statistics about the set of 
	 * ESTs encapsulated in this list. 
	 * 
	 * @return This method returns an array of doubles containing 
	 * the following information: estCount, minLen, maxLen, avgLen, lenSD.
	 */
	public double[] computeStatistics() {
		// Compute the basic sequence length statistics.
		int  minLen = Integer.MAX_VALUE, maxLen  = 0;
		long lenSum = 0, lenSqSum= 0;
		
		for(int idx = 0; (idx < ests.size()); idx++) {
			String seq = ests.get(idx).getSequence();
			int seqLen = seq.length();
			minLen    = Math.min(minLen, seqLen);
			maxLen    = Math.max(maxLen, seqLen);
			lenSum   += seqLen;
			lenSqSum += (seqLen * seqLen);
		}
		// Compute the statistics and store them in an array to return back
		double stats[] = new double[5];
		stats[0] = ests.size();
		stats[1] = minLen;
		stats[2] = maxLen;
		if (ests.size() > 0) {
			stats[3] = (lenSum / (double) ests.size());
			stats[4] = Math.sqrt((lenSqSum - (lenSum * stats[3])) / 
					(ests.size() -1 ));
		}
		return stats;
	}
	
	/**
	 * Determine if the EST list is current with work space classifiers.
	 * 
	 * This method can be used to determine if this EST list has already
	 * been classified using the current set of classifiers configured for
	 * this work space. 
	 * 
	 * @return This method returns true if the EST list has already been
	 * classified using the current set of classifiers.
	 */
	public synchronized boolean isClassified() {
		ClassifierList clasList = Workspace.get().getClassifierList();
		ArrayList<DBClassifier> classifiers = clasList.getClassifiers();
		// Intentional use of == to check references!
		return (classifiers == prevClusterList);
	}
	
	/**
	 * Method classify all the ESTs in this list.
	 * 
	 * <p>This method is typically invoked directly from a view (rather than from
	 * the model associated with the view) to classify all the ESTs in the list
	 * based on the current classifications set in the work space. This method
	 * iterates over all the ESTs and has the ESTs compute their classifications.</p>
	 * 
	 * <p>This method recomputes classifications only if the classifier list in the 
	 * work space has actually changed. Consequently, repeatedly calling this
	 * method does not have side effects. However, if classification is performed
	 * then it may be a long running process particularly, for large EST sets.
	 * Consequently, it is best to invoke this method from a separate thread
	 * so that the GUI does not appear to be hanging while classifications are
	 * computed.</p>
	 * 
	 * @param pm An optional progress monitor to be updated as ESTs are
	 * classified. The progress monitor provides feedback to the user about
	 * the progress. This parameter can be null.
	 */
	public synchronized void classify(ProgressMonitor pm) {
		ClassifierList clasList = Workspace.get().getClassifierList();
		ArrayList<DBClassifier> classifiers = clasList.getClassifiers();
		if (prevClusterList != classifiers) {
			for(int i = 0; (i < ests.size()); i++) {
				ests.get(i).classify(classifiers);
				// Update progress if monitor is not null.
				if (pm != null) {
					pm.setProgress(i); 
				}
			}
			// Track the latest classification information used.
			prevClusterList = classifiers;
		}
	}
	
	/**
	 * This instance variable contains the list of ESTs encapsulated in this
	 * ESTList class. The set of ESTs can be accessed via a call to the
	 * getESTs() method.
	 */
	private ArrayList<EST> ests;
	
	/**
	 * A unique name or identifier associated with this EST list. Typically
	 * the absolute path to the FASTA file from where the data was loaded
	 * is used as the identifier. This information helps to locate ESTs that
	 * were loaded from a given file. 
	 */
	private final String name;
	
	/**
	 * This reference is used to track the object used to hold the list
	 * of DBClassifiers in the work space. This information is used to
	 * decide if the classifications need to be recomputed. This works
	 * because the classifiers are modified as a single complete batch
	 * of entries. Each time a new object is set and if the objects have
	 * changed then the classifications need to be recomputed for this
	 * cluster file.
	 */
	private Object prevClusterList;
}
