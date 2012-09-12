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

package org.peace_tools.decagon.sdg;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PushbackInputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import org.peace_tools.data.EST;
import org.peace_tools.generic.BackgroundTask;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataFileStats;
import org.peace_tools.workspace.SeqTechType;

/**
 * Class to randomly merge multiple FASTA files into a single file. 
 *
 * This class is used in the last few steps of the Synthetic DataSet
 * Generator to merge one or more generated files into a single FASTA
 * file. This class uses multiple threads to reach each FASTA file to
 * try and make the process go faster, particularly on modern CPUs that
 * typically tend to have multiple cores. This class provides a
 * convenient interface to merge files so that i can be easily used
 * from the DataSetGenerationPage. This class can be used as shown
 * below:
 * 
 * <code>
 * FASTAFileMerger ffm = new FASTAFileMerger();
 * ffm.addSourceFiles("1.fna", -1);
 * ffm.addSourceFiles("2.fna", -1);
 * ffm.merge("TargetFile.fasta");
 * </code>
 */
public class FASTAFileMerger {
	/**
	 * The maximum number of FASTA sequences that will ever be cached from
	 * a single source file. Smaller values provide less random organizations
	 * but smaller memory footprint.
	 */
	private static final int MAX_QUEUE_SIZE = 100;
	
	/**
	 * An internal class used to encapsulate information about a source
	 * FASTA file.
	 * 
	 * This is an internal convenience class that is used to
	 * encapsulate related information about a given source file.
	 */
	private class SrcFileInfo implements Runnable, Comparator<EST> {
		/**
		 * Constructor to create a file information object.
		 * 
		 * @param srcFilePath The full path to the source FASTA file.
		 * 
		 * @param seqCountEst The estimate of the sequence count in the
		 * FASTA file.
		 * 
		 * @param techType The genomic-sequencing technology type associated
		 * with this data file. This information is used to create the
		 * aggregate statistics information object while merging data files.
		 */
		public SrcFileInfo(final String srcFilePath, 
				final long seqCountEst, final SeqTechType techType) {
			this.srcFilePath = srcFilePath;
			this.seqCountEst = seqCountEst;
			this.techType    = techType;
		}
		
		@Override
		public void run() {
			// The call gets pushed to the main class just to keep
			// this inner-class code reasonably streamlined.
			FASTAFileMerger.this.loadList(this);
		}
		
		@Override
		public int compare(EST est1, EST est2) {
			// This comparator provides some random ordering of the source
			// sequences in the estList.
			return est1.getSequence().compareTo(est2.getSequence());
		}
		
		/**
		 * The full path to the source FASTA file.
		 */
		public final String srcFilePath;
		
		/**
		 * The estimate of the sequence count in the FASTA file. If the
		 * estimate is not known then this value is set to -1.
		 */
		public final long seqCountEst;
		
		/**
		 * The sequencing technology type associated with this file. This
		 * value is used to create the DataFileStats information
		 * for this file as reads get merged.
		 */
		public final SeqTechType techType;
		
		/**
		 * The probability with which reads from this file should be chosen.
		 * This probability value is initialized based on the number
		 * of reads in this file ({@link #seqCountEst}) with respect to the 
		 * total number of reads to be merged. Files with more reads get 
		 * higher probability while files with fewer reads get lower 
		 * probability. This is done (in {@link FASTAFileMerger#initializeFileProbabilities()})
		 * to ensure a more uniform mix of reads from various
		 * files. Once the reads in a file have been exhausted, this
		 * value is set to zero and the probability is evenly 
		 * distributed to remaining files. 
		 */
		public double probability;
		
		/**
		 * List of cDNA fragments that have been read from the given FASTA
		 * file but have not yet been merged. Entries are added to this
		 * list by the {@link FASTAFileMerger#loadList(SrcFileInfo)}
		 * method. This entry intentionally uses a Priority queue in which
		 * cDNA fragments are randomly sorted based on {@link #compare(EST, EST)}
		 * method.
		 */
		public final PriorityBlockingQueue<EST> estList = 
			new PriorityBlockingQueue<EST>(MAX_QUEUE_SIZE, this);
		
		/**
		 * An exception (if any) that occurs when reading data from this
		 * FASTA file. This error is set in the {@link FASTAFileMerger#loadList(SrcFileInfo)}
		 * method. The value of this instance variable is meaningful only
		 * when the {@link #loaderThreadDone} boolean is true.
		 */
		public Exception error;
		
		/**
		 * Flag to indicate when the loader thread that reads sequences from
		 * the FASTA file has completed. This value is initialized to false.
		 * It is changed to true by the {@link FASTAFileMerger#loadList(SrcFileInfo)}
		 * method.
		 */
		public final AtomicBoolean loaderThreadDone = new AtomicBoolean(false);
		
		/**
		 * Aggregate statistics objects indicating the sequence statistics
		 * for this file. This value is set only after all the reads in
		 * the file have been successfully read. Until such time this
		 * value will be null. The value is computed and set in the
		 * {@link FASTAFileMerger#loadList(SrcFileInfo)} method.
		 */
		public DataFileStats stats;
	}

	/**
	 * The list of source FASTA files to be merged by this class. Entries 
	 * are added
	 * to this list via call to the {@link #loadSourceFiles(String, String)}
	 * or {@link #addSourceFiles(String...)} methods.
	 */
	private ArrayList<SrcFileInfo> sourceFileList = new ArrayList<SrcFileInfo>();
	
	/**
	 * The thread group that contains the background read threads created
	 * for each source FASTA file in {@link #sourceFileList}. Entries are
	 * added to this thread group by the {@link #merge(String)} method.
	 */
	private final ThreadGroup readThreadGroup = new ThreadGroup("FASTAFileMerger");
	
	/**
	 * The background task object to which logs (if any) generated by
	 * this class are to be written.
	 */
	private BackgroundTask bgTask = null;
	
	/**
	 * A simple constructor.
	 * 
	 * This is a simple constructor that more present to adhere to
	 * conventions.
	 */
	public FASTAFileMerger() {
		// Nothing much to be done here.
	}

	/**
	 * The background task to be used to generate log messages.
	 * 
	 * This method can be used to setup the background task to which
	 * log messages are to be written. This method must be invoked
	 * before the {@link #merge(String, WizardDialog)} method is called.
	 * Various method in this class use {@link BackgroundTask#log(String)}
	 * method to write log messages.
	 * 
	 * @param bgTask The background task object to which log messages
	 * are to be written by this class.
	 */
	public void setBackgroundTask(BackgroundTask bgTask) {
		this.bgTask = bgTask;
	}
	
	/**
	 * Add a source FASTA file to be merged by this class.
	 * 
	 * @parAam fileName The full path to the source FASTA file set to be 
	 * merged.
	 * 
	 * @param entryCountEstimate An estimate of the number of reads in the
	 * FASTA file. If the count is not known then this value must be -1.
	 * This value is used to provide a more uniform random distribution of
	 * the data from various files without incurring too much overhead to
	 * estimate file sizes.
	 * 
	 * @exception IllegalArgumentException This method throws this exception
	 * if the supplied file could not be read or is an invalid FASTA file.
	 */
	public void addSourceFile(final String fileName, 
			final long entryCountEstimate, final SeqTechType techType) throws IllegalArgumentException {
		// Validate the file to check if it has a ">" character first.
		int firstChar = 0; // First character in the given file.
		try {
			FileInputStream fis = new FileInputStream(fileName);
			firstChar = fis.read();
			fis.close();
		} catch (IOException ioe) {
			ProgrammerLog.log(ioe);
			throw new IllegalArgumentException("Unable to read FASTA file '" + fileName + "'");
		}
		if (firstChar != '>') {
			// The user says it should be an sff file but it is not.
			throw new IllegalArgumentException("The given file '" + fileName + 
					"' is not a valid FASTA file as it does not begin " + 
					"with a '>' character.");
		}
		sourceFileList.add(new SrcFileInfo(fileName, entryCountEstimate, techType));
	}

	/**
	 * Helper method to cut logs if background task is set.
	 * 
	 * This is an internal helper method that was introduced to streamline
	 * the process of generating progress logs. This method
	 * cuts logs only if the {@link #bgTask} member is not null.
	 * 
	 * @param logEntry The progress log information to be generated. This
	 * parameter cannot be null. A new-line character is appended to this
	 * string prior to logging.
	 */
	private void log(final String logEntry) {
		if (bgTask != null) {
			bgTask.log(logEntry + "\n");
		}
	}
	
	/**
	 * Method to initialize probabilities for each source FASTA file.
	 * 
	 * This method is called only once from the {@link #merge(String)}
	 * method to setup the probabilities with which cDNA entries are
	 * read from a given FASTA source file. This probability value is 
	 * initialized based on the estimated number of reads in each file 
	 * ({@link SrcFileInfo#seqCountEst}) with respect to the total number 
	 * of reads to be merged. Files with more reads get 
	 * higher probability while files with fewer reads get lower 
	 * probability. This is done by this method 
	 * to ensure a more uniform mix of reads from various
	 * files. Once the reads in a file have been exhausted, this
	 * value is set to zero and the probability is evenly 
	 * distributed to remaining files. 
	 */
	private void initializeFileProbabilities() {
		// Compute total number of reads to be processed to determine
		// relative probability for each file. We have to handle case
		// were the estimated number of reads is not really known as well.
		long totalSeqCount   = 0; // Total count of sequences
		long numFilesCounted = 0; // The number of files with seqCount > 0
		for(SrcFileInfo sfi: this.sourceFileList) {
			if (sfi.seqCountEst > 0) {
				totalSeqCount += sfi.seqCountEst;
				numFilesCounted++;
			}
		}
		// Handle case where none of the files had a useful
		// sequence count estimate value. In this case we assign
		// the same probabilities to all the files.
		if (numFilesCounted == 0) {
			totalSeqCount = numFilesCounted = sourceFileList.size();
		} else if (numFilesCounted < sourceFileList.size()){
			// At least one or more of the files had useful sequence
			// count values. For those files that did not have a useful
			// sequence count, assume the value to the average number of
			// reads.
			final long noCountEntries = sourceFileList.size() - numFilesCounted;
			// Add average sequence count for files without seqCount set
			totalSeqCount  += ((totalSeqCount / numFilesCounted) * noCountEntries);
			numFilesCounted = sourceFileList.size();
		}
		// Determine average number of reads in each file (which can be different
		// after the flow above).
		final long avgSeqCount = totalSeqCount / numFilesCounted;
		// Now assign initial probabilities to each file.
		for(SrcFileInfo sfi: this.sourceFileList) {
			final double seqCount = (sfi.seqCountEst > 0 ? sfi.seqCountEst : avgSeqCount);
			sfi.probability = seqCount / totalSeqCount;
		}
	}
	
	/**
	 * Get a random read from the given source FASTA file and handle 
	 * end-of-file.
	 * 
	 * <p><b>NOTE</b>:This method must be invoked from the same thread
	 * on which the {@link #merge(String)} method was invoked.</p>
	 * 
	 * <p>This is a helper method that returns the next read from the given
	 * source FASTA file (if a read is available). This method extracts
	 * the first read from {@link SrcFileInfo#estList}. The entries in
	 * {@link SrcFileInfo#estList} would be in random order (as it is a 
	 * priority queue and priority is randomly determined). If the 
	 * {@link SrcFileInfo#estList} is empty then it polls the queue until
	 * a valid read is obtained or until the loader thread signals it is
	 * done (by setting {@link SrcFileInfo#loaderThreadDone} flag).<p>
	 * 
	 * If all reads from the source file have been read then this method
	 * then this method performs the following tasks:
	 * 
	 * <ol> 
	 *   <li>If any exceptions were reported in {@link SrcFileInfo#error} then
	 *   this method throws this exception on this (main) thread.</li>
	 * 
	 *   <li>Next (when no exceptions were reported) this method redistributes
	 *   the probability set for this file to other source files that still
	 *   have some reads to be processed. This ensures reads from remaining
	 *   files are uniformly distributed.</li>
	 *   
	 *   <li>Finally, the probability for this file (the one in which EOF has
	 *   been reached) is set to zero so that this entry will be ignored 
	 *   in future processing in the {@link #getRandomRead()} method.</li>
	 * </ol>
	 * 
	 * @param sfi The source FASTA file object from where the next
	 * cDNA entry is to be extracted and returned. If all fragments from
	 * this source have been processed, then this entry's probability is
	 * set to zero to prevent it from being used in further processing.
	 * 
	 * @return This method returns null when all fragments from this
	 * source FASTA file have been successfully processed. Otherwise it
	 * returns a valid EST object.
	 * 
	 * @throws Exception This method exposes any exceptions that may 
	 * have occurred on the loader thread or on this main thread.
	 */
	private EST getRandomRead(final SrcFileInfo sfi) throws Exception {
		// Repeatedly try to get a read from the given source file
		// until we get a valid read or the background thread finishes.
		EST srcRead = null;
		do {
			// Try to get a read in 50 milliseconds (the background
			// thread may still be reading or has reached EOF).
			srcRead = sfi.estList.poll(50, TimeUnit.MILLISECONDS);	
		} while ((srcRead == null) && !sfi.loaderThreadDone.get());
		// If loader thread is done then we have to redistribute its
		// probability to other files giving them preference and clearing
		// out this files probability.
		if ((srcRead == null) && (sfi.loaderThreadDone.get())) {
			// Check to ensure no error occurred. If error occurred throw
			// it on this thread.
			if (sfi.error != null) {
				throw sfi.error;
			}
			// Report progress.
			log("    All reads merged from " + sfi.srcFilePath);
			// The loader thread is done. Reassign probabilities evenly.
			// For this we first need to determine number of other FASTA
			// files that are currently pending processing.
			int numPendingFiles = 0;
			for(SrcFileInfo otherSFI: sourceFileList) {
				if ((otherSFI != sfi) && (otherSFI.probability > 0)) {
					numPendingFiles++;
				}
			}
			// Now distribute probabilities evenly if we have any left.
			if (numPendingFiles > 0) {
				final double extraProb = sfi.probability / numPendingFiles;
				for(SrcFileInfo otherSFI: sourceFileList) {
					if ((otherSFI != sfi) && (otherSFI.probability > 0)) {
						otherSFI.probability += extraProb;
					}
				}
			}
			// Clear out the probability for the current entry
			sfi.probability = 0;
		}
		return srcRead;
	}
	
	/**
	 * Randomly obtain a read from one of the source FASTA files.
	 * 
	 * <p><b>NOTE</b>:This method must be invoked from the same thread
	 * on which the {@link #merge(String)} method was invoked.</p>
	 * 
	 * This is a top-level method that is repeatedly invoked from
	 * the {@link #merge(String)} method to obtain the next cDNA
	 * fragment to be written to the merged FASTA file. This method 
	 * generates a uniform random number and iterates over the entries 
	 * in the {@link #sourceFileList} to locate the source file that
	 * maps to this random number using {@link SrcFileInfo#probability}
	 * value. Once a valid file entry is located it uses the 
	 * {@link #getRandomRead(SrcFileInfo)} to try and obtain a cDNA read
	 * from it. If a valid cDNA read is obtained then this method returns
	 * that value. Otherwise it retries until a valid entry is obtained or
	 * if all entries have been processed then this method returns null. 
	 * 
	 * @return This method returns null if all cDNA reads from all the
	 * source FASTA files have already been read and processed. Otherwise
	 * this method returns a valid (not null) EST object be written to
	 * the target FASTA file.
	 * 
	 * @throws Exception This method exposes any exceptions that may 
	 * have occurred on the loader thread or on this main thread.
	 */
	private EST getRandomRead() throws Exception {
		// Generate a random number to determine the file from where the
		// next cDNA read is to be obtained.
		final double fileSelect = Math.random();
		// Try and locate the next random read. We need a loop here to
		// handle case when we reach end of FASTA file and we need to 
		// try again.
		int currSFI           = 0;
		EST srcRead           = null;
		double currFileSelect = fileSelect; // updated in while-loop below
		// Locate the file that contains the given probability.
		while (currSFI < sourceFileList.size()) {
			final SrcFileInfo sfi = sourceFileList.get(currSFI);
			if (currFileSelect < sfi.probability) {
				// Found a potential file entry. Try to get a
				// random read from the file.
				if ((srcRead = getRandomRead(sfi)) != null) {
					// Found a valid entry.
					break;
				} else {
					// We exhausted reads in this file. So
					// have to start over.
					currSFI        = 0;
					currFileSelect = fileSelect;
				}
			} else {
				// Onto the next valid file entry.
				currFileSelect -= sfi.probability;
				currSFI++;
			}
		}
		// Either we have a non-null entry or we have read
		// all reads from all source FATA files (in which case
		// this method will return null as expected).
		return srcRead;
	}
	
	/**
	 * Merge the set of source FASTA files into a single FASTA file.
	 * 
	 * This is the main method that performs the multi-file merge
	 * operation. The operation can take a while depending on the
	 * file sizes. To improve overall performance this method uses
	 * a multi-threaded approach for merging files. It creates 
	 * multiple background threads to read data from the source
	 * FASTA files. The background threads populate the 
	 * {@link SrcFileInfo#estList}. This method uses the
	 * {@link #getRandomRead()} method to obtain entries from the
	 * various {@link #sourceFileList} and write them to the target
	 * file until all entries have been merged into the target file.
	 * 
	 * @param targetFilePath The target path where the FASTA file is
	 * to be stored. This file should be writeable. Otherwise an
	 * exception will be thrown
	 * 
	 * @param wd An optional reference to a WizardDialog that is invoking
	 * this method. If this object is not null, then {@link WizardDialog#addThread(Thread)}
	 * and {@link WizardDialog#removeThread(Thread)} methods are used to
	 * register and unregister background reader threads to ensure that
	 * threads are cleaned upon user canceling the wizard.
	 * 
	 * @exception Exception This method exposes various exceptions that
	 * may occur when merging FASTA files.
	 * 
	 * @return This method returns an DataFileStats object that contains
	 * aggregate statistics about the target file.
	 */
	public ArrayList<DataFileStats> merge(String targetFilePath, WizardDialog wd) throws Exception {
		// Ensure we have at least one file to merge
		if (sourceFileList.size() < 1) {
			throw new IllegalArgumentException("No source FASTA files were specified for merging");
		}
		// Setup initial probabilities.
		initializeFileProbabilities();
		// Open the target file to be written.
		final PrintStream targetFile = 
			new PrintStream(new FileOutputStream(targetFilePath));
		// Kick-off background read threads for each of the source files.
		// Track threads for use later on.
		ArrayList<Thread> readerThreadList = new ArrayList<Thread>(sourceFileList.size());
		for(SrcFileInfo sfi: sourceFileList) {
			Thread bgReader = new Thread(readThreadGroup, sfi, sfi.srcFilePath);
			if (wd != null) {
				wd.addThread(bgReader);
			}
			readerThreadList.add(bgReader);
			bgReader.start();
		}
		// Repeatedly read a random read from a randomly chosen source
		// FASTA file and write it to the destination file. 
		EST srcEST    = null;
		int numReads  = 0;
		// Process until all reads are merged...
		while ((srcEST = getRandomRead()) != null) {
			srcEST.write(targetFile);
			// Track statistics
			numReads++;
			// Display intermediate results.
			if (numReads % 1000 == 0) {
				log("    Merged " + numReads + " reads...");
			}
		}
		targetFile.close();
		log("    Merged complete. Total reads merged = " + numReads);
		// Remove background threads we registered.
		for(Thread t: readerThreadList) {
			if (wd != null) {
				wd.removeThread(t);
			}
		}
		// Gather the statistics object from all worker threads.
		ArrayList<DataFileStats> statsList = new ArrayList<DataFileStats>();
		for(SrcFileInfo sfi: sourceFileList) {
			statsList.add(sfi.stats);
		}
		// Return aggregate statistics of the target file back to the caller.
		return statsList;
	}
	
	/**
	 * Helper method to load cDNA reads from a FASTA file.
	 * 
	 * This method is called from {@link SrcFileInfo#run()} method 
	 * which is called from different threads (one for each source
	 * FASTA file to be merged). This method is relatively 
	 * straightforward. It continuously loads reads from the
	 * given FASTA file and places them in {@link SrcFileInfo#estList}.
	 * However, this method will block if {@link SrcFileInfo#estList}
	 * size reaches {@link FASTAFileMerger#MAX_QUEUE_SIZE}. If
	 * errors are encountered then this method sets the exception
	 * in {@link SrcFileInfo#error}. This method always sets
	 * {@link SrcFileInfo#loaderThreadDone} to true just before
	 * it returns.
	 * 
	 * <p><b>NOTE</b>: Unlinke other methods in this class, this
	 * method is called from different threads. So use caution when
	 * inspecting this method and its operations.</p>
	 * 
	 * @param sfi The class that contains all the information about
	 * the source FASTA file on which this method is to operate.
	 */
	private void loadList(final SrcFileInfo sfi) {
		try {
			// Open the source file for reading.
			InputStream is = new BufferedInputStream(new FileInputStream(sfi.srcFilePath));
			// Wrap input stream into push back stream to load FASTA data
			final PushbackInputStream fasta = new PushbackInputStream(is);
			log("    Loading reads from file " + sfi.srcFilePath);
			// Track statistics as reads are loaded.
			int numReads = 0;
			int  minLen = Integer.MAX_VALUE, maxLen = 0;
			long lenSum = 0, lenSqSum = 0;
			// Load ESTs...
			while (fasta.available() > 0) {
				// Create dummy cDNA.
				EST est = new EST(sfi.estList.size(), null, null);
				// Get EST to load itself.
				est.read(fasta);
				// Add cDNA to the list of ESTs. The following put()
				// call will not block. However to keep memory in 
				// check we introduce a loop here to ensure that
				// the size never becomes too large.
				while (sfi.estList.size() > MAX_QUEUE_SIZE) {
					// The queue is large. Wait for some entries to
					// be removed from the queue.
					Thread.sleep(50);
				}
				// Add entry to the queue of things to be merged
				sfi.estList.put(est);
				numReads++;
				// Track statistics
				int seqLen = est.getSequence().length();
				minLen     = Math.min(minLen, seqLen);
				maxLen     = Math.max(maxLen, seqLen);
				lenSum    += seqLen;
				lenSqSum  += (seqLen * seqLen);
			}
			// Close the stream.
			is.close();
			// Print number of reads loaded from the given file.
			log("    Number of reads from file " + sfi.srcFilePath + 
					" = " + numReads);
			// Generate aggregate statistics object.
			float avgLen = (float) (lenSum / (float) numReads);
			float lenSD  = (float) (Math.sqrt((lenSqSum - (lenSum * avgLen)) / numReads));
			sfi.stats = new DataFileStats(sfi.srcFilePath, -1, numReads,
					avgLen, lenSD, minLen, maxLen);
			sfi.stats.setSeqInfo(sfi.techType, -1, -1, -1);
		} catch (Exception exp) {
			ProgrammerLog.log(exp);
			sfi.error = exp;
		} finally {
			// Set flag to indicate that this thread is done
			sfi.loaderThreadDone.set(true);
		}
	}
}	

// End of source code