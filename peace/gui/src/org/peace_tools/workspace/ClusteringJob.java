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

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.FWAnalyzer.FWAnalyzerType;
import org.peace_tools.workspace.Filter.FilterType;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * <p>This class corresponds to a "ClusteringJob" element in a 
 * PEACE work space XML data. This class encapsulates the core 
 * information associated with a Job. This class serves merely
 * as a read-only type class that is created when a new job is run.
 * The data is persisted in the XML work space for future reference
 * so that users can determine jobs that were scheduled and check on
 * status of long-running jobs.</p>
 * 
 * <p>In addition to encapsulating the data this class also provides
 * convenient interfaces for marshaling and un-marshaling XML data compatible
 * with the PEACE GUI configuration XML file format. </p>
 */
public class ClusteringJob extends Job {
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * Job entry. This method is typically  used to create a suitable
	 * Job entry when loading a Work space into the GUI.
	 * 
	 * @param jobNode The DOM element to be used for creating the job
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created job entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static ClusteringJob create(Element jobNode) throws Exception {
		// First extract the necessary information to create a minimally
		// populated clustering job object from the DOM tree.
		String jobID    = DOMHelper.getStringValue(jobNode, "JobID");
		String srvrID   = DOMHelper.getStringValue(jobNode, "ServerID");
		// Obtain the FWAnalyzer information into a separate object.
		Element fwData  = DOMHelper.getElement(jobNode, "FWAnalyzer");
		FWAnalyzer fwa  = FWAnalyzer.create(fwData);
		// Obtain the clustering threshold.
		int threshold   = DOMHelper.getIntValue(jobNode, "Threshold");
		
		// Now we have the information to create a minimal clustering
		// job class. So do it and use it for further operations.
		ClusteringJob job= new ClusteringJob(jobID, srvrID, fwa, threshold);
		
		// Now let the job class of job marshal the non-final values.
		job.unmarshal(jobNode);
		// Return the fully populated entry back to the caller
		return job;
	}
	
	@Override
	protected void unmarshal(Element jobNode) throws Exception {
		// Let base class un-marshal common elements
		super.unmarshal(jobNode);
		
		// Extract the heuristic and filter lists using helper methods.
		heuristics = parseHeuristicChain(jobNode);
		filters    = parseFilterChain(jobNode);
		parameters = parseParameters(jobNode);
	}

	/**
	 * Helper method to unmarshal the heuristic list into an 
	 * in-memory array list.
	 * 
	 * This is a helper method that is invoked from the the create() 
	 * method of a derived class to un-marshal the XML DOM tree 
	 * corresponding to the heuristic list node into
	 * suitable in-memory classes for further processing.
	 * 
	 * @param jobNode The top-level job node from where the 
	 * heuristic chain element is to be extracted and processed. 
	 * 
	 * @return An array list containing a list of Heuristic objects.
	 * If an heuristic chain entry is not found in the DOM tree,
	 * then this method returns null. 
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	private static ArrayList<Heuristic> parseHeuristicChain(Element jobNode) throws Exception {
		// Parse out the heuristic chain using a helper method.
		Element chainNode= DOMHelper.getElement(jobNode, "HeuristicChain");
		if (chainNode == null) {
			return null;
		}
		NodeList chain   = chainNode.getElementsByTagName("Heuristic");
		ArrayList<Heuristic> heuristics = new ArrayList<Heuristic>(); 
		for(int idx = 0; (idx < chain.getLength()); idx++) {
			Element node = (Element) chain.item(idx);
			// Create a heuristic entry using the DOM data
			Heuristic heuristic = Heuristic.create(node);
			heuristics.add(heuristic);
		}
		// Return the parsed-in heuristic chain for further use.
		return heuristics;
	}
	
	/**
	 * Helper method to un-marshal the filter list into an 
	 * in-memory array list.
	 * 
	 * This is a helper method that is invoked from the the create() 
	 * method to un-marshal the XML DOM tree corresponding to the filter
	 * list node into suitable in-memory classes for further processing.
	 * 
	 * @param jobNode The top-level job node from where the filter chain 
	 * element is to be extracted and processed. 
	 * 
	 * @return An array list containing the list of filter objects in the 
	 * filter chain. If a filter chain is not found then this method 
	 * return null. 
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	private static ArrayList<Filter> parseFilterChain(Element jobNode) throws Exception {
		// Parse out the heuristic chain using a helper method.
		Element chainNode= DOMHelper.getElement(jobNode, "FilterChain");
		if (chainNode == null) {
			return null;
		}
		NodeList chain   = chainNode.getElementsByTagName("Filter");
		ArrayList<Filter> filters = new ArrayList<Filter>(); 
		for(int idx = 0; (idx < chain.getLength()); idx++) {
			Element node = (Element) chain.item(idx);
			// Create a heuristic entry using the DOM data
			Filter heuristic = Filter.create(node);
			filters.add(heuristic);
		}
		// Return the parsed-in heuristic chain for further use.
		return filters;
	}
	
    /**
     * Helper method to create a default object representing the default
     * values used by PEACE clustering engine. This method is typically
     * used to create a suitable, default analyzer entry when performing
     * assembly with EAST.
     *
     * @return The newly created analyzer entry based on the default
     * parameters used by the PEACE clustering engine.
     */
	public static  ArrayList<Filter> createDefaultFilters() {
		// Create the length filter
		Filter lenFilter = new Filter(FilterType.LengthFilter);
		lenFilter.addParameter(new Param("minESTLen", "50"));
		// Create the low-complexity filter
		Filter lcFilter = new Filter(FilterType.LCFilter);
		lcFilter.addParameter(new Param("lcPatterns", "A C"));
		// Create the filter chain/list and add the filters to it
		ArrayList<Filter> filters = new ArrayList<Filter>();
		filters.add(lenFilter);
		filters.add(lcFilter);
		return filters;
	}

	/**
	 * Constructor to create a minimally populated clustering job
	 * object with just the fixed value fields initialized to 
	 * specific values.
	 * 
	 * This constructor is typically used to create a temporary job
	 * object into which additional values are going to be loaded
	 * from a XML work space. This constructor is currently used
	 * only by the {@link #create(Element)} method.
	 * 
	 * @param jobID The workspace-wide unique, generated job ID for this job.
	 * The jobID is generated via a call to JobList.reserveJobID() method.
	 * This is a valid string (typically in the form job####)
	 * 
	 * @param serverID The ID of the server on which this job is running.
	 * This is a cross reference ID of a Server configured in this workspace.
	 * 
	 * @param analyzer Information about the distance/similarity analyzer
	 * that was used to create this MST data set.
	 *            
	 * @param threshold The threshold value used to partition a MST.
	 */
	private ClusteringJob(String jobID, String serverID,
			FWAnalyzer analyzer, int threshold) {
		super(JobType.CLUSTERING, jobID, serverID);
		// Save information about analyzer and clustering parameters
		this.analyzer   = analyzer;
		this.threshold  = threshold;
		// Save information about heuristics and filters
		this.heuristics = null;
		this.filters    = null;
	}
	
	/**
	 * Constructor to create a Job object with the fixed value fields
	 * initialized to specific values.
	 * 
	 * @param jobID The workspace-wide unique, generated job ID for this job.
	 * The jobID is generated via a call to JobList.reserveJobID() method.
	 * This is a valid string (typically in the form job####)
	 * 
	 * @param description A user-supplied description for this job entry. This
	 * maybe an empty string (but cannot be null).
	 * 
	 * @param serverID The ID of the server on which this job is running.
	 * This is a cross reference ID of a Server configured in this workspace.
	 * 
	 * @param path The directory on the server where the data for this job is stored.
	 * 
	 * @param nodes The number of compute nodes that were requested on a cluster for
	 * running this job. This value must be at least 1.
	 * 
	 * @param cpusPerNode The number of CPUs on each node that were requested for this
	 * job. This value must be at least 1.
	 * 
	 * @param memory The total memory (sum of all memory used by all processes, 
	 * (in MB) that was allocated for this job.
	 * 
	 * @param maxRunTime The maximum runtime (in hours) that was assigned for 
	 * this job.
	 * 
	 * @param analyzer Information about the distance/similarity analyzer
	 * that was used to create this MST data set.
	 *            
	 * @param threshold The threshold value used to partition a MST.
	 * 
	 * @param heuristics The list of heuristics that were used to accelerate 
	 * the MST generation algorithm. This list may be null.
	 * 
	 * @param filters The list of filters that were used to filter out short ESTs or
	 * ESTs with low complexity sections to ensure that the resultant clustering is
	 * high quality. This list may be null.
	 * 
	 * @param parameters The list of general purpose parameters to be set
	 * for this job. These parameters are used to further customize the
	 * operations of the clustering engine.
	 */
	public ClusteringJob(String jobID, String description, String serverID,
			String path, int nodes, int cpusPerNode,
			int memory, int maxRunTime, FWAnalyzer analyzer, 
			int threshold, ArrayList<Heuristic> heuristics, 
			ArrayList<Filter> filters, ArrayList<Param> parameters) {
		super(JobType.CLUSTERING, jobID, serverID, description, path, 
			nodes, cpusPerNode, memory, maxRunTime);
		// Save information about analyzer and clustering parameters
		this.analyzer   = analyzer;
		this.threshold  = threshold;
		// Save information about heuristics and filters
		this.heuristics = heuristics;
		this.filters    = filters;
		this.parameters = parameters;
	}

	/**
	 * Command line for PEACE clustering engine to configure heuristics.
	 * 
	 * <p>This method can be used to obtain the heuristic information in the form
	 * of command line parameters that can be readily passed to the PEACE
	 * clustering engine. The command line parameters are used to configure the
	 * clustering engine to suit the configuration setup by the user for this
	 * job.</p>
	 * 
	 * <p><b>Note:<b> The return value of this method must be ignored if this
	 * call is being made in conjunction with the two pass d2 analyzer.</p> 
	 * 
	 * @return The command line (to correspondingly setup the heuristics) to
	 * be passed to the PEACE clustering engine. If no heuristics are 
	 * defined then this method returns an empty string. If an empty list
	 * of heuristics are specified then this method returns
	 * "--heuristics null" as the result.
	 */
	@Override
	public String getHeuristicsCmdLine() {
		if (heuristics == null) {
			return ""; // no heuristics at all.
		}
		String cmdLine = "";
		if (heuristics != null) {
			// Convert heuristic information to command line parameters.
			String heuristicParams = "";
			cmdLine += "--heuristics ";
			if (heuristics.size() > 0) {
				for(int i = 0; (i < heuristics.size()); i++) {
					Heuristic heur = heuristics.get(i);
					cmdLine += (i > 0 ? " " : "") + heur.getName();
					heuristicParams += heur.toCmdLine();
				}
			} else {
				// The command line parameter for no heuristics case is null
				cmdLine += "null";
			}
			cmdLine += heuristicParams;
		}
		// Return the cmd line for the heuristics
		return cmdLine;
	}

	/**
	 * Command line for PEACE clustering engine to configure filters.
	 * 
	 * This method can be used to obtain the filters information in the form
	 * of command line parameters that can be readily passed to the PEACE
	 * clustering engine. The command line parameters are used to configure the
	 * filters used by the clustering engine to mirror the configuration setup 
	 * by the user for this job.
	 * 
	 * @return The command line (to correspondingly setup the filters) to
	 * be passed to the PEACE clustering engine. If no filters are 
	 * defined then this method returns an empty string. If an empty list
	 * of filters are specified then this method returns
	 * "--filters null" as the result.
	 */
	public String getFiltersCmdLine() {
		if (filters == null) {
			return ""; // no filters at all
		}
		// Next covert filter information to command line parameters.
		String filterParams = "";
		String cmdLine = "--filters ";
		if (filters.size() > 0) {
			for(int i = 0; (i < filters.size()); i++) {
				Filter filter = filters.get(i);
				cmdLine += (i > 0 ? " " : "") + filter.getName();
				filterParams += filter.toCmdLine();
			}
		} else {
			// The command line parameter for no heuristics case is null
			cmdLine += "null";
		}
		cmdLine += filterParams;
		// Return the cmd line.
		return cmdLine;
	}
	
	/**
	 * Return the information in the form of a partial PEACE command line 
	 * parameter.
	 * 
	 * This method can be used to obtain the information needed to
	 * generate the MST and clusters based on the supplied information in the form
	 * of a command line parameter. 
	 * 
	 * @return Return the information as a command line to be passed to the
	 * PEACE clustering engine.
	 */
	public String toCmdLine() {
		// Use helper method to build the command line.
		return " " + getHeuristicsCmdLine() + " " + getFiltersCmdLine() +
			" " + analyzer.toCmdLine() + " " + toCmdLine(parameters);
	}
	
	@Override
	public String toString() {
		return jobID;
	}
	
	/**
	 * Obtain the list of general purpose parameters associated
	 * with this job.
	 * 
	 * This is a helper method creates the list of general purpose
	 * command-line parameters associated with this job. These
	 * parameters are not associated with filters or heuristics. However,
	 * they are used to further customize the operations of the 
	 * clustering engine.
	 * 
	 * @return An array list containing the set of parameters to
	 * be passed to the parallel clustering engine for this job.
	 * The return value is never null.
	 */
	public ArrayList<Param> getParameters() {
		return parameters;
	}
	
	/**
	 * Helper method to build the set of general purpose command-line
	 * parameters associated with this job.
	 * 
	 * This is a helper method creates the list of general purpose
	 * command-line parameters associated with this job. These
	 * parameters are not associated with filters or heuristics. However,
	 * they are used to further customize the operations of the 
	 * clustering engine.
	 */
	public void setupParameters(DataSet ds, GeneratedFileList gfl,
			boolean maskBases) {
		// Add command-line for input data set file type. Note that
		// here we map the path from local machine to remote machine
		String fileType = (ds.getFileType().equals(DataFileType.SFF) ?
				"--sffFile" : "--fastaFile");
		parameters.add(new Param(fileType, getServerESTFile(ds)));
		// Check and suppress base masking.
		if (!maskBases) {
			// Disable base masking. That is atcg is
			// equal to ATCG.
			parameters.add(new Param("--no-mask-bases", null));
		}
		// Cut the file names with appropriate parameter entries.
		FileEntry mstFE = gfl.findEntry(FileEntry.FileEntryType.MST);
		parameters.add(new Param("--output-mst-file", mstFE.getName()));
		// Add entry for clustering output in GUI compatible format
		parameters.add(new Param("--gui-print", null));
		FileEntry clsFE = gfl.findEntry(FileEntry.FileEntryType.CLS);
		parameters.add(new Param("--output-cls-file", clsFE.getName()));
		// Add command line option to generate progress information
		parameters.add(new Param("--progress", "progress.dat"));
	}
	
	/**
	 * Obtain the list of heuristics that were used to accelerate the
	 * process of constructing the MST. Specifically many of these
	 * heuristics accelerate the frame-word analyzer that was used to
	 * build the MST.
	 * 
	 * @return The list of heuristics that were used to accelerate the
	 * MST construction process. The return value can be null for
	 * certain types of jobs.
	 */
	public ArrayList<Heuristic> getHeuristicList() { return heuristics; }
	
	/**
	 * Obtain the list of filters that were used to improve quality of
	 * clustering. Specifically the filters weed out ESTs that are known
	 * to interfere with clustering and deteriorate overall quality of 
	 * results.
	 * 
	 * @return The list of filters hat were used to improve quality of
	 * clustering. The return value can be null for certain types of
	 * jobs.
	 */
	public ArrayList<Filter> getFilterList() { return filters; }
	
	/**
	 * Method to marshal the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent JobList node in the DOM tree.
	 * 
	 * @param jobList The DOM element corresponding to the "JobList"
	 * node that contains this entry.
	 */
	public final void marshall(Element jobList) {
		// Create a top-level entry for this "Job"
		Element job = DOMHelper.addElement(jobList, "ClusteringJob", null);
		// Let base class add common sub-elements
		super.marshall(job);
		// Marshal out the analyzer information
		analyzer.marshall(job);
		// Marshal out the clustering threshold
		DOMHelper.addElement(job, "Threshold", "" + threshold);

		// Now marshal the information regarding heuristics.
		if (heuristics != null) {
			Element chain = DOMHelper.addElement(job, "HeuristicChain", null);
			for(Heuristic heuristic : heuristics) {
				heuristic.marshall(chain);
			}
		}
		// Now marshal the information regarding filters.
		if (filters != null) {
			Element chain = DOMHelper.addElement(job, "FilterChain", null);
			for(Filter filter : filters) {
				filter.marshall(chain);
			}
		}
		// Now marshal the information regarding general parameters
		marshallParameters(parameters, job);
	}
	
	/**
	 * Method to marshal the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t";
		final String NUM_ELEMENT = Indent + "\t" + "<%1$s>%2$d</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<ClusteringJob>\n", Indent);
		// Let base class write out common elements first
		super.marshall(out);
		// Marshal out the analyzer information
		analyzer.marshall(out);
		// Marshal out the clustering threshold
		out.printf(NUM_ELEMENT, "Threshold", threshold);
		
		// Marshal out the heuristic chain (if present)
		if (heuristics != null) {
			out.printf("%s\t<HeuristicChain>\n", Indent);
			for(Heuristic heuristic : heuristics) {
				heuristic.marshall(out);
			}
			out.printf("%s\t</HeuristicChain>\n", Indent);
		}
		// marshal out the filter chain.
		if (filters != null) {
			out.printf("%s\t<FilterChain>\n", Indent);
			for(Filter filter : filters) {
				filter.marshall(out);
			}
			out.printf("%s\t</FilterChain>\n", Indent);
		}
		// marshal out any parameters we may have
		marshallParameters(parameters, out);
		// Close the job tag
		out.printf("%s</ClusteringJob>\n", Indent);
	}

	/**
	 * Method to write summary information about the job.
	 * 
	 * This method is a convenience method that is used by various 
	 * wizards to display summary information about the job.
	 * The summary information about the job include information
	 * about the filter, heuristics, server-configuration etc.
	 * 
	 * @param sw The summary writer to which the data is to be written.
	 * 
	 */
	@Override
	public void summarize(SummaryWriter sw) {
		// Cut summary information about filters (if applicable)
		if (filters != null) {
			sw.addSection("Filter summary");
			for(Filter f: filters) {
				f.summarize(sw);
			}
		}
		if (heuristics != null) {
			// Cut summary information about heuristic
			sw.addSection("Heuristic summary");
			if ((analyzer != null) && (FWAnalyzerType.TWOPASSD2.equals(analyzer.getType()))) {
				sw.addSummary("Automatic", "N/A", 
				"Automatic and dynamically configured u/v & t/v heuristics");
			} else {
				// Let each heuristic log summary
				for(Heuristic heur: heuristics) {
					heur.summarize(sw);
				}
			}
		}
		// Let the base class summarize common information
		super.summarize(sw);
	}
	
	/**
	 * Determine if this job uses automatically configured heuristics.
	 * 
	 * Some analyzers such as the two-pass-d2 analyzers are adaptive
	 * and they use a custom configured u/v and t/v heuristics.
	 * For such analyzers the heuristic chains are not meaningful
	 * and are typically not displayed.
	 *  
	 * @return This method returns true if this job relies on an
	 * analyzer that uses automatically configured heuristics.
	 * Otherwise this method returns false.
	 */
	public boolean isAutoHeuristic() {
		return FWAnalyzerType.TWOPASSD2.equals(analyzer.getType());
	}
	
	/**
	 * Obtain meta data about the frame/word analyzer that was used to generate
	 * metrics (may it be similarity or distance values) that were used to 
	 * build the MST.
	 * 
	 * @return This method returns information about the analyzer that was
	 * used to generate the MST. 
	 */
	public FWAnalyzer getFWAnalyzer() { return analyzer; }
	
	/**
	 * Information about the distance/similarity analyzer that was used
	 * to create this MST data set.
	 */
	private final FWAnalyzer analyzer;
	
	/**
	 * The threshold value used for cluster generating. The threshold defines
	 * the metric on a given edge of a MST that identifies transition to a 
	 * new cluster.
	 */
	private final int threshold;
	
	/**
	 * The list of heuristics configured by the user for this job. These
	 * heuristics are used to accelerate the MST generation algorithm.
	 * This array list can be null (depending on the job type).
	 */
	private ArrayList<Heuristic> heuristics;
	
	/**
	 * The list of filters configured by the user for this job. These
	 * filters are used to weed out fragments that may interfere with
	 * clustering and reduce overall quality of clustering.
	 * This array list can be null (depending on the job type).
	 */
	private ArrayList<Filter> filters;
	
	/**
	 * The list of parameters configured by the user for this job. These
	 * parameters are passed to corresponding assembly tools.
	 * This array list cannot be null (but can be empty).
	 */
	private ArrayList<Param> parameters;
}
