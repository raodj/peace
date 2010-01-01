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

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * <p>This class encapsulates all the data necessary to use a given Minimum
 * Spanning Tree (MST) data file. A MST data file essentially contains the MST
 * structure with closely related ESTs on the same subtree of the MST.
 * Constructing the MST is the most computationally intensive part of the
 * clustering algorithm supported by PEACE. This element contains meta data
 * about the MST file. The MSTData class is created each time a new Job is run
 * to build a MST on a given EST data set. The EST data is the primary input
 * from which the MST data is derived.</p>
 * 
 * <p>In addition to encapsulating the data this class also provides
 * convenient interfaces for marshalling and unmarshalling XML data compatible
 * with the PEACE GUI configuration XML file format. </p>
 */
public class MSTData {
	/**
	 * Different enumerations defining the different approaches that PEACE
	 * currently supports for building a MST.
	 */
	public enum MSTBuilderType {
		/**
		 * This is the default and standard MST building algorithm that
		 * is based on Prim's algorithm.
		 */
		MST,
		/**
		 * This method builds an approximately close MST -- it is not the 
		 * MST but is a tree that is close enough to MST.
		 */
		PMST,
		/**
		 * This approach uses the concept of transitivity to build a MST.
		 * Seems like this approach has some promise but has some 
		 * implementation related performance issues that need to be
		 * ironed out.
		 */
		TMST
	};
	
	/**
	 * Helper method to utilize data from a DOM tree to create a suitable
	 * MSTData entry. This method is typically  used to create a suitable
	 * entry when loading a Work space into the GUI.
	 * 
	 * @param mstData The DOM element to be used for creating the MSTData
	 * entry and populating with the needed data.
	 * 
	 * @return The newly created data entry based on the DOM data.
	 * 
	 * @throws Exception This method throws an exception when errors occur
	 * during reading and processing elements from the DOM node.
	 */
	public static MSTData create(Element mstData) throws Exception {
		// First extract the necessary information from the DOM tree.
		String id      = DOMHelper.getStringValue(mstData, "ID");
		String typeStr = DOMHelper.getStringValue(mstData, "Type");
		String path    = DOMHelper.getStringValue(mstData, "Path");
		String desc    = DOMHelper.getStringValue(mstData, "Description", true);
		// Obtain the FWAnalyzer information into a separate object.
		Element fwData  = DOMHelper.getElement(mstData, "FWAnalyzer");
		FWAnalyzer analy= FWAnalyzer.create(fwData);
		// Parse out the heuristic chain using a helper method.
		Element chainNode= DOMHelper.getElement(mstData, "HeuristicChain");
		NodeList chain   = chainNode.getElementsByTagName("Heuristic");
		ArrayList<Heuristic> heuristics = new ArrayList<Heuristic>(); 
		for(int idx = 0; (idx < chain.getLength()); idx++) {
			Element node = (Element) chain.item(idx);
			// Create a heuristic entry using the DOM data
			Heuristic heuristic = Heuristic.create(node);
			heuristics.add(heuristic);
		}
		// load up the job summary information.
		Element jobData    = DOMHelper.getElement(mstData, "JobSummary");
		JobSummary summary = JobSummary.create(jobData);
		// Now that we have all the information create the actual 
		// MSTData.
		typeStr = typeStr.toUpperCase();
		MSTBuilderType type = MSTBuilderType.valueOf(MSTBuilderType.class, typeStr);
		return new MSTData(id, type, path, desc, analy, heuristics, summary);
	}
	
	/**
	 * The constructor to create a fully populated MSTData object that contains
	 * all the meta data regarding a MST file that has been generated from a
	 * given EST file.
	 * 
	 * @param id
	 *            The unique generated data set ID value for this data.
	 * @param type
	 *            he type of MST building approach supported by PEACE that was
	 *            used to obtain the MST referred by this class.
	 * @param path
	 *            The path to the actual MST file (on the local machine) that is
	 *            referred by this entry.
	 * @param description
	 *            A user defined description for this file.
	 * @param analyzer
	 *            Information about the distance/similarity analyzer that was
	 *            used to create this MST data set.
	 * @param heuristics
	 *            The list of heuristics that were used to accelerate the MST
	 *            generation algorithm.
	 * @param summary
	 *            The core/useful information about the job that was run to
	 *            compute the MST.
	 */
	public MSTData(String id, MSTBuilderType type, String path,
			String description, FWAnalyzer analyzer, 
			ArrayList<Heuristic> heuristics, JobSummary summary) {
		this.id          = id;
		this.type        = type;
		this.path        = path;
		this.description = description;
		this.analyzer    = analyzer;
		this.heuristics  = heuristics;
		this.jobSummary  = summary;
	}
	
	/**
	 * Obtain the work space wide unique identifier set for this data set. The
	 * ID value is created when new jobs are scheduled to create MST data. These
	 * IDs are persisted in the work space configuration file and loaded when a
	 * work space is opened in the GUI.
	 * 
	 * @return This method returns the unique identifier set for this data set.
	 */
	public String getID() { return id; }
	
	/**
	 * The type of MST building approach supported by PEACE that was used to
	 * obtain the MST referred by this class.
	 * 
	 * @return The type of approach used to build the MST.
	 */
	public MSTBuilderType getType() { return type; }
	
	/**
	 * Obtain the complete file name and path where the actual MST file is
	 * located. This value is set when a new job to generate this MST data is
	 * scheduled. The path for the job is persisted in the work space
	 * configuration file and loaded when a work space is opened in the GUI.
	 * 
	 * @return This method returns the path to the MST file.
	 */
	public String getPath() { return path; }
	
	/**
	 * Obtain the user supplied description for this this data. The value is set
	 * when new jobs are scheduled. The description is persisted in the work
	 * space configuration file and loaded when a work space is opened in the
	 * GUI.
	 * 
	 * @return This method returns the user supplied description set for this
	 *         job.
	 */
	public String getDescription() { return description; }

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
	 * Obtain meta data about the job that was run to generate the MST.
	 * The JobID in the job summary can be used to look up additional
	 * information about the job. Similarly, the serverID value in the
	 * job summary can be used to look up information about the server
	 * on which the job was run.
	 * 
	 * @return This method returns information about the job that was
	 * run to generate the MST. 
	 */
	public JobSummary getJobSummary() { return jobSummary; }

	/**
	 * Obtain the list of heuristics that were used to accelerate the
	 * process of constructing the MST. Specifically many of these
	 * heuristics accelerate the frame-word analyzer that was used to
	 * build the MST.
	 * 
	 * @return The list of heuristics that were used to acclerate the
	 * MST construction process.
	 */
	public ArrayList<Heuristic> getHeuristicList() { return heuristics; }
	
	/**
	 * Method to marshall the data stored in this object to become part of
	 * a DOM tree element passed in. This method assumes that the element
	 * passed in corresponds to the parent DataSet node in the DOM tree.
	 * 
	 * @param dataset The DOM element corresponding to the "DataSet"
	 * node that contains this entry.
	 */
	public final void marshall(Element dataset) {
		// Create a top-level entry for this "Job"
		Element mstData = DOMHelper.addElement(dataset, "MSTData", null);
		// Add new sub-elements for each sub-element
		DOMHelper.addElement(mstData, "ID", id);
		DOMHelper.addElement(mstData, "Type", type.toString().toLowerCase());
		DOMHelper.addElement(mstData, "Path", path);
		DOMHelper.addElement(mstData, "Description", 
				(description != null) ? description : "");
		// Get the fw analyzer add its own information.
		analyzer.marshall(mstData);
		// Now marshall the information regarding heuristics.
		Element chain = DOMHelper.addElement(mstData, "HeuristicChain", null);
		for(Heuristic heuristic : heuristics) {
			heuristic.marshall(chain);
		}
		// Finally marshall the job summary information.
		jobSummary.marshall(mstData);
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a
	 * XML fragment. The XML fragment is guaranteed to be compatible
	 * with the PEACE work space configuration data. 
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		final String Indent = "\t\t";
		final String STR_ELEMENT = Indent + "\t" + "<%1$s>%2$s</%1$s>\n";
		
		// Create a top-level server entry for this server
		out.printf("%s<MSTData>\n", Indent); 
		// Add new sub-elements for each value.
		out.printf(STR_ELEMENT, "ID", id);
		out.printf(STR_ELEMENT, "Type", type.toString().toLowerCase());
		out.printf(STR_ELEMENT, "Path", path);
		out.printf(STR_ELEMENT, "Description", 
				(description != null) ? DOMHelper.xmlEncode(description) : "");
		// Get the f/w analyzer add its own information.
		analyzer.marshall(out);
		// marhsall out the heuristic chain.
		out.printf("%s\t<HeuristicChain>\n", Indent);
		for(Heuristic heuristic : heuristics) {
			heuristic.marshall(out);
		}
		out.printf("%s\t</HeuristicChain>\n", Indent);
		// Marshall the job summary out.
		jobSummary.marshall(out);
		// Close the MSTData tag
		out.printf("%s</MSTData>\n", Indent);
	}
	
	/**
	 * Provides a textual multi-line summary for this element.
	 * 
	 * This method is typically used by GUI to display a simple 
	 * representation of the data stored in this element.
	 * 
	 * @param indent The leading spaces to be used for indenting each
	 * line. This parameter cannnot be null.
	 * 
	 * @return A textual multi-line summary of the data in this element.
	 */
	public String getSummary(String indent) {
		StringBuilder sb = new StringBuilder(256);
		// First append standard information.
		sb.append(indent + "ID: " + id + "\n");
		sb.append(indent + "MST Builder: " + type.toString() + "\n");
		sb.append(indent + "Path: " + path + "\n");
		sb.append(indent + "Description: " + description + "\n");
		// Next append the heuristics information.
		sb.append(indent + "Heuristics: " + "\n");
		for(Heuristic heur: heuristics) {
			sb.append(indent + indent + heur + "\n");
		}
		return sb.toString();
	}
	
	/**
	 * Return the information in the form of a partial PEACE command line 
	 * parameter.
	 * 
	 * This method can be used to obtain the information needed to
	 * generate the MST based on the supplied information in the form
	 * of a command line parameter. 
	 * 
	 * @return Return the information as a command line parameter.
	 */
	public String toCmdLine() {
		String cmdLine = "--clusterMaker " + type.toString().toLowerCase();
		cmdLine += " " + analyzer.toCmdLine();
		// Convert heuristic information to command line parameters.
		String heuristicParams = "";
		cmdLine += " --heuristics ";
		for(int i = 0; (i < heuristics.size()); i++) {
			Heuristic heur = heuristics.get(i);
			cmdLine += (i > 0 ? "-" : "") + heur.getName();
			heuristicParams += heur.toCmdLine();
		}
		cmdLine += heuristicParams;
		// Add the mst file name (without path)
		File mstFile = new File(path);
		cmdLine += " --output-mst-file " + mstFile.getName();
		// Return the cmd line.
		return cmdLine;
	}
	
	/**
	 * Overrides the default implementation in the base class to 
	 * simply return the last part of the MST file associated
	 * with this entry.
	 * 
	 * @return A short string representation that is easy to display
	 * primarily in a tree view of the work space.
	 */
	@Override
	public String toString() {
		File tmpData = new File(path);
		return "MST Data [File: " + tmpData.getName() + "]";
	}
	
    /**
     * Set the data set that this object logically belongs to. 
     * 
     * <p><b>Note:</b>  This method is used by the DataSet. Use 
     * {@link DataSet#add(MSTData)} method to add a MSTData object 
     * to a data set.</p>
     * 
     * @param dataSet The data set to which this object has been added.
     */
    protected void setDataSet(DataSet dataSet) {
            this.dataSet = dataSet;
    }
    
    /**
     * Obtain the data set that logically contains this MST Data.
     * 
     * @return The data set that logically contains this object. This value
     * is set only after this object has been added to a data set.
     */
    public DataSet getDataSet() { return dataSet; }
    
    /**
     * Set revised updated information about the job associated with
     * this MSTData entry.
     * 
     * @param job The job from where the necessary information is
     * to be copied.
     */
    public void updateJobSummary(Job job) {
    	jobSummary = new JobSummary(job);
   		// Notify all listeners about the change.
    	Workspace ws = Workspace.get();
		ws.fireWorkspaceChanged(new WorkspaceEvent(this, 
				WorkspaceEvent.Operation.INSERT));
    }
    
    /**
     * Reference to the data set that logically contains this MSTData 
     * element. This element is not serialized to XML format 
     * (and consequently is flagged transient).
     */
    private DataSet dataSet;

	/**
	 * The unique generated data set ID value for this data. This value is
	 * typically created when a new job is scheduled. This value is persisted in
	 * the work space configuration.
	 */
	protected final String id;
	
	/**
	 * The type of MST building approach supported by PEACE that was used to
	 * obtain the MST referred by this class.
	 */
	private final MSTBuilderType type;

	/**
	 * The path to the actual MST file (on the local machine) that is 
	 * referred by this entry.
	 */
	private final String path;

	/**
	 * A user defined description for this file. This description is set when
	 * a new job is scheduled and persisted in the work space configuration for
	 * future references.
	 */
	private final String description;
	
	/**
	 * Information about the distance/similarity analyzer that was used
	 * to create this MST data set.
	 */
	private final FWAnalyzer analyzer;
	
	/**
	 * The list of heuristics that were used to accelerate the MST generation
	 * algorithm.
	 */
	ArrayList<Heuristic> heuristics;
	
	/**
	 * The core/useful information about the job that was run to compute
	 * the MST.
	 */
	private JobSummary jobSummary;
}
