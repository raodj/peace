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
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.GregorianCalendar;

import javax.xml.XMLConstants;
import javax.xml.datatype.DatatypeFactory;
import javax.xml.datatype.XMLGregorianCalendar;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Result;
import javax.xml.transform.Source;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;
import javax.xml.validation.Validator;

import org.peace_tools.generic.Utilities;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * <p>
 * The top-level Workspace class that stores and manages all the data associated
 * with a PEACE Workspace. Note that the purpose of the Workspace class is to
 * encapsulate and manage the data associated with the workspace and not to
 * display it. Display of Workspace data is delegated to other classes that
 * focus on presenting a suitable "view" of the data. This design follows the
 * Model-View-Controller (MVC) design pattern.
 * </p>
 * 
 * Currently, each PEACE GUI instance operates only with a single Workspace.
 * Consequently to enforce the use of a single, unique Workspace object, this
 * class has been designed using the Singleton design pattern.
 */
public class Workspace {
	/**
	 * Method to get a reference to the process-wide, unique instance of
	 * the PEACE workspace. The PEACE workspace object is a singleton (
	 * there is only one object). Making changes in one affects all of
	 * the views of the workspace.
	 * 
	 * @return The process-wide, unique instance of PEACE workspace, if
	 * one is available. If a workspace has not been created, then this
	 * method returns null.
	 * 
	 * @see Workspace#useWorkspace(String)
	 * @see Workspace#createDefault(String)
	 */
	public static Workspace get() {
		return workspace;
	}

	/**
	 * Method to create default work space met data in a given directory. 
	 * This method assumes that the file named "PEACEWorkspace.xml"  in the
	 * given directory can be overwritten if it exists. The created workspace
	 * does not have any data sets, servers, or jobs in it.  If errors occur
	 * during the creation the existing workspace remains untouched. On
	 * successful creation, the global workspace is changed to the newly
	 * created workspace (this tantamounts to calling useWorkspace method)
	 * 
	 * @param directory
	 *            The new workspace directory where an empty work space
	 *            needs to be created.
	 * 
	 * @throws Exception
	 *             This method generates exception if the workspace file could
	 *             not be created.
	 *             
	 * @return The reference to the newly created, global workspace. This
	 * will be the same object returned by the Workspace.get() method once
	 * this call is successful.
	 */
	public static Workspace createDefault(final String directory) throws Exception {
		// Create the actual file name to load.
		final String filename = getWorkspaceFile(directory);
		// Create the output stream to write initial empty meta data
		PrintWriter  outStream = new PrintWriter(new FileWriter(filename));
		
		// Now create an empty work space
		Workspace temp = new Workspace();
		// Setup the core data elements.
		GregorianCalendar now   = new GregorianCalendar();
		temp.workspaceDirectory = directory;
		DatatypeFactory codec   = DatatypeFactory.newInstance();
		temp.creationTimestamp  = codec.newXMLGregorianCalendar(now);
		temp.seqCounter         = 1;
		
		// Serialize the XML data out to the stream.
		temp.marshall(outStream);
		// Close the file.
		outStream.close();
		// Update to the new workspace
		workspace = temp;
		// Return a handy reference to the global workspace
		return workspace;
	}

	/**
	 * Method to load workspace data from a given directory. This method assumes
	 * that the given directory contains the "PEACEWorkspace.xml" file, loads
	 * the workspace file, and parses the data into convenient class hierarchy
	 * encapsulated by this class.
	 * 
	 * @param directory
	 *            The new workspace directory from where the necessary data is
	 *            to be loaded.
	 * 
	 * @throws Exception
	 *             This method generates exception if the workspace file does
	 *             not exist or the data in the workspace file is
	 *             corrupted/inconsistent with expectations
	 */
	public static Workspace useWorkspace(final String directory) throws Exception {
		// Create the actual file name to load.
		final String filename = getWorkspaceFile(directory);
		// Obtain a DOM parser from via suitable document builder factory.
		DocumentBuilderFactory fact = DocumentBuilderFactory.newInstance();
		fact.setNamespaceAware(true);
		DocumentBuilder builder = fact.newDocumentBuilder();
		// Create a file object and have the DOM builder build a DOM tree
		// for us. We would have done this with SAX API too. However, it
		// felt like DOM was going to be easier. Second, the Workspace
		// data is not not very big.
		File inputFile = new File(filename);
		Document workspaceData = builder.parse(inputFile);
		// Validate the document
		validate(workspaceData);
		// First load the standard top-level elements from the workspaceData
		// and create a new workspace object.
		workspace = loadCoreWorkspaceData(workspaceData, directory);
		// First serialize all the data sets in this work space.
		NodeList dataSets = workspaceData.getElementsByTagName("DataSet");
		for(int idx = 0; (idx < dataSets.getLength()); idx++) {
			Element dataSetNode = (Element) dataSets.item(idx);
			DataSet dataSet = DataSet.create(dataSetNode);
			workspace.dataSetList.add(dataSet);
		}
		// OK, next serialize the serverList. First set the appropriate element.
		Element serverListData = DOMHelper.getElement(workspaceData, "ServerList");	
		workspace.serverList = ServerList.create(serverListData);
		// Next serialize the job List. First set the appropriate element.
		Element jobListData = DOMHelper.getElement(workspaceData, "JobList");	
		workspace.jobList = JobList.create(jobListData);
		// Lastly serialize the DB classifier list. First set the appropriate element.
		Element dbClasListData = DOMHelper.getElement(workspaceData, "ClassifierList");	
		workspace.classifierList = ClassifierList.create(dbClasListData);
		// Return a handy reference to the global workspace
		return workspace;
	}

	/**
	 * Method to provide consistent name for the PEACE workspace file name
	 * so that the file names are consistently used throughout.
	 * 
	 * @param directory
	 *            The workspace directory from where the necessary data is
	 *            to be loaded.
	 */
	public static String getWorkspaceFile(final String directory) {
		return directory + "/PEACEWorkspace.xml";
	}
	
	/**
	 * This is a helper method that can be used to validate a Workspace DOM
	 * tree against the PEACE schema. If the DOM document is valid, this
	 * method silently returns. Otherwise it throws an exception to signal
	 * errors.
	 * 
	 * @param workspace The DOM document to be verified.
	 * 
	 * @throws Exception This method generates exceptions if the PEACE.xsd
	 * file could not be loaded or if the DOM document is invalid.
	 */
	protected static void validate(Document workspace) throws Exception {
		// Load the schema file into a string.
		InputStream is = Utilities.getStream("installFiles/PEACE.xsd");
		// Obtain schema factory to create a validator.
	    final SchemaFactory factory = 
	    	SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
	    // Load the PEACE schema data (from string) to create a in-memory schema.
	    Schema peaceSchema = factory.newSchema(new StreamSource(is));
	    // Crate a schema validator from the schema
	    Validator validator = peaceSchema.newValidator();
	    // Get the validator to verify the DOM document is valid.
	    validator.validate(new DOMSource(workspace));
	}
	
	/**
	 * Obtain a new DOM document representing the data associated with 
	 * the current work space configuration information. This method
	 * creates a new DOM document and populates it with all the configuration
	 * information associated with the workspace. The DOM tree can then be
	 * used for various transformations.
	 * 
	 * @return A newly created DOM Document that contains all the configuration
	 * data associated with this workspace.
	 * 
	 * @throws Exception This method throws a suitable exception if errors
	 * occur during the DOM creation process.
	 */
	public Document marshall() throws Exception {
		// First marshal all the work space information from various
		// helper classes into a single DOM tree. For that first we
		// create a DOM document via a suitable builder class.
		DocumentBuilderFactory fact = DocumentBuilderFactory.newInstance();
		fact.setNamespaceAware(true);
		DocumentBuilder builder = fact.newDocumentBuilder();
		Document workspaceData = builder.newDocument();
		// Add the main root element to the DOMTree.
		Element workspace = workspaceData.createElementNS("http://www.peace-tools.org/", "Workspace");
		workspaceData.appendChild(workspace);
		// Add the necessary namespace attributes to the top-level element.
		workspace.setAttributeNS(null, "Version", "0.2");
		workspace.setAttributeNS("http://www.w3.org/2000/xmlns/", "xmlns", DOMHelper.PEACE_NS);
		workspace.setAttributeNS("http://www.w3.org/2000/xmlns/", "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
		workspace.setAttributeNS("http://www.w3.org/2001/XMLSchema-instance", "xsi:schemaLocation", DOMHelper.PEACE_NS + " PEACE.xsd");
		// Now add the top-level non-list items.
		DOMHelper.addElement(workspace, "Directory", workspaceDirectory);
		DOMHelper.addElement(workspace, "CreationTimestamp", creationTimestamp.toString());
		DOMHelper.addElement(workspace, "SeqCounter", "" + seqCounter);
		// First serialize all the data sets in this workspace.
		for(DataSet data: dataSetList) {
			data.marshall(workspace);
		}
		// Next serialize all the server entries
		serverList.marshall(workspace);
		// Serialize all the jobs in this work space.
		jobList.marshall(workspace);
		// Finally marshall all the data base classifiers
		classifierList.marshall(workspace);
		// Return the DOM document created.
		return workspaceData;
	}
	
	/**
	 * Method to marshall the data stored in this object directly to a XML
	 * fragment. The XML fragment is guaranteed to be compatible with the
	 * PEACE work space configuration data. This method provides better control
	 * on the XML formatting to generate a more readable XML output when
	 * compared to the DOM tree.
	 * 
	 * @param out The stream to which the XML must be serialized.
	 */
	public final void marshall(PrintWriter out) {
		// Generate the top-level element with name space information.
		out.println("<?xml version=\"1.0\"?>");
		// Ideally, th common sub-strings used the two marshall methods must be
		// defined as static constants and used.
		out.println("<Workspace xmlns=\"http://www.peace-tools.org/\"\n" +
						"\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n" +
						"\txsi:schemaLocation=\"http://www.peace-tools.org/ PEACE.xsd\"\n" +
						"\tVersion=\"0.2\">\n");
		// Create a top-level server list entry for this class
		out.printf("\t<%1$s>%2$s</%1$s>\n",   "Directory", workspaceDirectory);
		out.printf("\t<%1$s>%2$s</%1$s>\n\n", "CreationTimestamp", creationTimestamp.toString());
		out.printf("\t<%1$s>%2$d</%1$s>\n\n", "SeqCounter", seqCounter);
		// First serialize all the data sets in this workspace.
		for(DataSet data: dataSetList) {
			data.marshall(out);
		}
		// Marshall the xml data for server list to the print writer
		serverList.marshall(out);
		// Next serialize all the jobs in this workspace.
		jobList.marshall(out, "");
		// Finally serialize all the data base classifiers
		classifierList.marshall(out);
		// Close the works pace element.
		out.println("</Workspace>");
	}
	
	/**
	 * This method can be used to save the work space configuration to disk. This
	 * method (unlike the saveWorkspace() method) uses a DOM tree to save the
	 * data. The DOM tree containing the work space configuration information is
	 * obtained from via the marhsal method.
	 *  
	 * @throws Exception This method throws a suitable exception if errors
	 * occur during the file saving process.
	 */
	public void saveWorkspaceViaDOM() throws Exception {
		// First marshall the configuration into a DOM document.
		Document workspace = marshall();
		// Validate the document
		validate(workspace);
		// Create a transformer object to actually write the DOM
		TransformerFactory tranFact = TransformerFactory.newInstance(); 
        Transformer transformer = tranFact.newTransformer();
        // Set properties on the transformer to make the generated DOM
        // to be more readable.
        // transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
        // Finally write the marshaled xml to disk using the transformer
        Source src  = new DOMSource(workspace); 
        Result dest = new StreamResult(System.out);
        transformer.transform(src, dest);
	}
	
	/**
	 * Save the current work space configuration information. This method
	 * saves the current work space configuration in a suitable XML format.
	 * The default work space directory and file name are used. This method
	 * does not create a DOM tree. Instead XML data is directly generated.
	 * This basically ensures that the XML is formatted to look a bit more
	 * organized and readable than otherwise.
	 * 
	 * @throws Exception This method throws a suitable exception if errors
	 * occur during the file saving process.
	 */
	public synchronized void saveWorkspace() throws Exception {
		// First marshal all the work space information from various
		// sub classes into a stream.
		StringWriter strStream = new StringWriter(1000);
		PrintWriter  outStream = new PrintWriter(strStream);
		// Serialize the XML data out to the stream.
		marshall(outStream);
		// Now write the data to the workspace file.
		final String filename = getWorkspaceFile(getDirectory());
		FileOutputStream fos  = new FileOutputStream(filename);
		fos.write(strStream.toString().getBytes());
		fos.close();
	}
	
	/**
	 * Helper method to simply load the top-level Directory for the Workspace
	 * and the creation time stamp for this workspace from a DOM document.
	 * 
	 * @param workspaceData
	 *            The DOM document from where the workspace information is to be
	 *            retrieved and stored in this class.
	 * 
	 * @return This method creates a new workspace object, populates it with
	 * the core information and returns it.
	 * 
	 * @throws Exception
	 *             This method generates exceptions if errors occur when reading
	 *             the workspace data.
	 */
	private static Workspace loadCoreWorkspaceData(final Document workspaceData, 
			final String directory) throws Exception {
		// The core workspace path and creation time from the DOM tree.
		String workspacePath = DOMHelper.getStringValue(workspaceData, "Directory");
		String timestamp     = DOMHelper.getStringValue(workspaceData,
				"CreationTimestamp");
		// First extract the sequence counter information as string
        String seqCounter = DOMHelper.getStringValue(workspaceData, "SeqCounter");
		// Double check to ensure that the workspace path in XML configuration
		// matches the directory we are loading the data from.
		if (!directory.equals(workspacePath)) {
			throw new IllegalArgumentException("The workspace file in '" +
					directory + "' is referring to a different workspace ('" +
					workspacePath + "')");
		}
		// OK, basic checks have passed. Create new workspace and set
		// the core information.
		workspace = new Workspace();
		workspace.workspaceDirectory= workspacePath;
		DatatypeFactory codec       = DatatypeFactory.newInstance();
		workspace.creationTimestamp = codec.newXMLGregorianCalendar(timestamp);
        // Convert the sequence number from string to integer
        workspace.seqCounter        = Integer.parseInt(seqCounter);
		// Return a reference to the global workspace.
		return workspace;
	}

	/**
	 * Obtain the directory where all the files pertaining to this workspace
	 * are stored by default.
	 * 
	 * @return The default directory where files for this workspace are stored.
	 */
	public String getDirectory() { return workspaceDirectory; }
	
	/**
	 * Obtain the ServerList that encapsulates the list of Servers that
	 * have been configured in this workspace.
	 * 
	 * @return The server list object that encapsulates all the servers configured
	 * on this workspace.
	 */
	public ServerList getServerList() { return serverList; }
	
	/**
	 * Obtain the JobList that encapsulates the list of Jobs that
	 * have been configured in this workspace.
	 * 
	 * @return The  list job object that encapsulates all the job run via
	 * this workspace.
	 */
	public JobList getJobList() { return jobList; }

	/**
	 * Obtain the list of data sets associated with this workspace.
	 * 
	 * @return The list of data set entries associated with this work space.
	 * The object returned by this method is never null. However, the 
	 * data set can be empty.
	 */
	public ArrayList<DataSet> getDataSets() { return dataSetList; }

	/**
	 * Find a data set entry given a ID.
	 * 
	 * This method can be used to find a data set with the given ID. This
	 * method iterates over the entries in {@link #dataSetList} and
	 * uses {@link DataSet#getID()} method to locate the desired entry. 
	 * 
	 * @param id The ID of the data set entry to be returned by this
	 * method. If this parameter is null, then this method immediately
	 * exits returning null.
	 * 
	 * @return This method returns the data set entry for the given ID.
	 * If the data set was not found then this method returns null.
	 */
	public DataSet findDataSet(final String id) {
		if (id == null) {
			// Nothing to search for.
			return null;
		}
		for(DataSet ds: dataSetList) {
			if (id.equals(ds.getID())) {
				return ds;
			}
		}
		return null;
	}
	
	/**
	 * Obtain the GeneratedFileList (GFL) entry for a given job ID.
	 * 
	 * @param jobID The ID of the job for which the generated/output
	 * files are to be searched and retrieved.
	 * 
	 * @return The GeneratedFileList corresponding to the given job 
	 * ID. If the entry was not found this method returns null.
	 */
	public GeneratedFileList getGFL(String jobID) {
		for(DataSet entry: dataSetList) {
			GeneratedFileList gfl = entry.getGFL(jobID);
			if (gfl != null) {
				return gfl;
			}
		}
		// Entry not found.
		return null;
	}
	
	/**
	 * Provides a string of the form Workspace (../path) to display
	 * some workspace related information to the user.
	 * 
	 * @return A simple string representation of this workspace.
	 */
	public String toString() {
		// Get the last entry in the path.
		return "Workspace";
	}
	
	/**
	 * Reserves the next sequence number for use in a data set. This method is
	 * used only when a data set object is added to this workspace.  Entries
	 * loaded from a persisted workspace file will have sequence numbers
	 * already filled in.
	 * 
	 * @return A new, unique (within the work space) ID for use with
	 * a new entry.
	 */
	public String reserveID() {
		// Use current sequence counter value to create ID
		String id = "x" + seqCounter;
		// Setup sequence counter for next server.
		seqCounter++;
		// Return the ID for use by the caller.
		return id;
	}

	/**
	 * Add a workspace listener to receive change notification events from
	 * this workspace. 
	 * 
	 * @param listener The listener to be added to the list of listeners to
	 * receive notifications.
	 */
	public synchronized void addWorkspaceListener(WorkspaceListener listener) {
		if (listener != null) {
			listeners.add(listener);
		}
	}

	/**
     * Remove a workspace listener from the listeners to receive notification 
     * events from this workspace. 
     * 
     * @param listener The listener to removed from the list of listeners to
     * receive notifications.
     */
    public synchronized void removeWorkspaceListener(WorkspaceListener listener) {
    	listeners.remove(listener);
    }
    
    /**
     * Utility method to notify all the listeners about a change to the Workspace.
     * 
     *  This method is used by various classes to report changes occurring to the
     *  workspace. This method is invoked by various classes constituting the
     *  workspace.
     *
     * @param event The event to be dispatched to various listeners informing them
     * about the change to the workspace.
     */
    public void fireWorkspaceChanged(WorkspaceEvent event) {
    	for(WorkspaceListener listener: listeners) {
    		listener.workspaceChanged(event);
    	}
    }
    
    /**
     * Method to add a new data set to the workspace.
     * 
     * This method adds the new data set to the workspace and notifies all the
     * listeners with a suitable event.
     * 
     * @param dataSet The data set to be added to the workspace.
     */
    public synchronized void addDataSet(DataSet dataSet) {
    	if (dataSet != null) {
    		dataSetList.add(dataSet);
    		// Notify all listeners about the change.
    		fireWorkspaceChanged(new WorkspaceEvent(dataSet, WorkspaceEvent.Operation.INSERT));
    	}
    }

    /**
     * Method to add a new data set to the workspace.
     * 
     * This method adds the new data set to the workspace and notifies all the
     * listeners with a suitable event.
     * 
     * @param dataSet The data set to be added to the workspace.
     */
    public synchronized void removeDataSet(DataSet dataSet) {
    	if (dataSet != null) {
    		dataSetList.remove(dataSet);
    		// Notify all listeners about the change.
    		fireWorkspaceChanged(new WorkspaceEvent(dataSet, WorkspaceEvent.Operation.DELETE));
    	}
    }
    
    /**
     * Obtain the classifier list set for this work space.
     * 
     * @return The classifier list set for this work space. 
     */
    public ClassifierList getClassifierList() {
    	return classifierList;
    }
    
	/**
	 * The default constructor. The constructor has been made private to ensure
	 * that this class is never directly instanted and there is only one unique
	 * instance of this class in each PEACE GUI process. The users of this class
	 * must use the get() method.
	 */
	private Workspace() {
		workspaceDirectory  = null;
		creationTimestamp   = null;
		dataSetList         = new ArrayList<DataSet>();
		serverList          = new ServerList();
		jobList             = new JobList();
		classifierList      = new ClassifierList();
	}

	/**
	 * The process-wide, unique instance of the workspace object. This object 
	 * contains all the necessary information regarding the workspace that
	 * is being currently used by the PEACE GUI. If a valid workspace is 
	 * currently unavailable, then this object will be null.
	 */
	private static Workspace workspace;

	/**
	 * The directory where all the files associated with this workpsace
	 * are stored. This value is set when a workspace is created or when
	 * a workspace is loaded for use.
	 */
	private String workspaceDirectory;

	/**
	 * The time stamp when this workspace was created. This is currently 
	 * unused and is present mostly for sanity checks to be done at a 
	 * later date.
	 */
	private XMLGregorianCalendar creationTimestamp;

	/**
	 * The list of data sets that have been configured for use in
	 * this work space. This list is either loaded from an XML file or
	 * populated by data set entries created by the user.
	 */
	private ArrayList<DataSet> dataSetList;
	
	/**
	 * The list of Server classes that have been configured for use in
	 * this work space. This list is either loaded from an XML file or
	 * populated by Server entries created by the user.
	 */
	private ServerList serverList;
	
	/**
	 * The list of Job classes that have been configured for use in
	 * this work space. This list is either loaded from an XML file or
	 * populated by Job entries created by the user.
	 */
	private JobList jobList;

	/**
	 * The list of data base classifiers that are used to classify and
	 * group ESTs obtained from different data bases. This information
	 * is persited and loaded from an XML work space file.
	 */
	private ClassifierList classifierList;
	
    /**
     * The list of listeners that are currently registered to receive 
     * notifications on changes that have occured to various entries
     * in the workspace. Listeners are managed via the addWorkspaceListener()
     * and removeWorkspaceListener() methods. The listeners are notified via
     * the fireWorkspaceChange() method.
     */
    private transient ArrayList<WorkspaceListener> listeners = 
    	new ArrayList<WorkspaceListener>();

	/**
	 * Sequence counter that is maintained on a per-work space basis to
	 * generate unique/valid IDs for each new entry that is added to this
	 * workspace.
	 */
	private long seqCounter;
}
