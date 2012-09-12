package org.peace_tools.decagon.helpers;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.xml.XMLConstants;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.SchemaOutputResolver;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.Result;
import javax.xml.transform.dom.DOMResult;
import javax.xml.transform.dom.DOMSource;
import javax.xml.validation.SchemaFactory;

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.core.job.ServerPanel;
import org.peace_tools.decagon.DecagonVariables;
import org.peace_tools.decagon.jaxb.Argument;
import org.peace_tools.decagon.jaxb.AssemblerDetails;
import org.peace_tools.decagon.jaxb.Executable;
import org.peace_tools.decagon.jaxb.Job;
import org.peace_tools.decagon.jaxb.OutputCheckData;
import org.peace_tools.decagon.jaxb.OutputFile;
import org.peace_tools.decagon.jaxb.OutputFileType;
import org.peace_tools.decagon.jaxb.Parameter;
import org.peace_tools.decagon.jaxb.Process;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Server.EXEKind;

/**
 * DECAGON Assembler Description XML (DADX) File Helper.
 * 
 * This is a top-level convenience class to help read parameters from an 
 * DECAGON Assembler Description XML (DADX) file and generate GUI 
 * components for/from a given DADX file. This class aids the remainder 
 * of PEACE and DECAGON GUI components with the following operations:
 * 
 * <ul>
 * 
 * <li>Unmarshal (load/read) an assembler description from a given XML
 * file. The XML file must be compatible with the Assembler.xsd schema.</li>
 * 
 * <li>Facilitate creation of menu components for the DECAGON menu in 
 * PEACE GUI.</li>
 * 
 * <li>Interface with the ParameterSetGUIHelper to facilitate creating GUI
 * components for the various parameters specified for the assembler.</li>
 * 
 * <li>Expose various features of the super class ParameterSetHelper
 * (such as: manage DECAGON variables, handle conditional statements,
 * and bind variable substitution).</li>
 * 
 * </ul>
 */
public class DADXHelper extends ParameterSetGUIHelper {
	/**
	 * The DECAGON assembler description that is currently being 
	 * handled by this parameter set.
	 */
	protected AssemblerDetails assemblerInfo = null;

	/**
	 * The list of DADX file from where this DADX helper has loaded
	 * data so far. Entries are added to this file by the
	 * {@link #unmarshal(String, boolean, boolean)} method.
	 */
	protected ArrayList<String> sourceDADXFilePaths =
		new ArrayList<String>();
	
	/**
	 * The default and only constructor for this class.
	 * 
	 * The constructor is relatively straightforward and merely initializes
	 * the various instance variables to their default initial value.
	 */
	public DADXHelper() {
		// Nothing else to be done here for now.
	}
	
	/**
	 * Obtain the top-level object that encapsulates assembler information.
	 * 
	 * @return The assembler details object encapsulated by this helper.
	 * If valid assembler details have not yet been loaded, then this
	 * method returns null.
	 */
	public AssemblerDetails getAssemblerDetails() {
		return assemblerInfo;
	}
	/**
	 * Method to create schemas to be used for validating a DADX file.
	 * 
	 * This method was introduced because JAXB does not seems to seamlessly
	 * handle using external schema files for validating DADX files.  Instead
	 * this method extracts the internal schema (created by JAXB using 
	 * annotations on the various classes generated from the DECAGON
	 * schema files) via JAXB API calls and returns that schema to be
	 * used for validation.
	 * 
	 * @param context The JAXBContext that has been created using the 
	 * top-level class (namely AssemblerDetails) for parsing DADX files.  
	 * 
	 * @return An array containing the source for the internal schema that
	 * can be used for validating DADX files.
	 * 
	 * @throws IOException This method merely exposes an exception that
	 * could be generated. Typically this exception is not generated because
	 * this method does not write any data to a file.
	 */
	private DOMSource[] getJaxbSchemas(JAXBContext context) throws IOException {
	    final List<DOMResult> results = new ArrayList<DOMResult>();
	    context.generateSchema(new SchemaOutputResolver() {
	        @Override
	        public Result createOutput(String ns, String file) throws IOException {
	            DOMResult result = new DOMResult();
	            result.setSystemId(file);
	            results.add(result);
	            return result;
	        }
	     });
		DOMSource[] schemaSources  = new DOMSource[results.size()];
		for(int i = 0; (i < results.size()); i++) {
			schemaSources[i] = new DOMSource(results.get(i).getNode());
		}
	    return schemaSources;
	}

	/**
	 * Overrides base class method to load parameters from an DADX file.
	 * 
	 * <p>This method provides a convenient mechanism to load an assembler
	 * description from a XML file. The XML file must be compatible
	 * with the Assembler.xsd schema file.</p>
	 * 
	 * <p>NOTE: Unmarshalling the a reasonable sized DADX file takes
	 * a few seconds. In order to ensure good user experience, it is
	 * preferable to invoke this method from a background worker thread.</p>
	 * 
	 * @param dadxFileName The name of the DADX file from where
	 * the description for an assembler is to be loaded. The
	 * file is opened via call to {@link Utilities#getStream(String)}
	 * method.
	 * 
	 * @param register If this flag is true, then a successfully loaded
	 * DADX description is added to the list of known genomic-assemblers
	 * via {@link DADXSummary#register(AssemblerDetails)}.
	 * 
	 * @param merge If this flag is true and an earlier DADX file
	 * was successfully loaded, then this method merges the jobs and
	 * parameters from the newly loaded DADX file into the existing one. 
	 * 
	 * @throws Exception This method throws different types of exceptions
	 * depending of the type of error that occurs when loading the
	 * parameters from the XML description file.
	 * 
	 * @see Utilities#getStream(String)
	 */
	public void unmarshal(final String dadxFileName,  
			final boolean register, boolean merge) throws Exception {
		// Obtain URL to the file to be loaded 
		final URL fileURL = resolveFilePath(dadxFileName, null);
		if (fileURL == null) {
			// Unable to resolve file name. Throw an exception
			throw new IllegalArgumentException("Unable to locate DADX file " + dadxFileName);
		}
		// create JAXB context and instantiate unmarshaller
		JAXBContext context = JAXBContext.newInstance(AssemblerDetails.class);
		Unmarshaller um = context.createUnmarshaller();
		SchemaFactory sf = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
		um.setSchema(sf.newSchema(getJaxbSchemas(context)));
		// Setup handler to intercept validation errors.
		JAXBValidationHandler pvh = new JAXBValidationHandler();
		um.setEventHandler(pvh);
		// Get the unmarshaller to load the data for us.
		final InputStream xml = Utilities.getStream(dadxFileName);
		// Process the files depending on whether it needs to be merged
		if ((assemblerInfo != null) && merge) {
			AssemblerDetails asmInfo = (AssemblerDetails) um.unmarshal(xml);
			merge(asmInfo);
		} else {
			assemblerInfo = (AssemblerDetails) um.unmarshal(xml);
			paramSet      = assemblerInfo.getAssemblyParams();
		}
		// Check and throw exception if any validation errors occur
		final String validationErrors = pvh.getErrorLog();
		if (validationErrors.length() != 0) {
			// Validation errors occurred
			throw new JAXBException(validationErrors);
		}
		// Check to ensure that there is at least one output file with 
		// the isContigFile flag set to true.
		if (this.getOutputContigFile() == null) {
			throw new JAXBException("The DADX file did not contain an OutputFile entry " +
					"tagged as the final generated contig file (via isContigFile attribute " + 
					"being set to true");
		}
		// DADX file loaded successfully.
		this.sourceDADXFilePaths.add(dadxFileName);
		if (register) {
			DADXSummary.register(fileURL.getPath(), assemblerInfo);
		}
	}

	/**
	 * Method to merge information from another AssemblerDetails
	 * into the currently loaded one.
	 * 
	 * This method merges the necessary information from another 
	 * assembler details object into the currently loaded one. 
	 * The following information is merged from the other 
	 * description into the current one:
	 * 
	 * <ul>
	 * <li>The set of parameters are merged.</li>
	 * <li>The set of executables are merged.</li>
	 * <li>The set of jobs are merged.</li>
	 * </ul>
	 *  
	 * <p><b>NOTE</b>: This method does not modify the identifiers and
	 * names of various entities during merging. Consequently, conflicts
	 * due to names being identical could arise.</p>
	 * 
	 * @param other The other assembler details to be merged.
	 * This object cannot be null and must be a valid assembler
	 * information object.
	 */
	public void merge(AssemblerDetails other) {
		// First merge all the parameters that the other job may 
		// be dependent on.
		final List<Parameter> currParamList  = assemblerInfo.getAssemblyParams().getParameter();
		final List<Parameter> otherParamList = other.getAssemblyParams().getParameter();
		for(Parameter param: otherParamList) {
			currParamList.add(param);
		}
		// Next merge the executables referred to in the jobs.
		final List<Executable> currExeList  = assemblerInfo.getExecutable();
		final List<Executable> otherExeList = other.getExecutable();
		for(Executable exe: otherExeList) {
			currExeList.add(exe);
		}
		// Finally merge the jobs from the other details with our current ones
		final List<Job> currJobList  = assemblerInfo.getJob();
		final List<Job> otherJobList = other.getJob();
		for(Job otherJob: otherJobList) {
			currJobList.add(otherJob);
		}
	}
	
	/**
	 * Obtain the path to the DADX file from where data for this helper
	 * has been loaded thus far.
	 * 
	 * This method essentially returns {@link #sourceDADXFilePaths}
	 * instance variable. Entries are added to this list by the
	 * {@link #unmarshal(String, boolean, boolean)} method once
	 * a valid DADX file is successfully processed and loaded.
	 * 
	 * @return This method returns the list of DADX files from
	 * where data has been successfully loaded into this class
	 * so far. This method always returns a valid array list
	 * but the list can be empty. 
	 */
	public ArrayList<String> getSourceDADXFilePaths() {
		return this.sourceDADXFilePaths;
	}
	
	/**
	 * Helper method to get command-line string from a given list of
	 * arguments.
	 * 
	 * This method is a helper method that can be used to obtain the 
	 * the command-line string to be used to run a given executable.
	 * This method processes the list of arguments supplied to this
	 * method and performs the following operations on each argument
	 * entry:
	 * 
	 * <ol>
	 * <li>It first checks any conditions specified on the argument 
	 * by calling {@link #checkCondition(org.peace_tools.decagon.jaxb.Condition)}
	 * method. If the condition evaluates to <code>false</code> then
	 * this method ignores the argument and moves onto the next entry.</li>
	 * 
	 * <li>If the argument is a parameter reference, then the 
	 * corresponding parameter is obtained (via call to 
	 * {@link #getParameter(String)} method) and a command-line entry is
	 * added for the corresponding parameter (via call to
	 * {@link #addToCmdLine(Parameter, List)} helper method).</li>
	 * 
	 * <li>If the argument directly provides a static string to be used
	 * then this method adds it to the final command-line string after
	 * processing any DECAGON variable references in it (via call
	 * to {@link #processVariables(String)} method).
	 * </ol>
	 * 
	 * @param argList The list of arguments to be processed and converted
	 * to a command-line string. This parameter cannot be null.
	 * 
	 * @return The command-line string to be used to run the executable
	 * using given list of arguments. The return value is never null
	 * but can be an empty string.
	 */
	public String getCmdLineArgs(final List<Argument> argList) {
		// The outgoing list of command-line entries
		ArrayList<String> cmdLine = new ArrayList<String>();
		// Process each argument in the arg list.
		for(Argument arg: argList) {
			// If argument has a condition, then check if condition returns true
			if ((arg.getCondition() != null) && 
				!checkCondition(arg.getCondition())) {
				// The condition indicates this argument is not to be used
				continue;
			}
			// Process either a static command-line string or a parameter
			// reference depending on what is specified in the DADX file.
			if (arg.getParameterRef() != null) {
				// Find the actual parameter entry for this parameter reference.
				final Parameter param = getParameter(arg.getParameterRef());
				if (param != null) {
					// Found the parameter. Use helper method to add necessary
					// entries to the cmdLine list.
					addToCmdLine(param, cmdLine);
				}
			} else if (arg.getCmdLineArg() != null) {
				// This is a string to be used as command-line argument.
				// Substitute DECAGON variable references using helper method.
				final String cmdLineEntry = processVariables(arg.getCmdLineArg());
				cmdLine.add(cmdLineEntry);
			}
		}
		// Convert the array list into a single string. Entries are
		// separated by one white space.
		return Utilities.toString(cmdLine, " ");
	}
	
	/**
	 * Convenience method to generate summary information for a given job.
	 * 
	 * This is a helper method that can be used to generate summary information
	 * (to be displayed to the user) for a given job.
	 *  
	 * @param sw The summary writer object to be used for generating the
	 * summary information. This object cannot be null.
	 * 
	 * @param job The job object for which the summary information is to be
	 * generated. This object cannot be null.
	 * 
	 * @param jobConfig The job configuration object that provides additional
	 * job configuration information. This object cannot be null.
	 */
	public void summarize(SummaryWriter sw, Job job, final ServerPanel jobConfig) {
		sw.addSection(job.getName() + " Job Summary");
		sw.addSummary("Description", null, job.getDescription());
		sw.addSummary("Job Type", job.getType().toString(), "");
		if (jobConfig != null) {
			jobConfig.summarize(sw, true);
		}
		// Summarize information about each process associated with the job(s)
		for(Process proc: job.getProcess()) {
			summarize(sw, proc, jobConfig.getSelectedServer());
		}
		// Summarize information about the output files generated by this job.
		sw.addSubSection("Output files", "", "List of files generated by the job");
		for(OutputFile of: job.getOutputFile()) {
			final String actualFilePath = processVariables(of.getPath());
			sw.addSubSummary(of.getFileType().toString(), actualFilePath, 
					of.getDescription());
		}
	}

	/**
	 * Obtain the Executable object with a given name.
	 * 
	 * This method iterates over the list of executables associated with
	 * the current assembler pipeline (via call to {@link AssemblerDetails#getExecutable()}
	 * method) and returns the executable whose name (obtained via
	 * {@link Executable#getName()} method) is the same as the given parameter.
	 * 
	 * @param exeName The name of the executable to search for. This value
	 * cannot be null.
	 * 
	 * @return The executable object whose name is the same as the given exeName.
	 * If a matching executable object is not found then this method return null.
	 */
	public Executable getExecutable(final String exeName) {
		for(Executable exe: assemblerInfo.getExecutable()) {
			if (exeName.equals(exe.getName())) {
				return exe;
			}
		}
		return null;
	}
	
	/**
	 * Obtain the actual executable name for a given executable on a given 
	 * server.
	 * 
	 * This is a convenience method that can be used to consistently obtain
	 * the final executable name for both internal and external executables.
	 * 
	 * @param exe The executable object whose final executable name is to
	 * be returned by this method. This object cannot be null.
	 * 
	 * @param srvr The server on which the executable is going to be run.
	 * The server information is used to decide on the actual executable
	 * name for internal executables.
	 * 
	 * @return The final executable name to be used for running the 
	 * actual process on the given server.
	 */
	public String getExeName(Executable exe, Server srvr) {
		if (exe.getFileName() != null) {
			// This is an external executable
			return exe.getFileName();
		}
		// Handle internal executable
		switch(exe.getInternal()) {
		case EAST:
			return srvr.getExecutable(EXEKind.EAST_EXE);
		case PEACE_CLUSTERING:
		case PEACE_ASSEMBLY:
			return srvr.getExecutable(EXEKind.PEACE_CLUSTER_MAKER);
		case DECAGON_ANALYZER:
			return srvr.getExecutable(EXEKind.DECAGON_ANALYZER);
		}
		return "";
	}
	
	/**
	 * Convenience method to generate the full command-line for a
	 * DECAGON process on a given server.
	 * 
	 * This is a helper method that can be used to generate the full
	 * command-line for executing a DECAGON process on a given server.
	 * This method handles both internal and external executables. 
	 * It handles various command-line arguments (specified for the
	 * process) while adhering to any conditions and handling any variables
	 * via call to {@link #getCmdLineArgs(List)} method.
	 * 
	 * @param proc The process for which this method is to create the 
	 * command-line for execution. This object cannot be null.
	 * 
	 * @param srvr The server on which the process is going to be executed.
	 * This object is used to determine the actual, platform-specific executable
	 * name for internal executables. This object cannot be null.
	 */	
	public String getProcessCmdLine(Process proc, Server srvr) {
		// Get the actual executable name for this process
		final Executable exe     = getExecutable(proc.getExecutable());
		final String     exeName = this.getExeName(exe, srvr);
		// Obtain the various command-line arguments while processing conditions
		// and handling any variables via call to parent-class method.
		final String cmdLine = exeName + " " + getCmdLineArgs(proc.getCmdLineArgs().getArgument());
		// Return the full command-line back to the caller
		return cmdLine;
	}
	
	/**
	 * Convenience method to setup DECAGON variable to refer to the
	 * generated output contig file.
	 * 
	 * This method can be used to set the DECAGON variable that 
	 * refers to the generated output contig file. This method
	 * essentially sets the value for the DECAGON variable
	 * {@link org.peace_tools.decagon.DecagonVariables#OUTPUT_CONTIG_FILE_NAME}.
	 * This method uses the {@link #getOutputContigFile()} method to
	 * locate the generated output file. 
	 * 
	 * @return This method returns true if the variable was successfully set.
	 * Otherwise (the {@link #getOutputContigFile()} returns false) this method
	 * returns false.
	 */
	public boolean setOutputContigVariable() {
		final OutputFile contigFile = getOutputContigFile();
		if (contigFile != null) {
			decVar.addVariable(DecagonVariables.OUTPUT_CONTIG_FILE_NAME, 
					processVariables(contigFile.getPath()));
		}
		return (contigFile != null);
	}
	
	/**
	 * Convenience method to search the set of output files in various
	 * jobs to locate an output contig file.
	 * 
	 * Each DADX description of a genomic-assembler software pipeline must
	 * contain at least one output file that contains the contigs
	 * generated by the genomic-assembler from the given set of input
	 * cDNA fragments. This output file is identified using the
	 * <code>isContigFile</code> attribute set for an output file entry
	 * in the DADX file. This method searches the output files for the
	 * jobs in the DADX file to locate the appropriate entry and return
	 * it.
	 *  
	 * @return The output file object associated with the generated
	 * contig file. If the output file is not found, then this method
	 * returns null.
	 */
	public OutputFile getOutputContigFile() {
		OutputFile contigFile = null;
		for(Job job: assemblerInfo.getJob()) {
			for(OutputFile of: job.getOutputFile()) {
				if (of.isIsContigFile()) {
					contigFile = of;
					break;
				}
			}
		}
		return contigFile;
	}
	
	/**
	 * Helper method to convert file mime-type from DECAGON to PEACE 
	 * mime-type.
	 * 
	 * This is a convenience method that can be used to convert DECAGON
	 * file mime-types defined in OutputFileType to the set of PEACE
	 * mime types defined in DataSet.DataFileType.
	 *  
	 * @param oft The DECAGON mime-type to be converted to PEACE mime-type.
	 * 
	 * @return The PEACE mime-type corresponding to the given DECAGON
	 * mime-type.
	 */
	public static DataSet.DataFileType toDataFileType(OutputFileType oft) {
		DataFileType retVal = DataFileType.valueOf(oft.name());
		if (retVal == null) {
			retVal = DataFileType.TXT;
		}
		return retVal;
	}
	
	/**
	 * Convenience method to generate summary information for a given process.
	 * 
	 * This is a helper method that can be used to generate summary information
	 * (to be displayed to the user) for a given process.
	 * 
	 * @param sw The summary writer object to be used for generating the
	 * summary information. This object cannot be null.
	 * 
	 * @param proc The process for which this method is to generate
	 * summary information. This object cannot be null.
	 * 
	 * @param srvr The server on which the process is going to be executed.
	 * This object is used to determine the actual, platform-specific executable
	 * name for internal executables. This object cannot be null.
	 */
	public void summarize(SummaryWriter sw, Process proc, Server srvr) {
		sw.addSubSection("Process Information", "", 
				"Information about a process in a job");
		// Get the actual executable name for this process
		final Executable exe     = getExecutable(proc.getExecutable());
		final String     exeName = this.getExeName(exe, srvr);
		final String cmdLine = exeName + " " + getCmdLineArgs(proc.getCmdLineArgs().getArgument());
		sw.addSubSummary("Command line", cmdLine, "Command line to be run");
		sw.addSubSummary("Expected exit code:", proc.getExitCode().toString(),
				"Expected exit code indicating successful processing");
		if (proc.getOutputCheck() != null) {
			final OutputCheckData ocd = proc.getOutputCheck();
			sw.addSubSummary("Stream(s) to check: ", (ocd.isStderr() ?
					" stderr" : "") + (ocd.isStdout() ? " stdout " : ""),
					"The output(s) from the process to checked");
			sw.addSubSummary("Additional output checks: ", ocd.getValue(),
					"Regular expression to search for in the output");
		}
	}
	/**
	 * Obtain the list of DADX files to be loaded when PEACE starts up.
	 * 
	 * This is an interface method between PEACE and DECAGON. It is invoked
	 * when a workspace is launched from {@link org.peace_tools.core.WorkspaceChooser}.
	 * This method searches for additional DADX files specified by the
	 * default PEACE directory (typically ~PEACE) and adds them to list
	 * of standard DADX files to be loaded. This information is used
	 * by the {@link org.peace_tools.core.WorkspaceChooser} to pre-load
	 * the necessary files (in a background worker thread).
	 * 
	 * @return A list of DADX files to be pre-loaded when a workspace is
	 * launched.
	 */
	public static List<String> getDADXFilesToLoad() {
		// Set standard set of DECAGON Assembler Description XML (DADX) files 
		// to be loaded.
		final ArrayList<String> dadxList = new ArrayList<String>(10);
		dadxList.add("decagon/EAST.dadx");
		// Determine the set of any additional
		// DECAGON Assembler Description XML (DADX) files to be loaded.
		final String defDirPath = Utilities.getDefaultDirectory();
		File defDir = new File(defDirPath);
		String[] additionalDADXFiles = defDir.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith("dadx");
			}			
		});
		if ((additionalDADXFiles != null) && (additionalDADXFiles.length > 0)) {
			// Append additional entries to the dadxList.
			for(int i = 0; (i < additionalDADXFiles.length); i++) {
				String addFilePath = defDirPath + File.separator + additionalDADXFiles[i];
				dadxList.add(addFilePath);
			}
		}
		return dadxList;
	}
	
	/**
	 * Convenience method to resolve the path to a given file.
	 * 
	 * This method is a helper method that was introduced to streamline the 
	 * code in the {@link #unmarshal(String, boolean)} method. This method
	 * does best effort to locate a given file in the following
	 * manner:
	 * 
	 * <ol>
	 * <li>First it checks to see if the given file is a resource
	 * that can be obtained directly from PEACE jar. If so, it returns
	 * a suitable URI to it.</li>
	 * 
	 * <li>If the file is not a resource in the JAR it checks to see if
	 * it is available in the default folders distributed with
	 * PEACE.</li>
	 * 
	 * <li>It checks to see if the path is a reference to a
	 * standard file on the machine and if so returns a URL to it.</li>
	 * 
	 * <li>Lastly, if a parent path is supplied (that is the second
	 * parameter to this method is not null) and all previous options
	 * did not yield success, then this method prefixes the parent path 
	 * to the given file path to determine if the file is found and 
	 * if so returns a URL to it.</li>
	 * 
	 * </ol>
	 * 
	 * @param path The path to a DADX file to which a canonical URL is
	 * desired.
	 * 
	 * @param parentPath An optional parent directory path in which to
	 * search for the given file. This parent path can be null and is 
	 * used as the last option to locate a given file.
	 * 
	 * @return A canonical URL to the given DADX file, if the path is
	 * valid. If the path is invalid then this method returns null.
	 * 
	 * @throws MalformedURLException This method exposes a possible
	 * exception that maybe generated by the methods used by this class.
	 * Generation of this exception is a rare case.
	 */
	protected static URL resolveFilePath(final String path, 
			final String parentPath) throws MalformedURLException {
		// Try to see if this a file in our JAR
		URL uri = Utilities.class.getResource("/" + path);
		// Try to see if it is available from other path in our jar (used when
		// running on desktop directory without a Jar).
		if (uri == null) {
			uri = Utilities.class.getResource(Utilities.PATH_PREFIX + path);
		}
		// Next simply try the path
		if (uri == null) {
			File tempFile = new File(path);
			if (tempFile.exists() && tempFile.canRead()) {
				uri = tempFile.toURI().toURL();
			}
		}
		// Next try the path with the parent path prefixed to it.
		if ((uri == null) && (parentPath != null)) {
			File tempFile = new File(parentPath + File.separator + path);
			if (tempFile.exists() && tempFile.canRead()) {
				uri = tempFile.toURI().toURL();
			}
		}
		return uri;
	}
}
