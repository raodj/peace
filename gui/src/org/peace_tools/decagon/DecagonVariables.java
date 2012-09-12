package org.peace_tools.decagon;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map.Entry;

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.core.Version;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.workspace.DataFileStats;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Param;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * A simple class to manage variables and values that can be 
 * used by DECAGON parameter definitions.
 * 
 * This class provides the various variable names and a hash table to hold
 * their values. This class contains variable definitions that are 
 * pertinent to the given context of operation.
 */
public class DecagonVariables {
	/**
	 * A DECAGON string variable that contains the PEACE GUI version. 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%PEACE_GUI_VERSION%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String PEACE_GUI_VERSION = "PEACE_GUI_VERSION";

	/**
	 * A DECAGON string variable that contains the full path to
	 * the working directory for a given job.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%WORK_DIRECTORY%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String WORK_DIRECTORY = "WORK_DIRECTORY";

	/**
	 * A DECAGON string variable that contains the name of the
	 * the target operating system on which the command is to be run. 
	 * Typically, this variable contains the values "Linux" or
	 * "Windows".
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%TARGET_OS%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String TARGET_OS = "TARGET_OS";

	/**
	 * A DECAGON string variable that contains the name of the
	 * the target server on which a job is to be run.  
	 * Typically, this variable contains the values such as:
	 * "redhawk.hpc.muohio.edu" or "glenn.osc.edu" etc.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%SERVER_NAME%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String SERVER_NAME = "SERVER_NAME";
	
	/**
	 * A DECAGON string variable that contains the path on a
	 * target server where PEACE/DECAGON has been installed.   
	 * Typically, this variable contains the values such as:
	 * "/home/raodm/PEACEtest"
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%INSTALL_PATH%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String INSTALL_PATH = "INSTALL_PATH";
	
	/**
	 * A DECAGON string variable that contains the user ID under
	 * which a server is being accessed. This user ID is the
	 * user ID specified by the user when a connection to a server
	 * is created. Typically, this variable contains the values 
	 * such as: "raodm". For local servers this value may not be 
	 * set.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%USER_ID%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String USER_ID = "USER_ID";

	/**
	 * A DECAGON string variable that contains the Job ID of the job
	 * that is currently being run by DECAGON. This job ID is 
	 * automatically created. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%JOB_ID%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String JOB_ID = "JOB_ID";
	
	/**
	 * A DECAGON string variable that contains the ID of the
	 * immediately previous job (if any). This value is useful
	 * if a DADX file uses multiple jobs and some files need
	 * to be passed from the previous job to the next job.
	 * This job ID is automatically created. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%PREV_JOB_ID%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */	
	public static final String PREV_JOB_ID = "PREV_JOB_ID";
	
	/**
	 * A DECAGON integer variable that contains the number of compute
	 * nodes (on a Supercomputing cluster) that have been reserved 
	 * for a job. Each DECAGON job can have a different value for 
	 * this variable depending on the user settings. Serial jobs 
	 * always have this value set to 1.  This variable contains only 
	 * integer values.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%COMPUTE_NODES%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String COMPUTE_NODES = "COMPUTE_NODES";

	/**
	 * A DECAGON integer variable that contains the number of CPUs
	 * per compute node reserved for a job. Each DECAGON job
	 * can have a different value for this variable depending on
	 * the user settings. Serial jobs always have this value set to
	 * 1.  This variable contains integer values.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%CPUS_PER_NODE%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String CPUS_PER_NODE = "CPUS_PER_NODE";
	
	/**
	 * A DECAGON integer variable that contains the maximum
	 * total memory (in Megabytes) reserved for the job on a 
	 * Supercomputing cluster that uses PBS for job scheduling.
	 * Each DECAGON job can have a different value for 
	 * this variable depending on the user settings. This 
	 * variable contains only integer values. This value is 
	 * not used if the Supercomputing cluster does not use
	 * PBS or other job scheduling infrastructure. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%MAX_MEMORY%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String MAX_MEMORY = "MAX_MEMORY";
	
	/**
	 * A DECAGON integer variable that contains the maximum
	 * run time (in hours) reserved for the job on a 
	 * Supercomputing cluster that uses PBS for job scheduling.
	 * Each DECAGON job can have a different value for 
	 * this variable depending on the user settings. This 
	 * variable contains only integer values. This value is 
	 * not used if the Supercomputing cluster does not use
	 * PBS or other job scheduling infrastructure. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%MAX_RUN_TIME%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String MAX_RUN_TIME = "MAX_RUN_TIME";

	/**
	 * A DECAGON string variable that contains the full path to
	 * the source genes corresponding to the cDNA reads in the 
	 * input file.  The relationship is that, once the cDNA reads
	 * in the input file have been assembled they should map to a
	 * source gene in the source gene file. This variable
	 * contains the path to the source gene file on the target
	 * server where the job is run.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%SRC_GENES_FILE_PATH%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String SRC_GENES_FILE_PATH = "SRC_GENES_FILE_PATH";

	/**
	 * A DECAGON string variable that contains the data storage
	 * format for the file referenced by {@link #SRC_GENES_FILE_PATH}.
	 * The value of this string is one of the following predefined
	 * strings: <code>FASTA</code>, <code>SFF</code>, <code>ACE</code>,
	 * <code>SAM</code>, <code>BAM</code>, <code>TSV</code>, or
	 * <code>TXT</code>.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%SRC_GENES_FILE_FORMAT%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String SRC_GENES_FILE_FORMAT = "SRC_GENES_FILE_FORMAT";

	/**
	 * A DECAGON string variable that contains simply the name of
	 * the file containing source cDNA fragments to be assembled.
	 * Unlike the {@link #INPUT_FILE_PATH} this variable does
	 * not contain any path information. This variable is useful
	 * for defining intermediate file names in DADX files.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%INPUT_FILE_NAME%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String INPUT_FILE_NAME = "INPUT_FILE_NAME";

	/**
	 * A DECAGON string variable that contains the full path to
	 * the input file from where data is to be read/processed. 
	 * Typically, this variable contains the path to the set
	 * of cDNA fragments to be processed/assembled.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%INPUT_FILE_PATH%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String INPUT_FILE_PATH = "INPUT_FILE_PATH";

	/**
	 * A DECAGON string variable that contains the data storage
	 * format for the file referenced by {@link #INPUT_FILE_PATH}.
	 * The value of this string is one of the following predefined
	 * strings: <code>FASTA</code>, <code>SFF</code>, <code>ACE</code>,
	 * <code>SAM</code>, <code>BAM</code>, <code>TSV</code>, or
	 * <code>TXT</code>.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%INPUT_FILE_FORMAT%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String INPUT_FILE_FORMAT = "INPUT_FILE_FORMAT";

	/**
	 * A DECAGON string variable that contains the full path to
	 * the quality file from where quality data for {@link #INPUT_FILE_PATH}
	 * is is to be read and used. Typically, this variable contains the 
	 * path to a FASTA quality style file for the set of cDNA fragments 
	 * to be processed/assembled.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%INPUT_QUALITY_FILE_PATH%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String INPUT_QUALITY_FILE_PATH = "INPUT_QUALITY_FILE_PATH";

	/**
	 * A DECAGON string variable that indicates if the data set 
	 * pointed by {@link #INPUT_FILE_PATH} contains Sanger-type 
	 * reads. This variable has the following three values:
	 * 
	 * <table border="0" cellSpacing="0">
	 * <tr><td>Undefined</td><td>This variable is not defined if the
	 * data set pointed by {@link #INPUT_FILE_PATH} is not a cDNA
	 * data set.</td></tr>
	 * 
	 * <tr><td><code>false</code></td><td>This variable is defined 
	 * and set to the string value <code>false</code> when the data 
	 * set pointed by {@link #INPUT_FILE_PATH} is a cDNA data set 
	 * but it does not contain any Sanger-type reads.</td></tr>
	 * 
	 * <tr><td><code>true</code></td><td>This variable is defined 
	 * and set to the string value <code>true</code> when the data 
	 * set pointed by {@link #INPUT_FILE_PATH} is a cDNA data set 
	 * and it contains Sanger-type reads.</td></tr> 
	 * </table>
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%HAS_SANGER_READS%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String HAS_SANGER_READS = "HAS_SANGER_READS";
	
	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * insertion-type error rate in Sanger sequences contained in the data set
	 * pointed by {@link #INPUT_FILE_PATH}, if it contains Sanger-type 
	 * reads. This variable is defined only when 
	 * {@link #HAS_SANGER_READS} has the value <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%SANGER_INS_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String SANGER_INS_ERR = "SANGER_INS_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * deletion-type error rate in Sanger sequences contained in the data set
	 * pointed by {@link #INPUT_FILE_PATH}, if it contains Sanger-type 
	 * reads. This variable is defined only when 
	 * {@link #HAS_SANGER_READS} has the value <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%SANGER_DEL_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String SANGER_DEL_ERR = "SANGER_DEL_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * substitution-type error rate in Sanger sequences contained
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains Sanger-type reads. This variable is defined 
	 * only when {@link #HAS_SANGER_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%SANGER_SUB_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String SANGER_SUB_ERR = "SANGER_SUB_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * total error rate in Sanger sequences contained
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains Sanger-type reads. This variable is defined 
	 * only when {@link #HAS_SANGER_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%SANGER_TOT_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String SANGER_TOT_ERR = "SANGER_TOT_ERR";

	/**
	 * A DECAGON string variable that indicates if the data set 
	 * pointed by {@link #INPUT_FILE_PATH} contains 454-type 
	 * reads. This variable has the following three values:
	 * 
	 * <table border="0" cellSpacing="0">
	 * <tr><td>Undefined</td><td>This variable is not defined if the
	 * data set pointed by {@link #INPUT_FILE_PATH} is not a cDNA
	 * data set.</td></tr>
	 * 
	 * <tr><td><code>false</code></td><td>This variable is defined 
	 * and set to the string value <code>false</code> when the data 
	 * set pointed by {@link #INPUT_FILE_PATH} is a cDNA data set 
	 * but it does not contain any 454-type reads.</td></tr>
	 * 
	 * <tr><td><code>true</code></td><td>This variable is defined 
	 * and set to the string value <code>true</code> when the data 
	 * set pointed by {@link #INPUT_FILE_PATH} is a cDNA data set 
	 * and it contains Roche 454-type reads.</td></tr> 
	 * </table>
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%HAS_R454_READS%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String HAS_R454_READS = "HAS_R454_READS";
	
	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * insertion-type error rate in 454 sequences contained in the data set
	 * pointed by {@link #INPUT_FILE_PATH}, if it contains Roche 454-type 
	 * reads. This variable is defined only when 
	 * {@link #HAS_R454_READS} has the value <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%R454_INS_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String R454_INS_ERR = "R454_INS_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * deletion-type error rate in Roche 454 sequences contained 
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains Roche 454-type reads. This variable is defined
	 * only when {@link #HAS_R454_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%R454_DEL_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String R454_DEL_ERR = "R454_DEL_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * substitution-type error rate in 454 sequences contained
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains Roche 454-type reads. This variable is defined 
	 * only when {@link #HAS_R454_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%R454_SUB_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String R454_SUB_ERR = "R454_SUB_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * total error rate in Roche 454 sequences contained
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains 454-type reads. This variable is defined 
	 * only when {@link #HAS_R454_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%R454_TOT_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String R454_TOT_ERR = "R454_TOT_ERR";

	/**
	 * A DECAGON string variable that indicates if the data set 
	 * pointed by {@link #INPUT_FILE_PATH} contains Illumina-type 
	 * reads. This variable has the following three values:
	 * 
	 * <table border="0" cellSpacing="0">
	 * <tr><td>Undefined</td><td>This variable is not defined if the
	 * data set pointed by {@link #INPUT_FILE_PATH} is not a cDNA
	 * data set.</td></tr>
	 * 
	 * <tr><td><code>false</code></td><td>This variable is defined 
	 * and set to the string value <code>false</code> when the data 
	 * set pointed by {@link #INPUT_FILE_PATH} is a cDNA data set 
	 * but it does not contain any Illumina-type reads.</td></tr>
	 * 
	 * <tr><td><code>true</code></td><td>This variable is defined 
	 * and set to the string value <code>true</code> when the data 
	 * set pointed by {@link #INPUT_FILE_PATH} is a cDNA data set 
	 * and it contains Illumina-type reads.</td></tr> 
	 * </table>
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%HAS_ILLUMINA_READS%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String HAS_ILLUMINA_READS = "HAS_ILLUMINA_READS";
	
	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * insertion-type error rate in Illumina sequences contained in the data set
	 * pointed by {@link #INPUT_FILE_PATH}, if it contains Illumina-type 
	 * reads. This variable is defined only when 
	 * {@link #HAS_ILLUMINA_READS} has the value <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%ILLUMINA_INS_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String ILLUMINA_INS_ERR = "ILLUMINA_INS_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * deletion-type error rate in Roche Illumina sequences contained 
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains Roche Illumina-type reads. This variable is defined
	 * only when {@link #HAS_ILLUMINA_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%ILLUMINA_DEL_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String ILLUMINA_DEL_ERR = "ILLUMINA_DEL_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * substitution-type error rate in Illumina sequences contained
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains Roche Illumina-type reads. This variable is defined 
	 * only when {@link #HAS_ILLUMINA_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%ILLUMINA_SUB_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String ILLUMINA_SUB_ERR = "ILLUMINA_SUB_ERR";

	/**
	 * A DECAGON numeric variable that indicates the estimated
	 * total error rate in Roche Illumina sequences contained
	 * in the data set pointed by {@link #INPUT_FILE_PATH}, if 
	 * it contains Illumina-type reads. This variable is defined 
	 * only when {@link #HAS_ILLUMINA_READS} has the value 
	 * <code>true</code>. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%ILLUMINA_TOT_ERR%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String ILLUMINA_TOT_ERR = "ILLUMINA_TOT_ERR";
	
	/**
	 * A DECAGON string variable that contains the full path to
	 * the working directory for another job that is dependent on
	 * this job. This information is used to appropriately start
	 * a dependent job (in a DADX pipeline) on the same server,
	 * once the current job has successfully completed. In addition,
	 * this information is used to copy the necessary dependent
	 * files to the appropriate folder prior to executing the job.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%WORK_DIRECTORY%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String DEP_WORK_DIRECTORY = "DEP_WORK_DIRECTORY";

	/**
	 * A DECAGON integer variable that contains the number of 
	 * processes associated with a given job. This value
	 * is set to the number of process entries specified in
	 * the DADX file for a given job.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%PROCESS_COUNT%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String PROCESS_COUNT = "PROCESS_COUNT";

	/**
	 * A DECAGON string variable that contains the full
	 * processed command lines for running each one of the processes
	 * associated with a given DECAGON job. Each of the command 
	 * lines is translated into a string by delimiting it with 
	 * double quotes. Any internal double quotes are suitably escaped.
	 * Each command line is placed on its own separate line. This style
	 * is useful for defining arrays in bash and other shells.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%PROCESS_CMD_LINE_LIST%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String PROCESS_CMD_LINE_LIST = "PROCESS_CMD_LINE_LIST";

	/**
	 * A DECAGON string variable that contains the name of the
	 * genome-assembler generated output contig file. This variable is
	 * set based based on the output file tagged as the contig file.
	 * This variable is used to process the generated output contig
	 * file appropriately. 
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%OUTPUT_CONTIG_FILE_NAME%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String OUTPUT_CONTIG_FILE_NAME = "OUTPUT_CONTIG_FILE_NAME";

	/**
	 * A DECAGON integer variable that contains the number of 
	 * output files generated by a given job. This value
	 * is set to the number of output entries specified in
	 * the DADX file for a given job.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%OUTPUT_FILE_COUNT%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String OUTPUT_FILE_COUNT = "OUTPUT_FILE_COUNT";

	/**
	 * A DECAGON string variable that contains the full
	 * path to output files to be generated by a given DECAGON job. 
	 * Each one of output file entries is translated into a string 
	 * by delimiting it with double quotes. Any internal double 
	 * quotes are suitably escaped. Each command line is placed 
	 * on its own separate line. This style is useful for 
	 * defining arrays in bash and other shells.
	 * 
	 * This variable can be used in parameters to define conditions
	 * or as meta variables (in the form 
	 * <code>%OUTPUT_FILE_LIST%</code>) in <code>cmdLine</code> or
	 * <code>value</code> XML elements in parameter definitions.  
	 */
	public static final String OUTPUT_FILE_LIST = "OUTPUT_FILE_LIST";

	/**
	 * Default constructor.
	 * 
	 * Default constructor to create a set of DecagonVariables. The
	 * constructor sets up the following default DECAGON variables
	 * and their values: {@link #PEACE_GUI_VERSION}
	 */
	public DecagonVariables() {
		varValues.put(PEACE_GUI_VERSION, Version.GUI_VERSION);
	}
	
	/**
	 * Setup various DECAGON variables corresponding to a given DataSet.
	 * 
	 * This method sets up the various DECAGON variables associated with
	 * sequencing technology and error rates using the data from 
	 * {@link DataSet#getStats()} methods. In addition, it sets the
	 * following DECAGON variables: {@link #INPUT_FILE_PATH},
	 * {@link #INPUT_FILE_FORMAT}.
	 * 
	 * @param ds The data set whose values are used to set up various
	 * DECAGON variables to be used in parameters.
	 */
	public void setVariables(final DataSet ds) {
		File srcFile = new File(ds.getPath());
		varValues.put(INPUT_FILE_PATH,   srcFile.getAbsolutePath());
		varValues.put(INPUT_FILE_FORMAT, ds.getFileType().toString());
		// Extract file name without path and extension information
		String fileName  = srcFile.getName();
		int dotPos       = fileName.lastIndexOf('.');
		dotPos           = (dotPos > 0 ? dotPos : fileName.length() - 1);
		fileName        = fileName.substring(0, dotPos);
		varValues.put(INPUT_FILE_NAME, fileName);
		
		// Setup the sequencing technology type and error variables.
		for (DataFileStats dfs: ds.getStats()) {
			switch (dfs.getTechType()) {
			case SANGER:
				setVariables(dfs, HAS_SANGER_READS, SANGER_DEL_ERR, SANGER_INS_ERR,
						SANGER_SUB_ERR, SANGER_TOT_ERR);
				break;
			case R454:
				setVariables(dfs, HAS_R454_READS, R454_DEL_ERR, R454_INS_ERR,
						R454_SUB_ERR, R454_TOT_ERR);
				break;
			case ILLUMINA:
				setVariables(dfs, HAS_ILLUMINA_READS, ILLUMINA_DEL_ERR, ILLUMINA_INS_ERR,
						ILLUMINA_SUB_ERR, ILLUMINA_TOT_ERR);
				break;
			default:
				ProgrammerLog.log("DecagonVariables.setVariables: " +
						"Unhandled technology type encountered: " + dfs.getTechType());
			}
		}
		// Update the path to point to the quality file on the server.
		// Update path to point to the source gene file on the server.
		final DataSet srcDS = Workspace.get().findDataSet(ds.getSourceID());
		if (srcDS != null) {
			File localSrcDSFile = new File(srcDS.getPath());
			varValues.put(SRC_GENES_FILE_PATH,   localSrcDSFile.getAbsolutePath());
			varValues.put(SRC_GENES_FILE_FORMAT, srcDS.getFileType().toString());
		}
	}
	
	/**
	 * Setup various DECAGON variables corresponding to a given Server.
	 * 
	 * This method sets up the various DECAGON variables associated with
	 * the operating system and other information associated with a given
	 * server. Specifically, this method sets up the following DECAGON
	 * variables:
	 * 
	 * <ul>
	 * 	<li> {@link #TARGET_OS} </li>
	 * 	<li> {@link #SERVER_NAME}</li>
	 *  <li> {@link #INSTALL_PATH}</li>
	 *  <li> {@link #USER_ID}</li>
	 * </ul>
	 * 
	 * @param ds The data set whose values are used to set up various
	 * DECAGON variables to be used in parameters.
	 */
	public void setVariables(final Server srvr) {
		varValues.put(TARGET_OS,    srvr.getOSType().name());
		varValues.put(SERVER_NAME,  srvr.getName());
		varValues.put(INSTALL_PATH, srvr.getInstallPath());
		if (srvr.getUserID() != null) {
			varValues.put(USER_ID,      srvr.getUserID());
		}
	}

	/**
	 * Setup various DECAGON variables corresponding to a full job
	 * configuration.
	 * 
	 * This method is a convenience method that can be used to 
	 * set various DECAGON variables associated with a job. The
	 * different constituents of a job are passed in as parameters.
	 * This method uses the following overloaded methods to setup
	 * some of the DECAGON variables: {@link #setVariables(DataSet)},
	 * {@link #setVariables(Server)}. In addition it sets up the
	 * following DECAGON variables associated with the job:
	 * 
	 * <ul>
	 * <li>{@link #COMPUTE_NODES}</li>
	 * <li>{@link #CPUS_PER_NODE}</li>
	 * <li>{@link #MAX_MEMORY}</li>
	 * <li>{@link #MAX_RUN_TIME}</li>
	 * <li>{@link #INPUT_FILE_NAME}</li>
	 * <li>{@link #INPUT_FILE_PATH}</li>
	 * </ul>
	 * 
	 * 
	 * @param ds The data set whose values are used to set up various
	 * DECAGON variables to be used in parameters. The variables
	 * are set by calling the {@link #setVariables(DataSet)} method.
	 * 
	 * @param srvr The server on which the job is to going to be run.
	 * The {@link #setVariables(Server)} method is used to set the
	 * variables associated with this server.
	 * 
	 * @param jobDir The full path on the server to the directory
	 * where the files for the job are to be stored.
	 * 
	 * @param jobID The PEACE job ID assigned to the current running job.
	 * 
	 * @param prevJobID The PEACE job ID assigned to the previous job 
	 * (if any). This value is used to set the {@link #PREV_JOB_ID} 
	 * DECAGON variable. If this parameter is NULL then this variable
	 * is cleared.
	 * 
	 * @param jobConfig This array has the elements in the order: 
	 * number of compute nodes, CPUs/node, total memory (in MB), 
	 * and runTime (in hours). This array is typically obtained via 
	 * call to 
	 * {@link org.peace_tools.core.job.ServerPanel#getPlatformConfiguration()} 
	 * method.
	 */
	public void setVariables(final DataSet ds, final Server srvr, 
			final String jobDir, final String jobID, 
			final String prevJobID, final int[] jobConfig,
			final String nextJobDir) {
		setVariables(ds);
		setVariables(srvr);
		varValues.put(JOB_ID,        jobID);
		varValues.put(COMPUTE_NODES, Integer.toString(jobConfig[0]));
		varValues.put(CPUS_PER_NODE, Integer.toString(jobConfig[1]));
		varValues.put(MAX_MEMORY,    Integer.toString(jobConfig[2]));
		varValues.put(MAX_RUN_TIME,  Integer.toString(jobConfig[3]));
		// Set/unset the value of PREV_JOB_ID
		if (prevJobID != null) {
			varValues.put(PREV_JOB_ID, prevJobID);
		} else {
			removeVariable(PREV_JOB_ID);
		}
		// Set/unset the value of DEP_WORK_DIRECTORY
		if (nextJobDir != null) {
			varValues.put(DEP_WORK_DIRECTORY, nextJobDir);
		} else {
			removeVariable(DEP_WORK_DIRECTORY);
		}
		// Update the INPUT_FILE_PATH to point to location on server.
		File localFile = new File(ds.getPath());
		varValues.put(INPUT_FILE_PATH, srvr.getServerPath(localFile));
		// Update the path to point to the quality file on the server.
		// Update path to point to the source gene file on the server.
		final DataSet srcDS = Workspace.get().findDataSet(ds.getSourceID());
		if (srcDS != null) {
			File localSrcDSFile = new File(srcDS.getPath());
			varValues.put(SRC_GENES_FILE_PATH,   srvr.getServerPath(localSrcDSFile));
			varValues.put(SRC_GENES_FILE_FORMAT, srcDS.getFileType().toString());
		}
		// Set the full working directory for the job
		varValues.put(WORK_DIRECTORY, jobDir);
	}
	
	/**
	 * Helper method to set DECAGON variables associated with a given 
	 * sequencing technology.
	 * 
	 * This method is invoked from the {@link #setVariables(DataSet)}
	 * method to set the error statistics variables associated with a given
	 * technology type. The variables are added to the {@link #varValues}
	 * hash table in this class.
	 * 
	 * @param dfs The data file statistics object from where the error
	 * statistics are to be obtained from the given sequencing technology.
	 * 
	 * @param techVarName The DECAGON variable name associated with the 
	 * technology. This variable is set to the string value <code>true</code>.
	 * 
	 * @param delVarName The DECAGON variable name associated with the
	 * deletion-type error rates for the given type of sequencing technology.
	 * 
	 * @param insVarName The DECAGON variable name associated with the
	 * insertion-type error rates for the given type of sequencing technology.
	 * 
	 * @param subVarName The DECAGON variable name associated with the
	 * substitution-type error rates for the given type of sequencing technology.
	 * 
	 * @param totErrVarName The DECAGON variable name associated with the
	 * total error rates for the given type of sequencing technology. This 
	 * value is computed by this method.
	 */
	private void setVariables(final DataFileStats dfs, final String techVarName,
			final String delVarName, final String insVarName, final String subVarName,
			final String totErrVarName) {
		varValues.put(techVarName,   "true");
		varValues.put(delVarName,    Float.toString(dfs.getDelErrorRate()));
		varValues.put(insVarName,    Float.toString(dfs.getInsErrorRate()));
		varValues.put(subVarName,    Float.toString(dfs.getSubErrorRate()));
		varValues.put(totErrVarName, Float.toString(dfs.getTotErrorRate()));
	}

	/**
	 * Add a variable to the list of variables.
	 * 
	 * This method can be used to add or update a variable and its value.
	 * If a variable does not exist then a new entry is created. If a
	 * variable with the same name already exists then its value is
	 * updated.
	 * 
	 * @param varName The name of the variable to be added to this list.
	 * This parameter cannot be null.
	 * 
	 * @param value The value associated with the variable. This parameter
	 * cannot be null.
	 */
	public void addVariable(final String varName, final String value) {
		varValues.put(varName, value);
	}
	
	/**
	 * Remove (un-define) a variable.
	 * 
	 * @param varName The name of the variable to be removed. If
	 * the variable does not exist then this method has no
	 * side effects.
	 */
	public void removeVariable(final String varName) {
		varValues.remove(varName);
	}
	
	/**
	 * Determine if a given variable is defined.
	 * 
	 * @param varName The name of the variable to be checked to verify if
	 * it is defined.
	 * 
	 * @return This method returns true if the variable is defined.
	 * Otherwise it returns false.
	 */
	public boolean isDefined(final String varName) {
		return varValues.containsKey(varName);
	}
	
	/**
	 * Get value for the variable as a string.
	 * 
	 * @param varName The name of the variable whose value is to be 
	 * returned as a string.
	 * 
	 * @return The string value of the variable if the variable is
	 * defined. Otherwise this method returns null.
	 */
	public String getStrValue(final String varName) {
		return varValues.get(varName);
	}

	/**
	 * Get value for the variable as a string.
	 * 
	 * @param varName The name of the variable whose value is to be 
	 * returned as a number.
	 * 
	 * @return The numeric value of the variable if the variable is
	 * defined. Otherwise this method returns null.
	 */
	public Double getNumValue(final String varName) {
		final String strVal = getStrValue(varName);
		Double retVal       = null;
		if (strVal != null) {
			try {
				retVal = Double.valueOf(strVal);
			} catch (NumberFormatException nfe) {
				// This exception is intentionally ignored.
			}
		}
		return retVal;
	}
	
	/**
	 * Convenience method to summary information about DECAGON variables.
	 * 
	 * This method can be used to summarize information about the 
	 * various DECAGON variables currently defined in this object.
	 * Currently, this method simply lists all variables and their
	 * values in a random order. Possibly this method can be revised
	 * to summarizes related variables in sections to make it more
	 * meaningful to the user.
	 * 
	 * @param sw The summary writer to be used to generate the summary
	 * of the variables.
	 * 
	 * @param subSummary If this flag is true then the summary information
	 * is generated using {@link SummaryWriter#addSubSummary(String, String, String)}
	 * method. If this flag is false, then the summary information is
	 * generated using {@link SummaryWriter#addSummary(String, String, String)}
	 * method.
	 */
	public void summarize(SummaryWriter sw, boolean subSummary) {
		for(Entry<String, String> entry: varValues.entrySet()) {
			if (subSummary) {
				sw.addSubSummary(entry.getKey(), entry.getValue(), "");
			} else {
				sw.addSummary(entry.getKey(), entry.getValue(), "");
			}
		}
	}
	
	/**
	 * Convenience method to copy all the variable -> value pairs
	 * in this object to a workspace parameter list.
	 * 
	 * @param paramList The target workspace parameter list to which
	 * all the DECAGON variables and values are to be copied.
	 */
	public void copyAll(ArrayList<Param> paramList) {
		for(Entry<String, String> entry: varValues.entrySet()) {
			Param wsEntry = new Param(entry.getKey(), entry.getValue());
			paramList.add(wsEntry);
		}
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(1024);
		for(Entry<String, String> entry: varValues.entrySet()) {
			sb.append(entry.getKey());
			sb.append(": ");
			sb.append(entry.getValue());
			sb.append("\n");
		}
		return sb.toString();
	}
	
	/**
	 * Hash table used to hold the currently defined set of variables
	 * and their values.
	 */
	private final Hashtable<String, String> varValues = 
		new Hashtable<String, String>();
}
