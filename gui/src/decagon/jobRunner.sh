#!/bin/bash

#--------------------------------------------------------------------
#
# This file is part of DECAGON and PEACE.
# 
# PEACE and DECAGON are free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# PEACE and DECAGON are distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PEACE and DECAGON. If not, see
# <http://www.gnu.org/licenses/>.
# 
# Miami University makes no representations or warranties about the
# suitability of the software, either express or implied, including
# but not limited to the implied warranties of merchantability,
# fitness for a particular purpose, or non-infringement.  Miami
# University shall not be liable for any damages suffered by licensee
# as a result of using, result of using, modifying or distributing
# this software or its derivatives.
#
# By using or copying this Software, Licensee agrees to abide by the
# intellectual property laws, and all other applicable laws of the
# U.S., and the terms of GNU General Public License (version 3).
#
# Authors:   Dhananjai M. Rao              raodm@muohio.edu
#
#---------------------------------------------------------------------

# This script is processed by DECAGON GUI before it is finally copied
# to the server on which the job is to be run. The processing of this
# script essentially involves substituting DECAGON variables with
# their appropriate values.  A more detailed description of the
# DECAGON variables is presented further below.  This script provides
# just two methods that are used by the various helper functions in
# the jobRunnerIncl.sh which is included at the end of this
# script. Refer to the jobRunnerIncl.sh script for details on the
# various operations supported by the job-runner script(s).
#
# This script is processed by the DECAGON sub-system in the following
# manner.  Strings/variables in the form %variable% (identifier
# surrounded by percentage sign) are meta variables whose values are
# substituted by DECAGON GUI when it creates a job. They may be used
# anywhere in the script and GUI simply does a string-search-replace
# operation on these variables. The following variables are made
# available by DECAGON for developing this script:
#
# PEACE_GUI_VERSION  Version of PEACE (string)
# WORK_DIRECTORY     Work directory for the current job (string) 
# TARGET_OS          The operating system on which this job is running (string)
# SERVER_NAME        Name or IP address of server (string)
# INSTALL_PATH       Absolute install path for PEACE on server (string)
# USER_ID            The user ID under which the job is running (string)
# JOB_ID             The unique ID assigned to the job (string)
# PREV_JOB_ID        The unique ID assigned to the previous job (string)
# COMPUTE_NODES      Number of compute nodes reserved for the job (int)
# CPUS_PER_NODE      Number of CPUs to request per node for job (int)
# MAX_MEMORY         The maximum memory in MB (int > 256)
# MAX_RUN_TIME       The maximum runtime in hours (int)

# INPUT_FILE_PATH    The path to the cDNA file to be processed (string)
# INPUT_FILE_NAME    Just the name (without path and suffix) of the input
#                    cDNA file to be processed (string)
# INPUT_FILE_FORMAT  File format for the source gene file (string)
# DEP_WORK_DIRECTORY Directory of follow-up job to be run (string)
# HAS_SANGER_READS   Defined if input file has sanger-type reads
# PROCESS_COUNT      Number of processes in this job (int)
#
# PROCESS_CMD_LINE_LIST   A new-line separated list of command-lines
# INPUT_QUALITY_FILE_PATH Path to the input quality file (string)
# SRC_GENES_FILE_PATH     Path to source genes (if any) (string)
# SRC_GENES_FILE_FORMAT   File format for the source gene file (string)

# The list of processes to be run for this job
procList=(%PROCESS_CMD_LINE_LIST%)

# The time command prefix and suffix used to measure the wall-clock
# runtime of the program.
TIME_PRE=' { time -p '
TIME_SUF=' ; } 2>>timings.txt'

# This function is used to create a PBS job file
function createPBSJobFile {
	# Create a simple script for PBS. This script uses DECAGON
	# variables that are subsituted by the DECAGON-GUI.
	echo "#!/bin/bash"                                        > job.sh
	echo "#PBS -N %JOB_ID%"                                  >> job.sh
	echo "#PBS -l walltime=%MAX_RUN_TIME%:00:00"             >> job.sh
	echo "#PBS -l nodes=%COMPUTE_NODES%:ppn=%CPUS_PER_NODE%" >> job.sh
	echo "#PBS -S /bin/bash"                                 >> job.sh
	echo "#PBS -l mem=%MAX_MEMORY%MB"                        >> job.sh
	echo "cd %WORK_DIRECTORY%"                               >> job.sh

	# A comment line to demarcate job outputs
	echo 'echo "# Timings for job %JOB_ID%" >> timings.txt'  >> job.sh

	# Generate command-lines for each one of the processes. The
	# function is defined in jobRunnerIncl.sh
	generateProcessCalls "$TIME_PRE mpiexec" "$TIME_SUF"
	# Generate the list of files to be copied to next job. This
	# array is used by genRunNextJob call below (which also copies
	# timings.txt to the next directory)
	echo "outFileList=(%OUTPUT_FILE_LIST% 'timings.txt')"    >> job.sh
	# Generate qsub command to schedule a follow up job (if any)
	genRunNextJob "%DEP_WORK_DIRECTORY%"
	# Everything went well
	return 0
}

# Function used to create a job.sh file that directly runs the job
# using mpiexec without using qsub.
function createCmdLineJobFile {
	# Compute number of CPUs requested.
	procs=$[%COMPUTE_NODES% * %CPUS_PER_NODE%]
	# Setup targets for redirecting output streams
	outFile="job.o$$"
	errFile="job.e$$"
	scriptFile="job.sh"
	# Check and set if we have mpicc to run in parallel via mpi
	# and change the procPrefix below
	procPrefix="$TIME_PRE"
	which mpicc 1>/dev/null 2>/dev/null
	if [ $? -eq 0 ]; then
		procPrefix="$TIME_PRE mpiexec -n $procs"
	fi
	
	# Run the job in the background and save exit status in
	# exit_status file. The following command line will be expanded
	# when this file is deployed and echoed to a script
	echo "#!/bin/bash" >  job.sh
	echo ""            >> job.sh

	# A comment line to demarcate job outputs
	echo 'echo "# Timings for job %JOB_ID%" >> timings.txt'  >> job.sh	
	
	# Generate command-lines for each one of the processes. The
	# function is defined in jobRunnerIncl.sh
	generateProcessCalls "$procPrefix" ">> $outFile 2>> $errFile $TIME_SUF"
	# Generate the list of files to be copied to next job
	echo "outFileList=(%OUTPUT_FILE_LIST%"  >> job.sh
	echo "             'timings.txt')"        >> job.sh	
	# Generate qsub command to schedule a follow up job (if any)
	genRunNextJob "%DEP_WORK_DIRECTORY%"
	# Everything went well
	return 0	
}


# The main part of the script that performs various tasks depending on
# the command-line arguments (supplied by GUI)
source %INSTALL_PATH%/decagon/scripts/jobRunnerIncl.sh

# Change to the job directory
cd %WORK_DIRECTORY%

# Call the main function in jobRunnerIncl.sh with all command-line args
main $*

# end of script
