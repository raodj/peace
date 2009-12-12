#!/bin/bash

#--------------------------------------------------------------------
#
# This file is part of PEACE.
# 
# PEACE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PEACE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
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

# Global flag to indicate if we have operational PBS system
haveQsub=0

# Check to see if we have operational PBS system by verifying
# qsub and qstat are working fine.
function checkPBS {
	which qsub > /dev/null 2> /dev/null
	if [ $? -eq 0 ]; then
	    # Check to ensure qstat etc. are working
		qstat -u $USER > /dev/null 2> /dev/null
		if [ $? -eq 0 ]; then
		    # Seems like PBS is up and running
			haveQsub=1
		fi
	fi
	# Everything went well
	return 0
}

# This function is used to create a PBS job file and submit the job
# via qsub.
function submitViaPBS {
	# First create a simple script for PBS.
	echo "#!/bin/bash"                               > job.sh
	echo "#PBS -N PEACE"                            >> job.sh
	echo "#PBS -l walltime=%maxRunTime%:00:00"      >> job.sh
	echo "#PBS -l nodes=%nodes%:ppn=%cpusPerNode%"  >> job.sh
	echo "#PBS -S /bin/bash"                        >> job.sh
	echo "#PBS -l mem=%memory%MB"                   >> job.sh
	echo "cd %workDir%"                             >> job.sh
	echo "time mpiexec %peace% %cmdLine%"           >> job.sh
	echo "echo $? > exit_status"                    >> job.sh

	# Ensure that the job.sh was created successfully.
	if [[ $? -ne 0 || !( -f job.sh ) ]]; then
		# Hmmm..some error occured
		return 1
	fi
	# Now finally run qsub and save job id
	qsub %workDir%/job.sh > job_id.qsub
	# Check to ensure there were no errors.
	if [[ $? -ne 0 || !( -f job_id.qsub ) ]]; then
		# Qsub job failed.
		return 2
	fi
	# Everything went well
	return 0
}

# Function to check status of a job submitted via PBS
function checkStatusViaPBS {
	# First check status of the job on PBS
	jobID=`cat job_id.qsub`
	qstat $jobID > /dev/null 2> /dev/null
	if [ $? -eq 0 ]; then
		echo "running"
		# Dummy exist status for consistency
		echo "-1"
	else
		echo "done"
		if [ -f exit_status ]; then
			cat exit_status
		else
			echo -1
		fi
	fi
	# Everything went well
	return 0
}

# Function used to directly run the job using mpiexec
# without using qsub.
function runCmdLine {
	# Compute number of CPUs requested.
	procs=$[%nodes% * %cpusPerNode%]
	# Start up a background job and redirect output streams
	outFile="job.o$$"
	errFile="job.e$$"
	scriptFile="job_$$.sh"
	# Run the job in the background and save exit status in
	# exit_status file. The following command line will be expanded
	# when this file is deployed and echoed to a script
	echo "#!/bin/bash" > $scriptFile
	echo "time mpiexec -n $procs %peace% %cmdLine% > $outFile 2> $errFile" >> $scriptFile
	echo 'echo $? > exit_status' >> $scriptFile
	chmod +x $scriptFile
	# Ensure that the runner script was created successfully.
	if [[ $? -ne 0 || !( -f $scriptFile ) ]]; then
		# Hmmm..some error occured
		return 1
	fi
	# Start up the job in background
	./$scriptFile < /dev/null 2>/dev/null 1>/dev/null &
	if [ $? -ne 0 ]; then
		echo "Unable to run job via command line" >&2
		return 2
	fi
	# Save PID of child process and create job_id file
	echo $! > job_id.pid
	# Everything went fine
	return 0
}

# Check to see if the non-PBS process is still running and if not
# print its exist status as well.
function checkStatus {
	pid=`cat job_id.pid`
	ps $pid > /dev/null 2> /dev/null
	if [ $? -eq 0 ]; then
		echo "running"
	    # Dummy exist status for output consistency
		echo "-1"
	else
		echo "done"
		if [ -f exit_status ]; then
			cat exit_status
		else
			echo -1
		fi
	fi
	# The task finished successfully
	return 0
}

# Function to print the contents of progress.dat in job directory
# if the file is present. Otherwise this function echo's -1,-1
function showProgress {
	# Check and print progress information if file is available.
	if [ -f progress.dat ]; then
		cat progress.dat
	else
		# Print invalid progress as we don't know progress.
		echo "-1,-1"
	fi
	# Everything went well
	return 0
}

#------------------------------------------------------------------------
#-------------------------------[ main ]---------------------------------
# Check and ensure we have exactly two command line arguments.
# The first one is the name of the script itself.
# The second argument must be "start", "status", or "output"
if [ $# -ne 1 ]; then
	echo "Usage: jobRunner.sh [start|status|output]"
	exit 1
fi

if [[ "$1" != "start" && "$1" != "status" && "$1" != "output" ]]; then
	echo "Usage: jobRunner.sh [start|status|output]"
	exit 1
fi

# Change to the job directory
cd %workDir%

if [ "$1" == "start" ]; then
    # Check to see if we have operational PBS system via a function
    # The following function sets the haveQsub flag.
	checkPBS
	
	if [ $haveQsub -eq 1 ]; then
	    # We have PBS. Submit job via PBS
		submitViaPBS
	else
	    # No pbs. Just run via command line
		runCmdLine
	fi
	# Return with exit status of last run command
	exit $?
elif [ "$1" == "status" ]; then
	# The user wants to check status of the job.
	# First print progress information.
	showProgress
	# Now do submission specific operation
	if [ -f job_id.qsub ]; then
		checkStatusViaPBS
	else
		checkStatus
	fi
	# Return with exit status of last run command
	exit $?
else
	# Echo the final results on standard out and standard error
	exit 3
fi

# end of script
