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

# This function is used to list jobs currently scheduled via
# PBS. This function takes one parameter, which is the user
# ID. If user ID is not specified then all jobs scheduled on
# the cluster are displayed.
function listPBSjobs {
	if [ -n "$1" ]; then
		# List jobs just for user where user ID is first parameter
		qstat -u $1
	else
		# List all jobs queued on this server.
		qstat
	fi
	# Return exit status of qstat back to caller
	return $?
}

# Function to simply list processes running on this machine.
function listJobs {
	if [ -n "$1" ]; then
		# List processes just for user where user ID is first parameter
		ps -f -u $1
	else
		# List all jobs queued on this server.
		ps -fea
	fi
	# Return exit status of ps back to caller
	return $?
}

#------------------------------------------------------------------------
#-------------------------------[ main ]---------------------------------
# Check and ensure we have exactly two command line arguments.
# The first one is the name of the script itself.
# The second argument is optional and must be the user-id

userID=""
if [ $# -gt 1 ]; then
	echo "Usage: listJobs.sh [userID]"
	exit 1
elif [ $# -eq 1 ]; then
	userID=$1
fi

# Check to see if PBS is up and running on this machine. The following
# function sets haveQsub variable to 1 if we have PBS.
checkPBS

if [ $haveQsub -eq 1 ]; then
	listPBSjobs $userID
else
	listJobs $userID
fi

# Exit with exit status of last run command
exit $?

# end of script
