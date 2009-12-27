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

echo "Starting site-specific configure and compilation"
echo "Default install path is $PWD"

# Unzip the tar file to get all the source files
echo "Unzipping peace.tar.gz..."
tar zxf peace.tar.gz
if [ $? -ne 0 ]; then
   echo "---------------------[ Error ]-------------------------"
   echo "Unzipping failed. This machine is not configured correctly"
   echo "or your environment is broken. This is not a PEACE problem."
   echo "Aborting install"
   exit 1
fi

# Change directory
cd peace
if [ $? -ne 0 ]; then
   echo "---------------------[ Error ]-------------------------"
   echo "Directory change failed!"
   exit 2
fi

echo "-----------------------------------------------------"
echo "Creating build system using autoconf..."
autoreconf -i -v 2>&1
if [ $? -ne 0 ]; then
   echo "-----------------------[ Error ]-----------------------------"
   echo "Autoreconf failed to complete correctly. This indicates that" 
   echo "your version of autoconf is too old. You need at least version"
   echo "1.9 of autoconf. You need to upgrade autoconf or contact your"
   echo "system administrator. This is not an issue with PEACE."
   echo "Aborting install"
   exit 3
fi

# Detect and workaround absence of mpicc
useMpi=""
which mpicc
if [ $? -ne 0 ]; then
	echo "Warning: mpicc not found in path."
	echo "Warning: Proceeding with default C++ compiler"
	useMpi="--without-mpi"
fi

echo "-----------------------------------------------------"
echo "Configuring build system"
./configure $useMpi
if [ $? -ne 0 ]; then
   echo "-----------------------[ Error ]-----------------------------"
   echo "Configure failed. This indicates that some of the software"    
   echo "required by PEACE could not be found on this machine. Refer"   
   echo "to the error logs above to install the necessary tools. If"    
   echo "you are on a cluster environment check with your administrator"  
   echo "to load suitable modules by default so the necessary software" 
   echo "is available by default. This is not a PEACE issue."           
   echo "Aborting install"
   exit 4
fi

echo "-----------------------------------------------------"
echo "Building PEACE"
cores=1
# If the machine has cpu info and machine has many cores do a
# parallel build to speedup installation.
if [ -f /proc/cpuinfo ]; then
	cores=`grep -c -i "^Processor" /proc/cpuinfo`
	# Don't exceed 4-way build. Not polite on shared clusters
	if [ $cores -gt 4 ]; then
		cores=4
	fi
fi

make -j$cores 2>&1
if [ $? -ne 0 ]; then
   echo "-----------------------[ Error ]-----------------------------"
   echo "Build Failed. This indicates that there is some unforeseen" 
   echo "incompatibility between PEACE and the software installed  " 
   echo "on this machine. If you would like to have this issue     " 
   echo "resolved, please email the complete output to developers. " 
   # Simply dump config.log to make life easier.
   echo "------------ contents of config.log -------------------"
   cat config.log
   echo "-------------------------------------------------------"
   echo "Aborting install"
   exit 5
fi

echo "-----------------------------------------------------"
echo "Verifying operational status of PEACE"
# ./src/peace --options
if [ $? -ne 0 ]; then
   echo "-----------------------[ Error ]-----------------------------"
   echo "PEACE executable is not operating as expected. This is some"
   echo "unforeseen situation. Most likely this is not a PEACE issue"
   echo "Aborting install"
   exit 6
fi

echo "PEACE built successfully."

echo "-----------------------------------------------------"
echo "Creating data directories."
mkdir estData
if [ $? -ne 0 ]; then
   echo "-----------------------[ Error ]-----------------------------"
   echo "Shared EST data directory could not be created. This is some" 
   echo "unforeseen situation. Most likely this is not a PEACE issue"  
   echo "Aborting install"
   exit 7
fi

echo "Creating directory for job entries"
mkdir jobs
if [ $? -ne 0 ]; then
   echo "-----------------------[ Error ]-----------------------------"
   echo "Jobs directory could not be created. This is some unforeseen" 
   echo "situation. This is not  PEACE issue."                         
   echo "Aborting install"
   exit 8
fi

echo "-----------------------------------------------------"
echo "**      Installation completed successfully.       **"
echo "-----------------------------------------------------"

# end of script
