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

echo "Starting site-specifing configure and compilation"
echo "Default install path is $PWD"

# Unzip the tar file to get all the source files
echo "Unzipping peace.tar.gz..."
tar zxf peace.tar.gz
if [ $? -ne 0 ]; then
   echo "Unzipping failed. Aborting install"
   exit 1
fi
# Change directory
cd peace
if [ $? -ne 0 ]; then
   echo "Directory change failed"
   exit 2
fi

echo "-----------------------------------------------------"
echo "Creating build system using autoconf..."
autoreconf -i -v 2>&1
if [ $? -ne 0 ]; then
   echo "Autoreconf failed to complete correctly"
   exit 3
fi

echo "-----------------------------------------------------"
echo "Configuring build system"
./configure
if [ $? -ne 0 ]; then
   echo "Configure failed."
   exit 4
fi

echo "-----------------------------------------------------"
echo "Building PEACE"
cores=`grep -c -i "Processor" /proc/cpuinfo`
make -j$cores 2>&1
if [ $? -ne 0 ]; then
   echo "Build Failed"
   exit 5
fi

echo "-----------------------------------------------------"
echo "Verifying operational status of PEACE"
# ./src/peace --options
if [ $? -ne 0 ]; then
   echo "PEACE executable is not operating as expected"
   exit 6
fi

echo "PEACE built successfully."

echo "-----------------------------------------------------"
echo "Creating data directoies."
mkdir estData
if [ $? -ne 0 ]; then
   echo "Shared EST data directory could not be created."
   exit 7
fi

echo "Creating directory for job entries"
mkdir jobs
if [ $? -ne 0 ]; then
   echo "Jobs directory could not be created."
   exit 8
fi

echo "-----------------------------------------------------"
echo "**      Installation completed successfully.       **"
echo "-----------------------------------------------------"

# end of script
