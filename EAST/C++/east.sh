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

# This script is a wrapper around EAST's Main executable. This
# script is used by the GUI. The primary feature that this script
# provides is that it automates the process of generating SAM
# file format from EAST's ACE output file format.

# The process of converting ACE-to-SAM is performed using
# AMOS tools (see http://sourceforge.net/apps/mediawiki/amos/)
# We assume AMOS tools are installed on the server on which
# EAST is run (via this script)

# We assume that EAST's Main executable is in the same path as
# this script. In addition, we assume that AMOS tools are
# available on the path.

# First detect if -VERSION command-line argument has been
# specified. If so, we just display VERSION information and exit
showVersion=`echo $* | grep -c -w "\-VERSION"`
if [ $showVersion -ne 0 ]; then
	echo "Expression fragment Assembly from Spanning Trees (EAST)"
	echo "EAST is a cDNA assembler."
	echo "Version 0.1. Released (under GPLv3) on December 12 2010"
	echo "Copyright (C) Miami University, Oxford, OH, USA (2009-)"
	exit 0
fi
		
# Ddetect if we need to do ACE-to-SAM conversion by searching
# for '-CONVERT_TO_SAM' command-line parameter.
convertToSAM=`echo $* | grep -c -w "\-CONVERT_TO_SAM"`

# Save the contig output file name for reference below.
# The contig file is the 3rd cmdline arg from the end
contig_pos=$(( $# - 2 ))
contigFile=`echo $* | cut -d' ' -f $contig_pos`

# First run east with the given command line parameters
# to generate ACE output format (if SAM conversion was selected)
# To do this we need to build full path to Main EAST executable
exePath=$0
exePathLen=$(( ${#exePath} - 7 ))
echo $exePathLen
cmdPath=${exePath:0:$exePathLen}
mainExe=$cmdPath"Main"
# Run EAST with full path
echo "$mainExe $*"
$mainExe $*

# Ensure EAST ran successfully. Otherwise bail out
eastExitCode=$?
if [ $eastExitCode -ne 0 ]; then
	# EAST run did not finish successfully.
	exit $eastExitCode
fi

# If SAM conversion is not needed then we are all done here!
if [ $convertToSAM -eq 0 ]; then
	# All done. Nothing further to be done.
	exit 0
fi

echo "Converting $contigFile to SAM file format..."
# Ensure that the output contig file is good.
if [ ! -f $contigFile -o ! -s $contigFile ]; then
	echo "EAST did not generate a valid contig file in $contigFile" >&2
	echo "Cannot convert to SAM file format" >&2
	exit 100
fi	

# Now copy the generated ACE file to a temporary file for further
# processing
srcFileName=${contigFile//\./\_}
aceFile=$srcFileName".ace"
mv $contigFile $aceFile
if [ $? -ne 0 ]; then
	echo "Unable to move $contigFile to $aceFile for ACE-to-SAM conversion" >&2
	exit 101
fi

# Next convert the generated ACE file to an AFG file
afgFile=$srcFileName".afg"
toAmos -ace $aceFile -o $afgFile
if [ $? -ne 0 ]; then
	echo "Unable to convert ACE file $aceFile to AFG file $afg" >&2
	exit 102
fi

# Create AMOS bank
bankFile=$srcFileName".bnk"
bank-transact -m $afgFile -b $bankFile -c
if [ $? -ne 0 ]; then
	echo "Unable to generate BANK file $bankFile from AFG file $afgFile" >&2
	exit 103
fi

# Create SAM output file.
bank2contig -i -s $bankFile > $contigFile
if [ $? -ne 0 ]; then
	echo "Unable to generate SAM file $contigFile from AMOS bank $bankFile" >&2
	exit 104
fi

# Successfully generated SAM file
exit 0

# End of script