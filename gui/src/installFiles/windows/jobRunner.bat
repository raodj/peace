@echo off

REM --------------------------------------------------------------------
REM 
REM  This file is part of PEACE.
REM  
REM  PEACE is free software: you can redistribute it and/or modify it
REM  under the terms of the GNU General Public License as published by
REM  the Free Software Foundation, either version 3 of the License, or
REM  (at your option) any later version.
REM  
REM  PEACE is distributed in the hope that it will be useful, but
REM  WITHOUT ANY WARRANTY; without even the implied warranty of
REM  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
REM  General Public License for more details.
REM  
REM  You should have received a copy of the GNU General Public License
REM  along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
REM  
REM  Miami University makes no representations or warranties about the
REM  suitability of the software, either express or implied, including
REM  but not limited to the implied warranties of merchantability,
REM  fitness for a particular purpose, or non-infringement.  Miami
REM  University shall not be liable for any damages suffered by licensee
REM  as a result of using, result of using, modifying or distributing
REM  this software or its derivatives.
REM 
REM  By using or copying this Software, Licensee agrees to abide by the
REM  intellectual property laws, and all other applicable laws of the
REM  U.S., and the terms of GNU General Public License (version 3).
REM 
REM  Authors:   Dhananjai M. Rao              raodm@muohio.edu
REM 
REM ---------------------------------------------------------------------

:main
    REM Check and ensure we have exactly two command line arguments.
    REM The first one is the name of the script itself.
    REM The second argument must be "start", "status", "output", "scripts"

    if "%1" == "" goto showUsage

    REM Change to the job directory
    cd %workDir%
    REM change current working directory to workDir
    %workDrive%
    
    if "%1" == "start" (
        start jobRunner.bat runPeace
        goto end
    ) 
    if "%1" == "status" (
        call :status
	goto end
    ) 
    if "%1" == "output" (
        call :output
		goto end
	if "%1" == "scripts" (
        call :scripts
        goto end
    )
    if "%1" == "runPeace" (
        call :runPeace
        goto end
    )
:showUsage
    echo Usage: jobRunner.bat [start|status|output]
    exit 1

:end
    REM all went well
    exit 0

REM --- Sub routine to start PEACE with a given set of parameters ---
:runPeace
    echo Start time: %TIME% > job.stdout
    start /wait %peace% %cmdLine%  1>> job.stdout 2> job.stderr
    echo %ERRORLEVEL% > exit_status
    echo End time: %TIME%  >> job.stdout
    goto :EOF

:status
    REM First display progress information
    if exist progress.dat (
       type progress.dat
    ) else (
       echo -1,-1
    )

    REM check if job is done and print exit code
    findstr /C:"End time: " job.stdout 1> nul 2> nul
    if %ERRORLEVEL% EQU 0 (
       REM command is done running
       echo done
       type exit_status
    ) else (
       REM job is still running
       echo running
       echo -1
    )
    goto :EOF

:output
    REM dump output to standard output
    type job.stdout
    type job.stderr 1>&2
    goto :EOF

:scripts
    REM dump out script information
	REM Echo files used and the scripts to run the job.
	echo List of files in %CD%:
	dir
	echo ---------------------------------------------------------------
	
	REM Display progress information if we have any
	if exist progress.dat (
	   echo Contents of progress data file in %CD%\progress.dat:
       type progress.dat
       echo ---------------------------------------------------------------
	)

	REM Display exit status if file exists
	if exist exit_status (
	   echo Contents of progress data file in %CD%\exit_status:
       type exit_status
       echo ---------------------------------------------------------------
	)
	
    goto :EOF
	
REM end of file
