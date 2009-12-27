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
        start %launcher% %launcher% %peace% %cmdLine%
        goto end
    ) 
    if "%1" == "status" (
        call :status
	goto end
    ) 
    if "%1" == "output" (
        call :output
	goto end
    )
    if "%1" == "error" (
        call :error
	goto end
    )
    if "%1" == "scripts" (
        call :scripts
        goto end
    )
    if "%1" == "abort" (
        call :abort
        goto end
    )
:showUsage
    echo Usage: jobRunner.bat [start|status|output|error|scripts|abort]
    exit 1

:end
    REM all went well
    exit 0

:status
    REM First display progress information
    if exist progress.dat (
       type progress.dat
    ) else (
       echo -1,-1
    )

    REM check if job is done and print exit code
    if exist exit_status (
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
    REM dump stdout to standard output
    type job.stdout
    goto :EOF
    
:error
	REM dump stderr to standard output
    type job.stderr
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

:abort
    if not exist job_id.pid (
	echo Unable to determine PID of PEACE process 1>&2
	exit 10
    )
    REM Read PID into variable.
    set PID=-1
    FOR /F %%A IN ('type job_id.pid') DO set PID=%%A

    REM validate PID
    IF %PID% EQU -1 (
	echo Unable to determine PID of PEACE process 1>&2
	exit 11
    )

    REM Ensure this PID is sane and valid and there is a peace 
    REM process running under this PID
    tasklist /FI "UserName eq %USERNAME%" /FI "PID eq %PID%" /FI "ImageName eq peace.exe" | FINDSTR "peace.exe"
    IF %ERRORLEVEL% NEQ 0 (
	echo The process for this job is no longer running and the 1>&2
	echo job cannot be aborted. 1>&2
	exit 12
    )
    REM OK, kill the PEACE process..
    TASKKILL /PID %PID% 
    IF %ERRORLEVEL% NEQ 0 (
	echo Unable to kill process with PID %PID% 1>&2
        echo Unable to abort PEACE job 1>&2
	exit 13
    )
    goto :EOF


REM end of file
