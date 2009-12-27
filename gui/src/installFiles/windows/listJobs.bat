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
    REM If we have a parameter then assume it is user ID and list only
    REM jobs for that user.

    if "%1" == "" (
    	tasklist
    ) else (
	tasklist /FI "UserName eq %USERNAME%"
    )

REM End of batch file