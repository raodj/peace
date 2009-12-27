//--------------------------------------------------------------------
//
// This file is part of PEACE.
// 
// PEACE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// PEACE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
// 
// Miami University makes no representations or warranties about the
// suitability of the software, either express or implied, including
// but not limited to the implied warranties of merchantability,
// fitness for a particular purpose, or non-infringement.  Miami
// University shall not be liable for any damages suffered by licensee
// as a result of using, result of using, modifying or distributing
// this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

/**
    This is a simple program that has been developed to launch peace
    from a batch file without holding up the console. The output from
    PEACE are redirected to job.stdout and job.stderr files.
*/

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>

// C Runtime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>

// C++ Runtime headers
#include <string>
#include <fstream>

#include "Resource.h"

#ifdef _UNICODE
#define tstring wstring
#else
#define tstring string
#endif

// Forward declarations of functions included in this code module:
HANDLE createOutputFile(const TCHAR* fileName);
HANDLE createProcess(const TCHAR *, const TCHAR*, const TCHAR *, HANDLE, HANDLE, DWORD&);

int APIENTRY _tWinMain(HINSTANCE hInstance,
                     HINSTANCE hPrevInstance,
                     LPTSTR    lpCmdLine,
		     int       nCmdShow)  {
	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(nCmdShow);
	UNREFERENCED_PARAMETER(hInstance);

	DWORD exitCode = 3;
	DWORD pid      = 0;

	if ((lpCmdLine == NULL) || (_tcslen(lpCmdLine) == 0)) {
	    // Missing command line arguments.
	    return 1;
	}

	// Extract the PEACE executable path and parameters from the first
	// Search for the word peace.exe in the command line parameter
	std::tstring cmdLine(lpCmdLine);
	int index = cmdLine.find(_T("peace.exe"));
	std::tstring exePath = cmdLine.substr(0, index + 9);
	std::tstring peaceParams = cmdLine.substr(index + 9);

	// Create files for redirecting stdout and stderr.
	HANDLE stdOut = createOutputFile(_T("job.stdout"));
	HANDLE stdErr = createOutputFile(_T("job.stderr"));
	if ((stdOut == NULL) || (stdErr == NULL)) {
	    // Could not create the necessary output files.
	    // Don't proceed further
	    return 2;
	}
	// Start up PEACE process with necessary parameters and
	// redirect stdout and stderr to files.
	HANDLE peaceProcess = createProcess(exePath.c_str(), peaceParams.c_str(), 
	    NULL, stdOut, stdErr, pid);
	if (peaceProcess != NULL) {
	    // Write the PID of the process to file.
	    std::wofstream pidFile(_T("job_id.pid"));
	    pidFile << pid << std::endl; 
	    pidFile.close();
	    // Wait until child process exits.
	    WaitForSingleObject(peaceProcess, INFINITE );
	    // Get the exit code to ensure things went fine.
	    GetExitCodeProcess(peaceProcess, &exitCode);
            // Close the process handle as we no longer need it.
	    CloseHandle(peaceProcess);
	}
	// Dump exit code out to appropriate file.
	std::wofstream statusFile(_T("exit_status"));
	statusFile << exitCode << std::endl; 
	statusFile.close();

	// Close output streams.
	CloseHandle(stdOut);
	CloseHandle(stdErr);
	// Exit with exit code of the process
	return (int) exitCode;
}

HANDLE
createOutputFile(const TCHAR* fileName) {
    SECURITY_ATTRIBUTES saAttr;
    saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
    saAttr.bInheritHandle = TRUE;
    saAttr.lpSecurityDescriptor = NULL;

    // Setup the handles for PEACE to use.
    return CreateFile(fileName, GENERIC_WRITE, 0, 
	&saAttr, CREATE_ALWAYS,	FILE_ATTRIBUTE_NORMAL, NULL);
}

HANDLE
createProcess(const TCHAR *appName, const TCHAR* arguments,
              const TCHAR *directory, HANDLE stdOut, HANDLE stdErr,
	      DWORD& pid) {
    PROCESS_INFORMATION processInfo;
    STARTUPINFO         startupInfo;
    
    ZeroMemory(&startupInfo, sizeof(STARTUPINFO));
    startupInfo.cb          = sizeof(startupInfo);
    startupInfo.wShowWindow = SW_SHOWMINIMIZED;
    startupInfo.dwFlags     = STARTF_USESHOWWINDOW;

    if ((stdOut != NULL) || (stdErr != NULL)) {
        startupInfo.dwFlags    |= STARTF_USESTDHANDLES;
	startupInfo.hStdOutput  = stdOut;
	startupInfo.hStdError   = stdErr;
	startupInfo.hStdInput   = INVALID_HANDLE_VALUE;
    }

    // Create process and note result.
    TCHAR *argCopy = _tcsdup(arguments);
    bool retVal = CreateProcess(appName, argCopy, NULL, NULL, TRUE,  
                                NORMAL_PRIORITY_CLASS, 
                                NULL,       // Environment
                                directory,  // CWD
                                &startupInfo,
                                &processInfo) == TRUE;
    // Extract child process handle if successful.
    HANDLE childProc = NULL;
    if (retVal) {
        // Sucessfully created process.  Clean up unused handles.
        childProc = processInfo.hProcess;
	pid       = processInfo.dwProcessId;
        CloseHandle(processInfo.hThread);
    }
    // Get rid of the memory for the duplicate arguments.
    free(argCopy);

    // Return handle to child process.
    return childProc;
}