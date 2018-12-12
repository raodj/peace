#ifndef MPI_HELPER_CPP
#define MPI_HELPER_CPP

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

#include "MPIHelper.h"
#include <cstring>

#ifndef HAVE_LIBMPI

#ifndef _WINDOWS
// A simple implementation for MPI_WTIME on linux
#include <sys/time.h>
double MPI_WTIME() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (tv.tv_usec / 1e6);
}

#else
// A simple implementation for MPI_WTIME on Windows
#include <windows.h>

double MPI_WTIME() {
    FILETIME st;
    GetSystemTimeAsFileTime(&st);
    long long time = st.dwHighDateTime;
    time <<= 32;
    time |= st.dwLowDateTime;
    return (double) time;
}


#endif

#else

// Convenience method to return the MPI rank when MPI library is
// available
int MPI_GET_RANK() {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

// Convenience method to return the number of processes when MPI
// library is available
int MPI_GET_SIZE() {
    int size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

#endif

int getOwnerProcess(const int estListSize, const int estIdx) {
    const int eSTsPerProcess = estListSize / MPI_GET_SIZE();
    const int extraESTs      = estListSize % MPI_GET_SIZE();
    const int extraESTsCutOff= (extraESTs * eSTsPerProcess) + extraESTs;
    
    // If the estIdx is less that the ExtraESTsCutOff then account for
    // the fact that each of these workers have one extra ESTs (as the
    // number of ESTs may not be evenly divisible by number of
    // processes).
    if (estIdx < extraESTsCutOff) {
        return estIdx / (eSTsPerProcess + 1);
    }
    // This est is in a process that does not have one extra...
    return (estIdx - extraESTs) / eSTsPerProcess;
}

void getLocallyOwnedESTidx(const int count, int& startIndex, int& endIndex) {
    const int eSTsPerProcess = count / MPI_GET_SIZE();
    const int extraESTs      = count % MPI_GET_SIZE();
    const int myRank         = MPI_GET_RANK();
    
    // First figure out the starting and ending EST this processs is
    // responsible for further use.
    startIndex = myRank * eSTsPerProcess;
    // Check for extra proceses as needed.
    if (myRank <= extraESTs) {
        // The previous processes have one extra ESTs as number of
        // ESTs are not evenly divisible by number of processes.  So
        // account for this.
        startIndex = ((eSTsPerProcess + 1) * myRank);
    } else {
        startIndex += extraESTs;
    }
    
    // Compute the last est index this process owns.
    endIndex = startIndex + eSTsPerProcess;
    if (myRank < extraESTs) {
        // This process owns one extra EST as ESTs are not evenly
        // divisible by the number of processes.
        endIndex++;
    }
}

#ifdef HAVE_LIBMPI
int MPI_GET_COUNT(MPI_STATUS& status, MPI_Datatype mpiType) {
    int count = 0;
    MPI_Get_count(&status, mpiType, &count);
    return count;
}
#else
int MPI_GET_COUNT(MPI_STATUS&, MPI_Datatype) {
    return 0;
}
#endif

#endif
