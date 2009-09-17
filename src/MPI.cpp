#ifndef MPI_CPP
#define MPI_CPP

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include "MPI.h"
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

#endif

#endif
