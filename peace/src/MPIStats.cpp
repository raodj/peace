#ifndef MPI_STATS_CPP
#define MPI_STATS_CPP

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

#include "MPIStats.h"

// Define and initialize various static stats. counters in the
// MPIStats class for use.

double MPIStats::idleTime   = 0.0;
int    MPIStats::sendCount  = 0;
int    MPIStats::recvCount  = 0;
int    MPIStats::bcastCount = 0;
int    MPIStats::probeCount = 0;

void
MPIStats::displayStats(std::ostream& os, const std::string& indent) {
    os << indent << "No. of calls to send     : " << sendCount  << "\n"
       << indent << "No. of calls to recv     : " << recvCount  << "\n"
       << indent << "No. of calls to bCast    : " << bcastCount << "\n"
       << indent << "No. of calls to probe    : " << probeCount << "\n"
       << indent << "Net amount of idle time  : " << idleTime
       << " (sec)" << std::endl;
}

#endif
