#ifndef PROGRESS_REPORTER_CPP
#define PROGRESS_REPORTER_CPP

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

#include "ProgressReporter.h"
#include <iostream>
#include <cstdlib>

// Create the static, singleton instance
ProgressReporter ProgressReporter::uniquePR;

ProgressReporter::~ProgressReporter() {
    if (progFile.is_open()) {
        progFile.close();
    }
}

void
ProgressReporter::initialize(const std::string& outFile, const int maxSteps) {
    // Save max steps for future reference.
    this->maxSteps = maxSteps;
    // Try and create the output file. If that fails, abort
    progFile.open(outFile.c_str());
    if (!progFile.good()) {
        std::cerr << "Unable to open progress file " << progFile
                  << " for writing. Aborting." << std::endl;
        abort();
    }
    // Report initial progress
    reportProgress(0);
}

void
ProgressReporter::reportProgress(const int progress) {
    if (progFile.good()) {
        progFile << progress << "," << maxSteps << "\n"
                 << std::flush;
        progFile.seekp(0);
    }
}

#endif
