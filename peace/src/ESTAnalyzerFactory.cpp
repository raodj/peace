#ifndef EST_ANALYZER_FACTORY_CPP
#define EST_ANALYZER_FACTORY_CPP

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

#include "ESTAnalyzerFactory.h"
#include "arg_parser.h"

#include "FMWSCA.h"
#include "CLU.h"
#include "D2.h"
#include "D2Zim.h"
#include "AdaptiveTwoPassD2.h"
#include "TwoPassD2.h"
#include "OldTwoPassD2.h"
#include "MatrixFileAnalyzer.h"

void
ESTAnalyzerFactory::displayList(std::ostream &os) {
    // Create dummy command-line args to make display prettier and
    // easier.
    arg_parser::arg_record dummy_args[] = {
        {"fmwsca", "Framed, Multi-Word String Compare Analyzer",
         NULL, arg_parser::STRING},
        {"clu", "CLU's similarity metric generation algorithm",
         NULL, arg_parser::STRING},
        {"matrixFile", "Use distance data stored in a matrix file",
         NULL, arg_parser::STRING},
        {"d2", "Use D2 (wcd-based) distance metric generation algorithm",
         NULL, arg_parser::STRING},
        {"d2zim", "Use D2 (Zimmerman) distance metric generation algorithm",
         NULL, arg_parser::STRING},
        {"twopassD2", "Use two-pass asymmetric/bounded symmetric D2",
         NULL, arg_parser::STRING},
        {"twopassD2adapt", "Use two-pass asymmetric/bounded symmetric D2, adaptive mode",
         NULL, arg_parser::STRING},
        //        {"d2sim", "Use special version of D2 for simulation project",
        //         NULL, arg_parser::STRING},
        {NULL, NULL, NULL, arg_parser::BOOLEAN}
    };
    arg_parser dummyParser(dummy_args);
    os << dummyParser;
}

ESTAnalyzer*
ESTAnalyzerFactory::create(const char* name, const int refESTidx,
                           const std::string& outputFileName) {
    if (name == NULL) {
        return NULL;
    }
    if (refESTidx < 0) {
        std::cerr << "A reference EST index has not been specified.\n"
                  << "Use --estIdx command line option." << std::endl;
        return NULL;
    }
    
    if (!strcmp("fmwsca", name)) {
        return new FMWSCA(refESTidx, outputFileName);
    } else if (!strcmp("clu", name)) {
        return new CLU(refESTidx, outputFileName);
    } else if (!strcmp("matrixFile", name)) {
        return new MatrixFileAnalyzer(refESTidx, outputFileName);
    } else if (!strcmp("d2", name)) {
        return new D2(refESTidx, outputFileName);
    } else if (!strcmp("d2zim", name)) {
        return new D2Zim(refESTidx, outputFileName);        
    } else if (!strcmp("twopassD2", name)) {
        return new TwoPassD2(refESTidx, outputFileName);
    } else if (!strcmp("oldTwopassD2", name)) {
        return new OldTwoPassD2(refESTidx, outputFileName);
    } else if (!strcmp("twopassD2adapt", name)) {
        return new AdaptiveTwoPassD2(refESTidx, outputFileName);
        //    } else if (!strcmp("d2sim", name)) {
        //        return new SimulationD2(refESTidx, outputFileName);
    }
    
    // invalid analyzer name!
    std::cerr << "Invalid analyzer name." << std::endl;
    return NULL;
}

#endif
