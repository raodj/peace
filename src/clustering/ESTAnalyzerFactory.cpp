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
#include "ArgParser.h"
#include "config.h"

#include "FMWSCA.h"
#include "CLU.h"
#include "D2.h"
#include "D2Zim.h"
#include "AdaptiveTwoPassD2.h"
#include "TwoPassD2.h"
#include "OldTwoPassD2.h"
#include "MatrixFileAnalyzer.h"
//#include "BatonAnalyzer.h"
#include "PrimesESTAnalyzer.h"

#ifdef HAVE_CUDA
#include "D2Cuda.h"
#endif

#include "AlignmentAnalyzer.h"
#include "TMAlignAnalyzer.h"

void
ESTAnalyzerFactory::addCommandLineInfo(ArgParser& argParser) {
    // Create dummy command-line args to make display prettier and
    // easier.
    const ArgParser::ArgRecord DummyArgs[] = {
        {"", "fmwsca: Framed, Multi-Word String Compare Analyzer",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "clu   : CLU's similarity metric generation algorithm",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "matrixFile: Use distance data stored in a matrix file",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "d2  : Use D2 (wcd-based) distance metric generation algorithm",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "d2zim : Use D2 (Zimmerman) distance metric generation algorithm",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "twoPassD2: Use two-pass asymmetric/bounded symmetric D2",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "twoPassD2adapt: Use two-pass asymmetric/bounded symmetric D2, adaptive mode",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "baton: Baton-based similarity metric generation algorithm",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "primes: Primes-based similarity metric",
         NULL, ArgParser::INFO_MESSAGE},
#ifdef HAVE_CUDA
        {"", "d2cuda: CUDA-based D2 distance metric generation algorithm",
         NULL, ArgParser::INFO_MESSAGE},
#endif
        {"", "align: Alignment-based distance metric generation",
         NULL, ArgParser::INFO_MESSAGE},        
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(DummyArgs);
}

ESTAnalyzer*
ESTAnalyzerFactory::create(const std::string& name) {    
    if (name == "fmwsca") {
        return new FMWSCA();
    } else if (name == "clu") {
        return new CLU();
    } else if (name == "matrixFile") {
        return new MatrixFileAnalyzer();
    } else if (name == "d2") {
        return new D2();
    } else if (name == "d2zim") {
        return new D2Zim();        
    } else if (name == "twoPassD2") {
        return new TwoPassD2();
    } else if (name == "oldTwoPassD2") {
        return new OldTwoPassD2();
    } else if (name == "twoPassD2adapt") {
        return new AdaptiveTwoPassD2();
    } else if (name == "baton") {
        // return new BatonAnalyzer();
    } else if (name == "primes") {
        return new PrimesESTAnalyzer();
    } else if (name == "align") {
        return new AlignmentAnalyzer();
    } else if (name == "tmalign") {
        return new TMAlignAnalyzer();
    }

#ifdef HAVE_CUDA
    if (name == "d2cuda") {
        return new D2Cuda();
    }
#endif

    // invalid analyzer name!
    std::cerr << "Invalid analyzer name." << std::endl;
    return NULL;
}

#endif
