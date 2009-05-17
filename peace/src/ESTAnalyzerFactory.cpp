#ifndef EST_ANALYZER_FACTORY_CPP
#define EST_ANALYZER_FACTORY_CPP

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

#include "ESTAnalyzerFactory.h"
#include "arg_parser.h"

#include "FMWSCA.h"
#include "CLU.h"
#include "D2.h"
#include "D2Zim.h"
#include "TwoPassD2.h"
#include "SimulationD2.h"
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
        {"d2sim", "Use special version of D2 for simulation project",
         NULL, arg_parser::STRING},
        {NULL, NULL}
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
    } else if (!strcmp("d2sim", name)) {
       return new SimulationD2(refESTidx, outputFileName);
    }
    
    // invalid analyzer name!
    std::cerr << "Invalid analyzer name." << std::endl;
    return NULL;
}

#endif
