#ifndef HEURISTIC_FACTORY_CPP
#define HEURISTIC_FACTORY_CPP

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
// Authors: James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "HeuristicFactory.h"
#include "arg_parser.h"

void
HeuristicFactory::displayList(std::ostream &os) {
    // Create dummy command-line args to make display prettier and
    // easier.
    arg_parser::arg_record dummy_args[] = {
        //{"fmwsca", "Framed, Multi-Word String Compare Analyzer",
        // NULL, arg_parser::STRING},        
        {NULL, NULL}
    };
    arg_parser dummyParser(dummy_args);
    os << dummyParser;
}

Heuristic*
HeuristicFactory::create(const char* name, const int refESTidx,
                           const std::string& outputFileName) {
    if (name == NULL) {
        return NULL;
    }
    if (refESTidx < 0) {
        std::cerr << "A reference EST index has not been specified.\n"
                  << "Use --estIdx command line option." << std::endl;
        return NULL;
    }
    
    // if (!strcmp("fmwsca", name)) {
    //    return new FMWSCA(refESTidx, outputFileName);
    //    } else if (!strcmp("clu", name)) {
    //        return new CLU(refESTidx, outputFileName);
    //    } else if (!strcmp("matrixFile", name)) {
    //        return new MatrixFileAnalyzer(refESTidx, outputFileName);
    //    } else if (!strcmp("d2", name)) {
    //        return new D2(refESTidx, outputFileName);
    //    }

    // At present we do not have any heuristics, so any heuristic name
    // will be invalid.
    
    // invalid heuristic name!
    std::cerr << "Invalid heuristic name." << std::endl;
    return NULL;
}

#endif
