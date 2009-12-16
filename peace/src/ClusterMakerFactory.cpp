#ifndef CLUSTER_MAKER_FACTORY_CPP
#define CLUSTER_MAKER_FACTORY_CPP

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

#include "ClusterMakerFactory.h"
#include "ClusterMaker.h"
#include "arg_parser.h"

#include "MSTClusterMaker.h"
#include "TransMSTClusterMaker.h"

void
ClusterMakerFactory::displayList(std::ostream &os) {
    // Create dummy command-line args to make display prettier and
    // easier.
    arg_parser::arg_record dummy_args[] = {
        {"mst", "MST-based Cluster Maker",
         NULL, arg_parser::STRING},
        {"tmst", "MST-based Cluster Maker with Transitivity",
         NULL, arg_parser::STRING},
        {NULL, NULL, NULL, arg_parser::BOOLEAN}
    };
    arg_parser dummyParser(dummy_args);
    os << dummyParser;
}

ClusterMaker*
ClusterMakerFactory::create(const char* name, ESTAnalyzer *analyzer,
                            const int refESTidx,
                            const std::string& outputFileName) {
    if (name == NULL) {
        return NULL;
    }
    if (refESTidx < 0) {
        std::cerr << "A reference EST index has not been specified.\n"
                  << "Use --estIdx command line option." << std::endl;
        return NULL;
    }
    if (analyzer == NULL) {
        std::cerr << "A EST analyzer has not been specified.\n"
                  << "Use --analyzer command line option." << std::endl;
        return NULL;
    }
    
    if (!strcmp("mst", name)) {
        return new MSTClusterMaker(analyzer, refESTidx, outputFileName);
    } else if (!strcmp("tmst", name)) {
        return new TransMSTClusterMaker(analyzer, refESTidx, outputFileName);
    }
    
    // invalid analyzer name!
    std::cerr << "Invalid analyzer name." << std::endl;
    return NULL;
}

#endif

