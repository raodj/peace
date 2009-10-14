#ifndef CLUSTER_MAKER_CPP
#define CLUSTER_MAKER_CPP

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

#include "ClusterMaker.h"
#include "Utilities.h"

// The static instance variables for command line arguments.


// The common set of arguments for all EST analyzers
arg_parser::arg_record ClusterMaker::commonArgsList[] = {
    {NULL, NULL}
};

ClusterMaker::ClusterMaker(const std::string& myName, ESTAnalyzer *myAnalyzer,
                           const int estIdx, const std::string& outputFile) 
    : name(myName), refESTidx(estIdx), outputFileName(outputFile),
      analyzer(myAnalyzer) {
    // Nothing else to be done for now.
}

ClusterMaker::~ClusterMaker() {}

void
ClusterMaker::showArguments(std::ostream& os) {
    if (commonArgsList[0].arg_text == NULL) {
        // No valid arguments are currently used for all
        // ClusterMakers.  So skip displaying any arguments.
        return;
    }
    
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(ClusterMaker::commonArgsList);
    os << "Common options for all cluster makers are:\n";
    os << ap;
}

bool
ClusterMaker::parseArguments(int& argc, char **argv) {
    arg_parser ap(ClusterMaker::commonArgsList);
    // Process the arguments
    ap.check_args(argc, argv, false);
    // Everything went well.
    return true;
}

ClusterMaker& 
ClusterMaker::operator=(const ClusterMaker&) {
    return *this;
}

#endif
