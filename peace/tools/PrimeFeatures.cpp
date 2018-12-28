#ifndef PRIME_FEATURES_CPP
#define PRIME_FEATURES_CPP

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
// Authors:   Dhananjai M. Rao          raodm@miamiOH.edu
//
//---------------------------------------------------------------------

#include "PrimeFeatures.h"
#include "Common.h"
#include "EST.h"
#include "ESTList.h"
#include "ArgParser.h"

int
PrimeFeatures::main(int argc, char *argv[]) {
    std::string estFileName; // File with ESTs to be processed.
    int numFeatures = 3;     // features/dimensions
    int atPrime     = 71;    // prime value for 'a', 't'
    int cgPrime     = 113;   // prime value for 'c', 'g'
    int refEST      = -1;    // Reference EST for distance computation
    bool showOptions= false; // Flag to indicat if help is to be displayed

    // Create the list of valid arguments to be used by the arg_parser.
    ArgParser::ArgRecord arg_list[] = {
        {"--estFile", "File with ESTs to be processed",
         &estFileName, ArgParser::STRING},
        {"--pri-heur-features", "Number of features for primes heuristic",
         &numFeatures, ArgParser::INTEGER},
        {"--pri-heur-at", "Prime value for A/T in primes heuristic",
         &atPrime, ArgParser::INTEGER},
        {"--pri-heur-cg", "Prime value for C/G in primes heuristic",
         &cgPrime, ArgParser::INTEGER},
        {"--pri-dist", "Distance value from a given read",
         &refEST, ArgParser::INTEGER},        
        {"--options", "Lists options for this tool",
         &showOptions, ArgParser::BOOLEAN},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    ArgParser ap;
    // Let base class add some defaults for us.
    Tool::addCmdLineArgs("ShowAlignment", ap);
    // Add our custom arguments
    ap.addValidArguments(arg_list);
    // Validate the command-line options.
    ap.parseArguments(argc, argv, false);    
    if (showOptions) {
        std::cout << ap << std::endl;
        // Nothing further to be done.
        return 0;
    }
    // Ensure we have all the necessary parameters.
    std::string tool = "Align";
    CHECK_ARGS(tool, estFileName.empty(), "FASTA file to load ESTs was not " \
               "specified.\nUse --estFile option\n");
    // OK, create an object for further processing.
    PrimeFeatures pf(atPrime, cgPrime);
    // Try and load the source EST file
    if (!pf.loadFastaFile(estFileName)) {
        // All the sequences could not be loaded.
        return 2;
    }
    // OK, generate the features for each read in the FASTA file
    pf.generateFeatures(numFeatures, refEST);
    // Everything went well.
    return 0;
}

// Construtor
PrimeFeatures::PrimeFeatures(const int atPrime, const int cgPrime) {
    // Set primes to be used in the base helper class
    setPrimes(atPrime, cgPrime);
}

void
PrimeFeatures::generateFeatures(int numFeatures, int refEST, std::ostream &os) {
    // Get the list of reads to process.
    ESTList& estList = getESTList();
    // If a reference EST is specified get its features so that we can
    // compute distances.
    LongVec refFeatures;
    if (refEST != -1) {
        const EST *est = estList.get(refEST, true);
        ASSERT( est != NULL );
        // Get n-dimensional features
        refFeatures = extractFeatures(est->getSequenceString(), numFeatures);
    }
    // Process each read and print n-dimensional features.    
    for (int i = 0; (i < estList.size()); i++) {
        // Get the read.
        const EST *est = estList.get(i, true);
        ASSERT( est != NULL );
        // Get n-dimensional features
        const LongVec features = extractFeatures(est->getSequenceString(),
                                                 numFeatures);
        ASSERT( numFeatures == (int) features.size() );
        // Print the information.
        os << i; //  << ", " << est->getInfo();
        for (int f = 0; (f < numFeatures); f++) {
            os << ", " << features[f];
        }
        // Print distance information from reference read, if specified
        if (refEST != -1) {
            const double dist = getDistance(refFeatures, features);
            os << ", " << dist;
        }
        os << std::endl;
    }
}

#endif
