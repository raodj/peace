#ifndef ASSEMBLER_CPP
#define ASSEMBLER_CPP

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

#include "Assembler.h"
#include "FileLoader.h"

// The static instance variables for command line arguments.
char* Assembler::fastaFileName   = NULL;
char* Assembler::sffFileName     = NULL;
bool  Assembler::noMaskBases     = false;
bool  Assembler::randomizeNbases = false;

// The common set of arguments for all assemblers
arg_parser::arg_record Assembler::commonArgsList[] = {
    {"--fastaFile", "Name of EST file (in FASTA format) to be processed",
     &Assembler::fastaFileName, arg_parser::STRING},
    {"--sffFile", "Name of SFF (Standard Flowgram Format) file to be processed",
     &Assembler::sffFileName, arg_parser::STRING},    
    {"--no-mask-bases", "Don't mask out all lower case neucleotides in reads",
     &Assembler::noMaskBases, arg_parser::BOOLEAN},
    {"--rand-N", "Randomly converts Ns in input nt to A,T,C, or G",
     &Assembler::randomizeNbases, arg_parser::BOOLEAN},
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

Assembler::Assembler(const std::string& myName, const std::string& outputFile)
    : name(myName), outputFileName(outputFile) {
    // Nothing else to be done for now.
}


Assembler::~Assembler() {
    // An empty destructor, symmetric with an empty constructor.
    // Symmetry -- Hall mark of a good design!
}

void
Assembler::showArguments(std::ostream& os) {
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(Assembler::commonArgsList);
    os << "Common options for all assemblers are:\n";
    os << ap;
}

int
Assembler::initialize() {
    if ((fastaFileName != NULL) &&
        (!FileLoader::loadFASTAFile(name, fastaFileName, !noMaskBases,
                                    randomizeNbases))) {
        // Fasta file loading failed. Can't do any assembly.
        return 1;
    }
    if ((sffFileName != NULL) &&
        (!FileLoader::loadSFFFile(name, fastaFileName, !noMaskBases,
                                  randomizeNbases))) {
        // SFF file loading failed. Can't do any assembly.
        return 2;
    }
    // Initialization successful.
    return 0;
}

bool
Assembler::parseArguments(int& argc, char **argv) {
    // Create a argument parser with custom list of arguments
    arg_parser ap(Assembler::commonArgsList);
    // Process the arguments.  There could be some leftover for other
    // classes to process and use.
    ap.check_args(argc, argv, false);
    // Check if necessary arguments have been specified for processing.
    if ((fastaFileName == NULL) && (sffFileName == NULL)) {
        // Necessary arguments have not been specified.
        std::cerr << name
                  << ": Input data file not specified (use either --fastaFile "
                  << "or --sffFile option)\n";
        return false;
    }
    // Ensure two inputs have not been specified.
    if ((fastaFileName != NULL) && (sffFileName != NULL)) {
        // Too many  arguments have been specified.
        std::cerr << name
                  << ": Specify only one input data file (use either "
                  << "--fastaFile or --sffFile option; but not both)\n";
        return false;
    }
    // Everything went well.
    return true;
}

Assembler& 
Assembler::operator=(const Assembler&) {
    return *this;
}

#endif
