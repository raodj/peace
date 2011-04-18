#ifndef FMWSCA_CPP
#define FMWSCA_CPP

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

#include "FMWSCA.h"
#include "ArgParser.h"

FMWSCA::FMWSCA() : FWAnalyzer("fmwsca") {
    caseSensitive = false;
    // Nothing else to be done for now.
}

FMWSCA::~FMWSCA() {
    // Nothing else to be done for now.
}

void
FMWSCA::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add its parameters first    
    FWAnalyzer::addCommandLineArguments(argParser);
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--case", "Make comparisons of nucleotides case sensitive",
         &caseSensitive, ArgParser::INVALID},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Setup valid arguments for this analyzer
    argParser.addValidArguments(ArgsList);    
}

float
FMWSCA::getMetric(const std::string& refFrame,
                 const std::string& otherFrame, const int wordSize) {
    const int LastWordPos = (int) refFrame.length() - wordSize;
    int matchCount        = 0;
    for(int wordPos = 0; (wordPos < LastWordPos); wordPos++) {
        const std::string word = refFrame.substr(wordPos, wordSize);
        int searchPos = 0;
        do {
            if ((searchPos = (int) otherFrame.find(word, searchPos)) != -1) {
                matchCount++;
                searchPos++;
            }
        } while (searchPos != -1);
    }

    return ((float) matchCount) / (LastWordPos * LastWordPos);
}

#endif
