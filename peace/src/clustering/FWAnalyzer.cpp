#ifndef FW_ANALYZER_CPP
#define FW_ANALYZER_CPP

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

#include "ArgParser.h"
#include "FWAnalyzer.h"
#include "ResultLog.h"
#include "ESTList.h"

#include <algorithm>
#include <time.h>

FWAnalyzer::FWAnalyzer(const std::string& analyzerName)
    : ESTAnalyzer(analyzerName) {
    // Initialize instance variables
    frameSize = 100;
    wordSize  = 6;
}

FWAnalyzer::~FWAnalyzer() {
    // Nothing to be done
}

void
FWAnalyzer::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add common parameters.
    ESTAnalyzer::addCommandLineArguments(argParser);    
    const ArgParser::ArgRecord CommonArgsList[] = {
        {"--frame", "Frame size (in base pairs)",
         &frameSize, ArgParser::INTEGER},
        {"--word", "Word size (in base pairs)",
         &wordSize, ArgParser::INTEGER},    
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(CommonArgsList);
}

bool
FWAnalyzer::initialize() {
    // Validate command-line parameters first.
    if (frameSize < 0) {
        // Necessary parameters not specified.
        std::cerr << getName()
                  << ": Frame size must be greater than zero "
                  << "(use --frame option)\n";
        return false;
    }
    if ((wordSize < 0) || (wordSize > frameSize)) {
        // Something is wrong with word size information.
        std::cerr << getName()
                  << ": Word size (greater than frame size) not specified? "
                  << "(use --word option)\n";
        return false;
    }
    // OK. basic stuff seems fine. Let base class initialize first
    if (!ESTAnalyzer::initialize()) {
        // Base class initialization failed.
        return false;
    }
    // Everything went well.
    return true;    
}

std::string
FWAnalyzer::getFrame(const EST* est, bool start) {
    std::string result(est->getSequence());
    if (start) {
        return result.substr(0, frameSize);
    }
    // Get base pair from end of est sequence.
    return result.substr(result.size() - frameSize);
}

void
FWAnalyzer::dumpEST(ResultLog& log, const EST* est, const bool isReference) {
    const std::string frame = getFrame(est, !isReference);

    if (isReference) {
        const char* start=(htmlLog ? "<font face=\"courier\" color=red>" : "");
        const char* end  =(htmlLog ? "</font>" : "");
        log.report(est->getInfo(), "%s%s%s", "n/a", start, frame.c_str(), end);
    } else {
        const char* start = (htmlLog ? "<font face=\"courier\">" : "");
        const char* end   = (htmlLog ? "</font>" : "");
        log.report(est->getInfo(), "%s%s%s", "%f",
                   start, frame.c_str(), end, est->getSimilarity());
    }
}

void
FWAnalyzer::dumpESTList(const ESTList& estList,
                        const EST* refEST,
                        ResultLog& log) {
    // Start up a Table in the log.
    const char* Titles[] = {"Name", "EST Frame", "Metric", NULL};
    log.startTable(Titles);
    // Dump each est out.
    for(int id = 0; (id < (int) estList.size()); id++) {
        if (id % 20 == 0) {
            // Dump reference EST out.
            dumpEST(log, refEST, true);
        }
        if (estList[id] == refEST) {
            // Don't dump reference EST again;
            continue;
        }
        // Dump a non-Reference EST out.
        dumpEST(log, estList[id], false);
    }
}


int
FWAnalyzer::setReferenceEST(const EST* est) {
    ASSERT( est != NULL );
    // Clear out the previous reference EST information.
    referenceFrame = "";
    // Setup the new reference EST.
    refEST  = est;
    // Obtain list of ESTs for further processing.
    if (refEST->getID() >= estList->size()) {
        // Invalid reference EST id.
        std::cerr << "Reference EST index is greater than number of ESTs.\n"
                  << "Cannot continue further processing.\n";
        return 2;
    }
    
    // Obtain the reference frame from end of the reference EST sequence.
    referenceFrame = getFrame(refEST, false);
    // Return 0 to indicate things proceeded successfully.
    return 0;
}

float
FWAnalyzer::getMetric(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    if (otherEST->getID() == refEST->getID()) {
        // Comparing reference to itself simply results in 0 all the
        // time, in all FWAnalyzers
        return 0;
    }
    return  analyzeFrame(referenceFrame, getFrame(otherEST), wordSize);
}

int
FWAnalyzer::analyze() {
    int result = initialize();
    if (result != 0) {
        // Error occured during initialization...
        return result;
    }
    // Set the reference EST.
    EST* refEST = estList->get(0, true);
    if ((result = setReferenceEST(refEST)) != 0) {
        // Error occured when setting up the reference EST...
        return result;
    }

    // Perform the core analysis and track total similarity to compute
    // mean value at the end of the loop.
    double total = 0;
    for(int id = 0; (id < estList->size()); id++) {
        // Analyze and update similarity value.
        EST* est = estList->get(id, true);
        const float similarity = ESTAnalyzer::analyze(est);
        est->setSimilarity(similarity);
        total += similarity;
    }
    // Compute mean and then compute variance by reiterating over the
    // list and summing up the square of deviations.
    const double mean = total / estList->size();
    double deviations = 0;
    for(int id = 0; (id < estList->size()); id++) {
        if (id != refEST->getID()) {
            const double diff = estList->get(id)->getSimilarity() - mean;
            deviations += (diff * diff);
        }
    }
    const double variance = deviations / estList->size();
    //Sort the resulting ESTs.
    // std::sort(estList.begin(), estList.end(), EST::LessEST());

    // ----------- Start logs and dump header.----------------
    ResultLog log(outputFileName, htmlLog);
    // Dump header information.
    dumpHeader(log, mean, variance);
    // Dump information about each EST.
    dumpESTList(*estList, refEST, log);
    
    // All the processing went on successfully.
    return 0;
}

void
FWAnalyzer::dumpHeader(ResultLog& log, const double mean,
                       const double variance) {
    const char   *HTMLTags[] = {"<b>", "</b>", "<i>", "</i>", "<u>", "</u>"};
    const char   *TextTags[] = {"*  ", "  *", "", "", "_", "_"};
    const char   **Tags      = (htmlLog ? HTMLTags : TextTags);
    const char   *Title      = (htmlLog ? "E S T &nbsp;&nbsp;&nbsp; A N A L Y S I S &nbsp;&nbsp;&nbsp; R E P O R T" : "E S T     A N A L Y S I S     R E P O R T");
    // Get curren time for reporting
    char  now_str[128];
    getTime(now_str);

    // First log some general information one, line after another with
    // some HTML decorations to make text bold or underlined.  In text
    // mode, many of the HTML decorations are absent.
    log.reportLine("");
    log.reportLine("%s %s %s", Tags[0], Title, Tags[1]);
    log.reportLine("Analysis conducted on %s", now_str);
    log.reportLine("%sEST Analyzer used: %s%s", Tags[0], getName().c_str(),
                   Tags[1]);
    log.reportLine("EST data read from file: %s%s%s", Tags[2], "", // estFileName,
                   Tags[3]);
    log.reportLine("");

    // Now report the parameters used for analysis.  This information
    // is very useful for comparisons later on.
    log.report("Frame Size: %d", "  ", "Word Size: %d", frameSize, wordSize);
    log.report("Reference EST index: %d", "  ",
               "Number of ESTs: %d", refEST->getID(), estList->size());
    log.report("Similarity Metric mean: %lf", "  ",
               "Similarity Metric variance: %lf", mean, variance);
    log.endTable();
    log.reportLine("");
}

float 
FWAnalyzer::analyzeFrame(const std::string&,
                         const std::string&,
                         const int) {
    return 0;
}

#endif
