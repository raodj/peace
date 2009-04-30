#ifndef FW_ANALYZER_CPP
#define FW_ANALYZER_CPP

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

#include "FWAnalyzer.h"
#include "ResultLog.h"
#include "EST.h"

#include <algorithm>
#include <time.h>

// The static instance variables for command line arguments.
int FWAnalyzer::frameSize = 20; // 20 is default from CLU implementation
int FWAnalyzer::wordSize  = 6;  // 6  is default from CLU implementation.

// The common set of arguments for all FW EST analyzers
arg_parser::arg_record FWAnalyzer::commonArgsList[] = {
    {"--frame", "Frame size (in base pairs, default=20)",
     &FWAnalyzer::frameSize, arg_parser::INTEGER},
    {"--word", "Word size (in base pairs, default=6)",
     &FWAnalyzer::wordSize, arg_parser::INTEGER},    
    {NULL, NULL}
};

FWAnalyzer::FWAnalyzer(const std::string& analyzerName, const int refESTidx,
                       const std::string& outputFile)
    : ESTAnalyzer(analyzerName, refESTidx, outputFile) {
    // Nothing else to be done here for now.
}

FWAnalyzer::~FWAnalyzer() {
    // Clear out all EST information.
    EST::deleteAllESTs();
}

void
FWAnalyzer::showArguments(std::ostream& os) {
    ESTAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(FWAnalyzer::commonArgsList);
    os << "Options for " << analyzerName << " are:\n";
    os << ap;
}

bool
FWAnalyzer::parseArguments(int& argc, char **argv) {
    arg_parser ap(FWAnalyzer::commonArgsList);
    ap.check_args(argc, argv, false);
    // Now let the base class do processing and return the result.
    if (!ESTAnalyzer::parseArguments(argc, argv)) {
        // There are invalid/missing arguments!
        return false;
    }
    if (frameSize < 0) {
        // Necessary parameters not specified.
        std::cerr << analyzerName
                  << ": Frame size must be greater than zero "
                  << "(use --frame option)\n";
        return false;
    }
    if ((wordSize < 0) || (wordSize > frameSize)) {
        // Something is wrong with word size information.
        std::cerr << analyzerName
                  << ": Word size (greater than frame size) not specified? "
                  << "(use --word option)\n";
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
FWAnalyzer::dumpESTList(const std::vector<EST*>& estList,
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
FWAnalyzer::initialize() {
    if (!loadFASTAFile(estFileName)) {
        // Loading EST's was not successful.  Can't do much further.
        return 1;
    }
    // Now initialize our heuristic chain to prep heuristics for
    // analysis.
    if (chain->initialize()) {
        // Error occured during initialization. Bail out.
        return 2;
    }
    // Initializawtion successful.
    return 0;
}

int
FWAnalyzer::setReferenceEST(const int estIdx) {
    // Clear out the reference EST information.
    referenceFrame = "";
    // Setup the new reference EST idx.
    refESTidx      = estIdx;
    // Obtain list of ESTs for further processing.
    std::vector<EST*>& estList = EST::getESTList();
    if (refESTidx >= (int) estList.size()) {
        // Invalid reference EST id.
        std::cerr << "Reference EST index is greater than number of ESTs.\n"
                  << "Cannot continue further processing.\n";
        return 2;
    }

    // Obtain the reference frame from end of the reference EST sequence.
    referenceFrame = getFrame(estList[refESTidx], false);
    // Return 0 to indicate things proceeded successfully.
    return 0;
}

float
FWAnalyzer::analyze(const int estIdx) {
    if (estIdx == refESTidx) {
        // Comparing reference to itself simply results in 0 all the
        // time, in all FWAnalyzers
        return 0;
    }

    const EST* est = EST::getESTList()[estIdx];
    return  analyze(referenceFrame, getFrame(est), wordSize);
}

int
FWAnalyzer::analyze() {
    int result = initialize();
    if (result != 0) {
        // Error occured during initialization...
        return result;
    }
    
    // Set the reference EST.
    if ((result = setReferenceEST(refESTidx)) != 0) {
        // Error occured when setting up the reference EST...
        return result;
    }

    // Obtain list of ESTs for further processing.
    std::vector<EST*>& estList = EST::getESTList();    
    // Perform the core analysis and track total similarity to compute
    // mean value at the end of the loop.
    double total = 0;
    for(int id = 0; (id < (int) estList.size()); id++) {
        // Analyze and update similarity value.
        const float similarity = analyze(id);
        estList[id]->setSimilarity(similarity);
        total += similarity;
    }
    // Compute mean and then compute variance by reiterating over the
    // list and summing up the square of deviations.
    const double mean = total / estList.size();
    double deviations = 0;
    for(int id = 0; (id < (int) estList.size()); id++) {
        if (id != refESTidx) {
            const double diff = estList[id]->getSimilarity() - mean;
            deviations += (diff * diff);
        }
    }
    const double variance = deviations / estList.size();

    // Note the reference EST for future reference.
    const EST* refEST = estList[refESTidx];
    //Sort the resulting ESTs.
    std::sort(estList.begin(), estList.end(), EST::LessEST());

    // ----------- Start logs and dump header.----------------
    ResultLog log(outputFileName, htmlLog);
    // Dump header information.
    dumpHeader(log, mean, variance);
    // Dump information about each EST.
    dumpESTList(estList, refEST, log);
    
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
    log.reportLine("%sEST Analyzer used: %s%s", Tags[0], analyzerName.c_str(),
                   Tags[1]);
    log.reportLine("EST data read from file: %s%s%s", Tags[2], estFileName,
                   Tags[3]);
    log.reportLine("");

    // Now report the parameters used for analysis.  This information
    // is very useful for comparisons later on.
    log.report("Frame Size: %d", "  ", "Word Size: %d", frameSize, wordSize);
    log.report("Reference EST index: %d", "  ",
               "Number of ESTs: %d", refESTidx, EST::getESTList().size());
    log.report("Similarity Metric mean: %lf", "  ",
               "Similarity Metric variance: %lf", mean, variance);
    log.endTable();
    log.reportLine("");
}

float 
FWAnalyzer::analyze(const std::string&,
                    const std::string&,
                    const int) {
    return 0;
}

#endif
