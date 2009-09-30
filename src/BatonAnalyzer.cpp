#ifndef BATON_ANALYZER_CPP
#define BATON_ANALYZER_CPP

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
//          James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "BatonAnalyzer.h"
#include "ResultLog.h"
#include "EST.h"

#include <algorithm>
#include <cmath>

// Some defaults (taken from TestESTAnalyzerMain.java)
// Should be parameterized.
int BatonAnalyzer::nMer = 3;
int BatonAnalyzer::sectionWidth = 100;
int BatonAnalyzer::maxErrorAllowed = 4;

// Common arguments for baton-based analyzers
// TODO - parameters
arg_parser::arg_record BatonAnalyzer::commonArgsList[] = {
    {NULL, NULL}
};

BatonAnalyzer::BatonAnalyzer(const int refESTidx,
                             const std::string& outputFile)
    : ESTAnalyzer("baton", refESTidx, outputFile) {
    refESTdata = NULL;
    otherESTdata = NULL;
    // Nothing else to be done here for now.
}

BatonAnalyzer::~BatonAnalyzer() {
    // Clear out all EST information.
    EST::deleteAllESTs();
    // Clean up memory
    if (lastPos != NULL) {
        delete[] lastPos;
    }
    if (batonArr != NULL) {
        delete[] batonArr;
     }
    if (refESTdata != NULL) {
        delete refESTdata;
    }
    if (otherESTdata != NULL) {
        delete otherESTdata;
    }
}

void
BatonAnalyzer::showArguments(std::ostream& os) {
    ESTAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(BatonAnalyzer::commonArgsList);
    os << "Options for " << analyzerName << " are:\n";
    os << ap;
}

bool
BatonAnalyzer::parseArguments(int& argc, char **argv) {
    arg_parser ap(BatonAnalyzer::commonArgsList);
    ap.check_args(argc, argv, false);
    // Now let the base class do processing and return the result.
    if (!ESTAnalyzer::parseArguments(argc, argv)) {
        // There are invalid/missing arguments!
        return false;
    }

    // Baton-specific argument validation goes here
    
    // Everything went well.
    return true;    
}

int
BatonAnalyzer::initialize() {
    if ((estFileName != NULL) && (!loadFASTAFile(estFileName))) {
        // Loading EST's was not successful.  Can't do much further.
        return 1;
    }
    // Now initialize our heuristic chain to prep heuristics for
    // analysis.
    if ((chain != NULL) && (chain->initialize())) {
        // Error occured during initialization. Bail out.
        return 2;
    }

    // Initialize important variables/structures
    maxMerValue = (int) pow(4.0, nMer);
    lastPos = new int[maxMerValue];
    batonArr = new std::vector<Baton>[maxMerValue];

    return 0;
}

int
BatonAnalyzer::setReferenceEST(const int estIdx) {
    // Clean up old reference EST, if one exists.
    if (refESTdata != NULL) {
        delete refESTdata;
        refESTdata = NULL;
    }
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

    // Now do baton-specific stuff
    refESTdata = new BatonESTData(EST::getEST(estIdx)->getSequence(), sectionWidth, nMer, maxMerValue);
    buildESTData(&refESTdata, refESTidx);   

    // Return 0 to indicate things proceeded successfully.
    return 0;
}

void
BatonAnalyzer::buildESTData(BatonESTData** data, const int
			    UNREFERENCED_PARAMETER(estIdx)) {
    // extract the batons
    resetBatonCollector();

    int* estMers = (*data)->intEST;
    for (int k = 0; k < (*data)->size; k++) {
        if (estMers[k] == maxMerValue) {
            continue;
        } else if (lastPos[estMers[k]] < 0){   // Simply add the first position
            lastPos[estMers[k]]=k;
        } else {
            int width = k-lastPos[estMers[k]];
            // If the baton size is larger than nMer, create the baton.
            if (width >= nMer) {
                batonArr[estMers[k]].push_back(Baton(lastPos[estMers[k]],
                                                     width));
                lastPos[estMers[k]]=k;
            }
        }
    }
    (*data)->batons = new std::vector<Baton>[maxMerValue];
    for (int i = 0; i < maxMerValue; i++) {
        //        if (lastPos[i] < 0) {
        //            (*data)->batons[i] = NULL;
        //        } else {
        //            (*data)->batons[i] = batonArr[i]; // will copy
        //        }
        if (lastPos[i] >= 0) {
            (*data)->batons[i] = batonArr[i]; // will copy
        }
    }
}

void
BatonAnalyzer::resetBatonCollector() {
    for (int i = 0; i < maxMerValue; i++) {
        batonArr[i].clear();
        lastPos[i] = -1;
    }
}

float
BatonAnalyzer::getMetric(const int estIdx) {
    if (estIdx == refESTidx) {
        // Comparing reference to itself - how to return max similarity?
        // This case shouldn't occur if the clustermaker is correct...
        return 100;
    }
    otherESTdata = new BatonESTData(EST::getEST(estIdx)->getSequence(), sectionWidth, nMer, maxMerValue);
    buildESTData(&otherESTdata, estIdx);

    // code here comes from ESTAnalyzer -> findSimilarityBetweenTwoEsts
    int n1 = refESTdata->nSections;
    int n2 = otherESTdata->nSections;
    int cbfLength = n1+1;
    int cbfSubLength = n2+1;
    int** commonBatonFreqs = new int*[cbfLength];
    for (int i = 0; i < cbfLength; i++) {
        commonBatonFreqs[i] = new int[cbfSubLength];
        memset(commonBatonFreqs[i], 0, sizeof(int) * cbfSubLength);
    }

    commonBatons(commonBatonFreqs, n1, n2);
    
    int highFreqLength = 0;
    int** highFreqSections = highestFrequencySection(commonBatonFreqs, 7,
                                                     &highFreqLength,
                                                     cbfLength,
                                                     cbfSubLength);
    
    for (int i = 0; i < cbfLength; i++) {
        delete [] commonBatonFreqs[i];
    }
    delete [] commonBatonFreqs;

    if (highFreqSections == NULL) {
        //printf("high freq sections == NULL, returning 0\n");
        
        // Clean up.  I suspect this will be a major issue in the future;
        // right now it is just a hack and we're not saving information.
        delete otherESTdata;
        otherESTdata = NULL;
    
        return 0; // is this right? his code returns NULL
    } else {
        int maxAlignment[3];
        maxAlignment[0] = 0;
        int alignment[3];
        
        for (int i = 0; i < highFreqLength; i++) {
            alignment[0] = alignBetweenSections(highFreqSections[i][0],
                                                highFreqSections[i][1],
                                                maxErrorAllowed,
                                                &alignment[1],
                                                &alignment[2]);
            if (alignment[0] > maxAlignment[0]) {
                maxAlignment[0] = alignment[0];
                maxAlignment[1] = alignment[1];
                maxAlignment[2] = alignment[2];
            }
        }

        //if (maxAlignment[0] != 0) {
        //    printf("%d %d %d\n", refESTidx, estIdx, maxAlignment[0]);
        //}
        
        // Clean up.  I suspect this will be a major issue in the future;
        // right now it is just a hack and we're not saving information.

        for (int i = 0; i < highFreqLength; i++) {
            delete [] highFreqSections[i];
        }
        delete [] highFreqSections;    
        delete otherESTdata;
        otherESTdata = NULL;
        
        // Finished
    
        return (float) maxAlignment[0];
    }
}

int
BatonAnalyzer::alignBetweenSections(int sect1, int sect2,
                                    int totalAllowedErrors,
                                    int* begin1, int* begin2) {
    std::vector<Baton>* refBatons = refESTdata->getSectionBatons(sect1);
    std::vector<Baton>* otherBatons = otherESTdata->getSectionBatons(sect2);
    int alignmentNumber = 0;
    for (int k = 0; k < maxMerValue; k++) {
        if (refBatons[k].empty() || otherBatons[k].empty()) {
            continue;
        }
        Baton* b2s = new Baton[otherBatons[k].size()];
        for (size_t j = 0; j < otherBatons[k].size(); j++) {
            b2s[j] = otherBatons[k][j];
        }
        for (size_t i = 0; i < refBatons[k].size(); i++) {
            Baton b1 = refBatons[k][i];
            int w1 = b1.getWidth();
            int beg1 = b1.getBegin();
            for (size_t j = 0; j < otherBatons[k].size(); j++) {
                if (w1 == b2s[j].getWidth()) {
                    int beg2 = b2s[j].getBegin();
                    int newAlignmentNumber;
                    if (quickCheckAlignment(refESTdata->intEST,
                                            otherESTdata->intEST,
                                            beg1, beg2)) {
                        newAlignmentNumber = alignESTs(refESTdata->intEST,
                                                       otherESTdata->intEST,
                                                       refESTdata->size,
                                                       otherESTdata->size,
                                                       beg1, beg2,
                                                       totalAllowedErrors);
                        // TODO:
                        // Did not implement case where totalAllowedErrors = 0
                        if (newAlignmentNumber > 15) { // Found it
                            *begin1 = beg1;
                            *begin2 = beg2;
                            delete [] b2s;
                            delete [] refBatons;
                            delete [] otherBatons;
                            return newAlignmentNumber;
                        } else {
                            if (newAlignmentNumber > alignmentNumber) {
                                alignmentNumber = newAlignmentNumber;
                            }
                        }
                    }
                }
            }
        }
        delete [] b2s;
    }
    delete [] refBatons;
    delete [] otherBatons;
    return alignmentNumber;
}

bool
BatonAnalyzer::quickCheckAlignment(int* est1, int* est2, int beg1, int beg2) {
    bool furtherCheck = true;
    if (est1[beg1+3] != est2[beg2+3]) {
        if (beg1 > 3 && beg2 > 3 && est1[beg1-3] != est2[beg2-3]) {
            furtherCheck = false;
        }
    }
    return furtherCheck;
}

int
BatonAnalyzer::alignESTs(int* est1, int* est2, int est1size, int est2size,
                         int beg1, int beg2, int nTotalMutationError) {
    int leftStep = 0;
    int rightStep = 0;
    int end1 = est1size - nMer;
    int end2 = est2size - nMer;

    while (nTotalMutationError > 0) {
        rightStep = rightExtension(nMer, est1, est2, beg1, beg2, end1, end2,
                                   rightStep);
        rightStep += nMer;
        nTotalMutationError--;
        leftStep = leftExtension(nMer, est1, est2, beg1, beg2, leftStep);
        leftStep += nMer;
        nTotalMutationError--;
        // TODO Possible bug here? could end up with < 0 if it's 1 starting
        // Test later, after I've checked with java version for same scores
        //assert(nTotalMutationError >= 0);
    }

    return leftStep+rightStep - 2*nMer;
}

int
BatonAnalyzer::rightExtension(int nMer, int* est1, int* est2, int beg1,
                              int beg2, int end1, int end2, int rStepSoFar) {
    while ((beg1+rStepSoFar) < end1 && (beg2+rStepSoFar) < end2) {
        rStepSoFar++;
        if (est1[beg1+rStepSoFar] != est2[beg2+rStepSoFar]) {
            break;
        }
    }
    return rStepSoFar;
}

int
BatonAnalyzer::leftExtension(int nMer, int* est1, int* est2, int beg1,
                              int beg2, int lStepSoFar) {
    while ((beg1-lStepSoFar) >= nMer && (beg2-lStepSoFar) >= nMer) {
        lStepSoFar++;
        if (est1[beg1-lStepSoFar] != est2[beg2-lStepSoFar]) {
            break;
        }
    }
    return lStepSoFar;
}

int**
BatonAnalyzer::highestFrequencySection(int** commonBatonFreqs, int criticalV,
                                       int* arrLength, const int cbfLength,
                                       const int cbfSubLength) {
    int* overlappingSectionCandidates[50];
    int k = 0;
    int m1 = cbfLength;
    int m2 = cbfSubLength;

    for (int i = 0; i < m1; i++) {
        for (int j = 0; j < m2; j++) {
            if (commonBatonFreqs[i][j] >= criticalV) {
                int* pair = new int[2];
                pair[0] = i;
                pair[1] = j;
                overlappingSectionCandidates[k++] = pair;
            }
        }
    }
    if (k == 0) {
        return NULL;
    } else {
        int** candidateSections = new int*[k];
        for (int i = 0; i < k; i++) {
            candidateSections[i] = overlappingSectionCandidates[i];
        }
        *arrLength = k;
        return candidateSections;
    }
}

void
BatonAnalyzer::commonBatons(int** identicalBatonNumber, int n1, int n2) {
    int s1 = (refESTdata->section)/2;  // this is the half width.
    int s2 = (otherESTdata->section)/2;
    int l1 =0;
    int beg1 =0;

    for (int m=0; m < maxMerValue; m++){
        if (refESTdata->batons[m].empty() || otherESTdata->batons[m].empty())
            continue;
        for (size_t i = 0; (i < refESTdata->batons[m].size()); i++){
            l1 = refESTdata->batons[m][i].getWidth();
            beg1 = refESTdata->batons[m][i].getBegin();
            int k1 = beg1/s1; // Section number.

            for (size_t j = 0; j < (otherESTdata->batons[m].size()); j++){
                if (l1 == otherESTdata->batons[m][j].getWidth()){
                    // if the baton widths are the same
                    // then count them in the frequencies.
                    int k2 = otherESTdata->batons[m][j].getBegin()/s2;

                    if (refESTdata->batons[m][i].sectionCount >= k2)
                        // Already counted.
                        continue ;
                    if (otherESTdata->batons[m][j].sectionCount >= k1)
                        // This prevents double-counting of a baton of b2
                        // for the same section of b1.
                        continue;

                    refESTdata->batons[m][i].sectionCount = k2;
                    otherESTdata->batons[m][j].sectionCount = k1;


                    if (k1==0 || k2==0){
                        // One or both sections are the first.
                        if (k1==0&& k2==0)
                            identicalBatonNumber[0][0]++;
                        else if ( k1 == 0){
                            if (k2 <= n2) { // not last section
                                identicalBatonNumber[0][k2]++;
                            }
                            identicalBatonNumber[0][k2-1]++;
                        }
                        else {
                            //if (k1 ==n1+1 )
                            //    printf("\n");
                            if (k1 <= n1) { // not last section
                                identicalBatonNumber[k1][0]++;
                            }
                            identicalBatonNumber[k1-1][0]++;
                        }
                    }
                    else {
                        if (k1 >= n1 || k2 >= n2){
                            // One or both section are the last.
                            if (k1 >= n1 && k2 >= n2)
                                identicalBatonNumber[k1-1][k2-1]++;
                            else if (k1 >= n1){
                                identicalBatonNumber[k1-1][k2]++;
                                identicalBatonNumber[k1-1][k2-1]++;
                            }
                            else {
                                identicalBatonNumber[k1][k2-1]++;
                                identicalBatonNumber[k1-1][k2-1]++;
                            }
                        }
                        else {
                            identicalBatonNumber[k1-1][k2-1]++;
                        }
                    }
                }
            }
        }
    }

    refESTdata->resetBatons();
    otherESTdata->resetBatons();
}

int
BatonAnalyzer::analyze() {
    int result = initialize();
    if (result != 0) {
        // Error occured during initialization...
        return result;
    }
    
    // Obtain list of ESTs for further processing.
    std::vector<EST*>& estList = EST::getESTList();
    
    // Perform the core analysis and track total similarity to compute
    // mean value at the end of the loop.
    for (int id1 = 0; (id1 < (int) estList.size()); id1++) {
        // Set the reference EST.
        setReferenceEST(id1);
        for (int id2 = id1+1; (id2 < (int) estList.size()); id2++) {
            // Analyze and update similarity value.
            getMetric(id2);
        }
    }

    // All the processing went on successfully.
    return 0;
}

#endif
