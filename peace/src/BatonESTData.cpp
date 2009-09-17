#ifndef BATON_EST_DATA_CPP
#define BATON_EST_DATA_CPP

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

#include "BatonESTData.h"
#include <cstring>

BatonESTData::BatonESTData(const char* estSequence, const int sectionWidth,
                           const int nMer, const int maxMer) {
    maxMerValue = maxMer;
    section = sectionWidth;
    size = strlen(estSequence) - (nMer - 1);
    nSections = (size + (nMer - 1)) / section;
    int hangover = size - (nSections * section);
    section += hangover/nSections;
    nSections = 2*nSections - 1;
    intEST = new int[size];
    int v = 0, exp = 1, cv = 0;
    for (int i = 0; i < size; i++) {
        v = 0;
        exp = 1;
        for (int j = 0; j < nMer; j++) {
	    char c = estSequence[i+j];
	    cv = codeBaseLetter(c);
	    if (cv == -1) {
                v = maxMerValue;
                //printf("Error: something other than ACTG\n");
                break;
	    }
	    v += cv*exp;
	    exp*=4;
        }
        intEST[i] = v;
    }
    batons = NULL;
}

BatonESTData::~BatonESTData() {
    if (intEST != NULL) {
	delete[] intEST;
    }
    if (batons != NULL) {
	delete[] batons;
    }
}

std::vector<Baton>*
BatonESTData::getSectionBatons(int nSec) {
    int beg = nSec*section/2;
    int end = beg+section;
    std::vector<Baton>* sectionBatons = new std::vector<Baton>[maxMerValue];
    for (int k = 0; k < maxMerValue; k++) {
        if (!batons[k].empty()) {
            sectionBatons[k] = std::vector<Baton>();
            for (size_t i = 0; i < batons[k].size(); i++) {
                int bBeg = batons[k][i].getBegin();
                if (bBeg >= beg && bBeg < end) {
                    sectionBatons[k].push_back(batons[k][i]);
                }
            }
        }
    }
    return sectionBatons;
}

void
BatonESTData::resetBatons() {
    for (int i = 0; i < maxMerValue; i++) {
        if (!batons[i].empty()) {
            for (size_t j = 0; j < batons[i].size(); j++) {
                batons[i][j].sectionCount = -1;
            }
        }
    }
}

#endif
