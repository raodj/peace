#ifndef BATON_EST_DATA_H
#define BATON_EST_DATA_H

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

#include "Baton.h"
#include <vector>

/** A class to encapsulate data regarding an EST and the
    characteristics of that EST needed for the baton-based analysis.

    \note This class is meant to be used only by the BatonAnalyzer.
*/
class BatonESTData {
public:
    /** The constructor.

        This constructor simply sets all the instance variables to
	corresponding values specified by the parameters.

        \param[in] estSequence Sequence of the EST.
          
        \param[in] sectWidth   Section width of the BatonAnalyzer.
    */
    BatonESTData(const char* estSequence, const int sectionWidth,
		   const int nMer, const int maxMer);
    
    /** The destructor.
    */
    virtual ~BatonESTData();

    int section;

    int size;

    int nSections;

    int* intEST;

    int maxMerValue;

    std::vector<Baton>* batons;

    std::vector<Baton>* getSectionBatons(int nSec);

    void resetBatons();


protected:
    // Currently this class does not have any protected members.
    
private:
    int codeBaseLetter(char c) const {
      switch (c) {
      case 'A':
	return 0;
      case 'C':
	return 1;
      case 'G':
	return 2;
      case 'T':
	return 3;
      default:
	return -1;
      }
    }
};

#endif
