#ifndef BATON_EST_DATA_H
#define BATON_EST_DATA_H

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
