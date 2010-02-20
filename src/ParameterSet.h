#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H

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

/** Simple class to encapsulate a set of parameters for analyzers
    and/or heuristics.
*/
class ParameterSet {
  friend class ParameterSetManager;
public:

    // Two-Pass D2 parameters
    int minLength;
    int maxLength;
    int frameSize;
    int frameShift;
    int threshold;
    int maxThreshold;

    // TV Heuristic parameters (also uses frameSize)
    int t;

    // UV Heuristic parameters
    int u;
    int wordShift;
    
    /** The destructor.
        
        The destructor for the parameter set.
    */
    virtual ~ParameterSet() {}

protected:
    /** The constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead the class
        must be instantiated via the ParameterSetManager.
    */
    ParameterSet(int min, int max, int fsz, int fsh, int thresh,
		 int maxThresh, int tIn, int uIn,  int ws) :
      minLength(min), maxLength(max), frameSize(fsz), frameShift(fsh),
      threshold(thresh), maxThreshold(maxThresh), t(tIn), u(uIn), 
      wordShift(ws) {}
    
private:

};

#endif
