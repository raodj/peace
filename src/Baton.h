#ifndef BATON_H
#define BATON_H

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

/** A simple class to encapsulate data regarding a baton.

    \note This class is meant to be used only by the BatonAnalyzer.
*/
class Baton {
public:
    /** The constructor.

        This constructor simply sets all the instance variables to
	corresponding values specified by the parameters.

        \param[in] beginValue
          
        \param[in] widthValue
    */
    Baton(const int beginValue, const int widthValue) :
        width(widthValue), begin(beginValue), sectionCount(-1),
        isCounted(false) {}

    // Default constructor.
    Baton() :
      begin(0), width(0), isCounted(false), sectionCount(-1) {}
    
    /** The destructor.

        This class does not contain any members that utilize dynamic
        memory to hold data.  Consequently, the destructor has no
        special tasks to perform and is merely present to adhere to
        coding conventions.
    */
    virtual ~Baton() {}

    virtual int getBegin() const
      { return begin; }

    virtual int getWidth() const
      { return width; }

    int sectionCount;

protected:
    // Currently this class does not have any protected members.
    
private:
    int begin;

    int width;

    bool isCounted;
};

#endif
