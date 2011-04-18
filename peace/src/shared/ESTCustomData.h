#ifndef EST_CUSTOM_DATA_H
#define EST_CUSTOM_DATA_H

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

/** A dummy class that serves as a place holder.

    <p> This is a standard decorator pattern that serves as a place holder
    for associating some user defined values with a given EST.
    Another class can be derived from this class to hold some
    information pertinent to a given EST and associated with an EST
    via the EST::setCustomData() and EST::getCustomData() methods. </p>

    Note that this class does not contain any data members and is
    merely a decorator.  The constructor and destructor are present
    merely to adhere to coding conventions.
*/
class ESTCustomData {
public:
    /** The default constructor.

        The default constructor does nothing and is present merely to
        adhere to coding conventions.
    */
    ESTCustomData() {}

    /** The destructor.

        The destructor does nothing and is present merely to adhere to
        coding conventions.
    */    
    virtual ~ESTCustomData() {}
    
protected:
    // This class has no protected members.
    
private:
    // This class has no private members.
};

#endif
