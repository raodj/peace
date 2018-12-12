#ifndef HASH_MAP_H
#define HASH_MAP_H

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
// Authors:   Dhananjai M. Rao          raodm@miamiOH.edu
//
//---------------------------------------------------------------------

#include <cstring>
#include <string>
#include <memory>
#include <unordered_map>

/** As of Dec 2018, we now require the use of c++11 to compile PEACE.
    With the C++11 standards, we have a consistent unordered_map that
    we can use.  Accordingly, this header has been cleaned-up to use
    unordered_map on all (Windows, Linux, and Mac) platforms
*/
#define HashMap std::unordered_map

/** \typedef A hash_map<std::string, int>

    A typedef for a hash map whose key is std::string and contains
    integers.
    
    The following typedef provides a short cut for using a hash map
    whose key is a std::string and contains integers.
*/
using StringIntMap = std::unordered_map<std::string, int>;

/** \typedef A hash_map<std::string, std::string>

    A typedef for a hash map whose key is std::string and contains
    std::string objects.
    
    The following typedef provides a short cut for using a hash map
    whose key is a std::string and contains std::string objects..
*/
using StringStringMap = std::unordered_map<std::string, std::string>;

#endif
