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
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include <cstring>
#include <string>
#include <memory>

#ifndef _WINDOWS

// The following compile guard is needed so that the hash map can be
// used for a test when running configure (because when configure is
// running, there is no config.h)
#ifdef REGULAR_COMPILE
#include "config.h"
#endif

#ifdef HAVE_TR1_UNORDERED_MAP

// We must be compiling with gcc 4.2+ in which hash map has been moved
// to unordered_map in preperation for the Cx0 standard compliance. So
// use the latest and greatest standard
#include <tr1/unordered_map>
#define GLIBC_NAMESPACE std::tr1
#define HashMap std::tr1::unordered_map
#define Hash    std::tr1::hash

#else // Not gcc 4.2+ (but linux)

#ifdef HAVE_EXT_HASH_MAP
// Use standard hash map from the extended name space.  With GCC 3.2
// and above the hash map data structure has been moded to a extended
// directory under a non-std namespace.  The following defines
// attempts to put the hash map in the standard context making it much
// easier to work with.
#include <ext/hash_map>
#define GLIBC_NAMESPACE __gnu_cxx
#define HashMap __gnu_cxx::hash_map
#define Hash    __gnu_cxx::hash

#else
// No TR1 and no EXT. Hopefully HAS_HASH_MAP is true
#include <hash_map>
#define GLIBC_NAMESPACE std
#define HashMap std::hash_map
#define Hash    std::hash

#endif

#endif

// Hasher for std::string. This hash function is needed to use std::string
// objects as key value in hash_map
struct StringHasher  {
    inline size_t operator()(const std::string& s) const {
        Hash<const char*> hasher;
        return hasher(s.c_str());

    }
};

/** \typedef hash_map<std::string, int, StringHasher> StringIntMap

	\brief A typedef for a hash map whose key is std::string and
    contains integers.
    
    The following typedef provides a short cut for using a hash map
    whose key is a std::string and contains integers.
*/
typedef HashMap<std::string, int> StringIntMap;

#else // Windows code path begins 

#include <hash_map>

// In windows the following mappings are used.
#define HashMap stdext::hash_map
#define Hash    stdext::hash

/** String comparison structure for const char *.

    The following structure essentially provides the comparison
    operator needed by the hash_map for comparing hash_key values.
    This structure specifically provides comparison for hash_map's
    whose key values are C strings.
*/
struct LessString {
    inline bool operator()(const std::string& s1, const std::string& s2) const {
        return s1.compare(s2) < 0;
    }
};

/** \typedef A hash_map<std::string, int>

    A typedef for a hash map whose key is std::string and contains
    integers.
    
    The following typedef provides a short cut for using a hash map
    whose key is a std::string and contains integers.
*/
typedef HashMap<std::string, int, stdext::hash_compare<std::string, LessString> > StringIntMap;

#endif

/** String comparison structure for HashMap.
    
    The following structure essentially provides the functor for
    comparison operator needed by the hash_map for comparing hash_key
    values.  This structure specifically provides comparison for
    hash_maps' whose key values are standard c strings.
*/
struct EqualStr {
    bool operator()(const char* s1, const char* s2) const {
        return strcmp(s1, s2) == 0;
    }
};

/** Integer comparison structure for HashMap.
    
    The following structure essentially provides the functor for
    comparion operator needed by the hash_map for comparing hash_key
    values.  This structure specifically provides comparison for
    hash_map's whose key values are standard integers.
*/
struct EqualInteger {
    inline bool operator()(const int i1, const int i2) const {
        return (i1 == i2);
    }
};

#ifdef ICC
namespace __gnu_cxx {
    template<> struct hash<std::string> {
        inline size_t operator()(const std::string& str) const {
            return hash<const char*>()(str.c_str());
        }
    };
}
#endif

#endif
