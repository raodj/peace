#ifndef HASH_MAP_H
#define HASH_MAP_H

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OH.
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

#include "config.h"

#ifdef HAVE_TR1_UNORDERED_MAP

// We must be compiling with gcc 4.2+ in which hash map has been moved
// to unordered_map in preperation for the Cx0 standard compliance. So
// use the latest and greatest standard
#include <tr1/unordered_map>
#define GLIBC_NAMESPACE std::tr1
#define HashMap std::tr1::unordered_map
#define Hash    std::tr1::hash

#else // Not gcc 4.2+

// Use standard hash map from the extended name space.  With GCC 3.2
// and above the hash map data structure has been moded to a extended
// directory under a non-std namespace.  The following defines
// attempts to put the hash map in the standard context making it much
// easier to work with.

#include <ext/hash_map>
#define GLIBC_NAMESPACE __gnu_cxx
#define HashMap __gnu_cxx::hash_map
#define Hash    __gnu_cxx::hash

#endif

#include <cstring>

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

/** \typedef A hash_map<std::string, int>

    A typedef for a hash map whose key is std::string and contains
    integers.
    
    The following typedef provides a short cut for using a hash map
    whose key is a std::string and contains integers.
*/
typedef HashMap<std::string, int> StringIntMap;

#endif
