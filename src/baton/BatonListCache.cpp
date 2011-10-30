#ifndef BATON_LIST_CACHE_CPP
#define BATON_LIST_CACHE_CPP

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

#include "BatonListCache.h"
#include "BatonList.h"
#include "ArgParser.h"
#include "Utilities.h"

BatonListCache::BatonListCache() {
    windowSize = 100;
    nMerSize   = 3;
    blCache.clear();
}

BatonListCache::~BatonListCache() {
    // Free up the list entries in the cache.
    std::map<int, BatonList*>::iterator entry = blCache.begin();
    while (entry != blCache.end()) {
        delete entry->second;
        entry++;
    }
    // Clear out the map itself
    blCache.clear();
}

void
BatonListCache::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord ArgsList[] = {
        {"--nMers", "The number of bases (1 to 3) to be used for baton ends",
         &nMerSize, ArgParser::INTEGER},
        {"--window", "The average size (in nt) of a window for matching identical batons",
         &windowSize, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(ArgsList);
}

BatonList*
BatonListCache::getBatonList(const EST* est, const bool getRC) {
    ASSERT( est != NULL );
    // Determine index position for the cached entry.  The normal and
    // RC baton lists for an cDNA fragment with index k are stored in
    // at position k*2 and k*2+1 respectively.
    const int cachePos = est->getID() * 2 + (getRC ? 1 : 0);
    // The following pointer is populated in the code below and
    // returned.
    BatonList* retVal = NULL;
    // Check to see if we already have an entry.
    std::map<int, BatonList*>::iterator entry = blCache.find(cachePos);
    // If we don't have an entry at the corresponding position, then
    // we need to construct one now.
    if (entry == blCache.end()) {
        // Don't have a baton list here. So create one and save it for
        // further use.
        retVal = new BatonList(est, nMerSize, getRC, windowSize);
        // Put newly created entry into cache for future use
        blCache[cachePos] = retVal;
    } else {
        // Use existing entry from the cache.
        retVal = entry->second;
    }
    // When control drops here we always have a valid entry
    return retVal;
}

void
BatonListCache::unpopulate(const EST* est) {
    ASSERT( est != NULL );
    // Determine index position for the cached entry.  The normal and
    // RC baton lists for an cDNA fragment with index k are stored in
    // at position k*2 and k*2+1 respectively.
    const int cachePos = est->getID() * 2;
    // Check to see if we have regular entry
    std::map<int, BatonList*>::iterator regEntry = blCache.find(cachePos);
    if (regEntry != blCache.end()) {
        delete regEntry->second;
        blCache.erase(regEntry);
    }
    // Check and remove reverse-complementary (rc) entry
    std::map<int, BatonList*>::iterator rcEntry = blCache.find(cachePos + 1);
    if (rcEntry != blCache.end()) {
        delete rcEntry->second;
        blCache.erase(rcEntry);
    }
}

#endif
