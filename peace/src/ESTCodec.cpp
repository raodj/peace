#ifndef EST_CODEC_CPP
#define EST_CODEC_CPP

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
//
//---------------------------------------------------------------------------

#include "ESTCodec.h"
#include "Utilities.h"

// The statically allocated array to translate characters to their
// normal encoding.
char ESTCodec::charToInt[255];

// The statically allocated array to translate characters to their
// complementary encoding.
char ESTCodec::charToIntComp[255];

// The globally unique instance of ESTCodec.
ESTCodec ESTCodec::estCodec;

ESTCodec::ESTCodec() {
    // Initialize the array CharToInt for mapping A, T, C, and G to 0,
    // 1, 2, and 3 respectively.  Initialize values for specific base
    // pairs. Leave rest intentionally uninitialized. So that way if
    // invalid entires are accessed, valgrind will hopefully report an
    // "uninitialized memory access" error.
    charToInt[(int) 'A'] = charToInt[(int) 'a'] = 0;
    charToInt[(int) 'G'] = charToInt[(int) 'g'] = 1;
    charToInt[(int) 'C'] = charToInt[(int) 'c'] = 2;
    charToInt[(int) 'T'] = charToInt[(int) 't'] = 3;
    // Now initialize the complementary array.
    charToIntComp[(int) 'A'] = charToIntComp[(int) 'a'] = 3;
    charToIntComp[(int) 'G'] = charToIntComp[(int) 'g'] = 2;
    charToIntComp[(int) 'C'] = charToIntComp[(int) 'c'] = 1;
    charToIntComp[(int) 'T'] = charToIntComp[(int) 't'] = 0;
    // Initialize pointers
    revCompTable = NULL;
}

ESTCodec::~ESTCodec() {
    // Free up all the reverse complement table that was constructed.
    HashMap<int, int*>::iterator curr = revCompTables.begin();
    while (curr != revCompTables.end()) {
        int *rcTable = curr->second;
        delete [] rcTable;
        // Onto the next entry
        curr++;
    }
    // Clear out all the entires in the hash map
    revCompTables.clear();
}

void
ESTCodec::setRevCompTable(const int wordSize) {
    if ((revCompTable = revCompTables[wordSize]) == NULL) {
        // A table does not exist. So create one
        revCompTable = addRevCompTable(wordSize);
    }
    ASSERT ( revCompTable != NULL );
}

int*
ESTCodec::addRevCompTable(const int wordSize) {
    // Create a reverse complement table with 4^wordSize entries.
    const int EntryCount = 1 << (wordSize * 2); // 4^wordSize
    // Create the translation table.
    int *rcTable = new int[EntryCount];
    // Populate the rcTable now.
    for(int entry = 0; (entry < EntryCount); entry++) {
        // Obtain complement
        int word    = ~entry; 
        // Reverse 2-bits at a time to obtain
        int rcValue = 0;
        for(int bp = 0; (bp < wordSize); bp++) {
            // Shift two words over and OR-in the 2 bits.
            rcValue = (rcValue << 2) | (word & 3);
            // Get rid of the two bits we used
            word >>= 2;
        }
        // Update the reverse-complement entry
        rcTable[entry] = rcValue;
    }
    // Add newly created reverse-complement entry to the hash map for
    // future lookups and use.
    revCompTables[wordSize] = rcTable;
    // Return newly created table back to caller as per API contract
    return rcTable;
}

#endif
