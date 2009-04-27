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
}

ESTCodec::~ESTCodec() {
    // Nothing to be done for now.
}

#endif
