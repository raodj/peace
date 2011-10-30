#ifndef NT_DISTR_H
#define NT_DISTR_H

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

#include "ESTCodec.h"
#include "Utilities.h"

#include <vector>
#include <ostream>
#include <algorithm>

/** \file NtDistr.h

    \brief Class to encapsulate nucleotide occurrence frequencies at a
    given column position in a contig.

    This file contains the definition for a simple class used to
    encapsulate nucleotide occurrence frequencies.  This class is used
    by both the assembly sub-system and output sub-system constitution
    PEACE.  Consequently, it is in the shared area of PEACE.  Data is
    maintained in a coded form for memory efficiency.  Encoding and
    decoding is performed via the ESTCodec class.
*/

/** 
    <p>A simple class to encapsulate the nucleotide occurrence
    frequencies for each nucleotide/base position in this
    contig. Each value indicates the number of occurences of a
    specific nucleotide -- that is, it counts number of 'A', 'T',
    'C', 'G', or 'N' bases that occur in a single column of the
    contig to facilitate creation of the consensus sequence.</p>
    
    <p>The information is stored as an array to facilitate
    processing the distribution information.  The array indexes for
    a given base pair is computed via ESTCodec::encode()
    method.</p>
    
    \warning This class is sent over the wire via MPI
    calls. Consequently, this class cannot have any pointers or
    other objects.  All members must be primitive data types only.
*/
class NtDistr {
    friend std::ostream& operator<<(std::ostream&os, const NtDistr&);
public:
    /** Default constructor.

        The default constructor initializes the distribution of
        the nucleotides to zero.
    */
    inline NtDistr() {
        freq[0] = freq[1] = freq[2] = freq[3] = freq[4] = 0;
    }
        
    /**
       Helper method to encode and add one more occurrence of a
       given nucleotide to the set of bases constituting a given
       column in this contig.

       \param[in] nucleotide The base pair to be accounted for an
       occurence.  The nucleotide is encoded (to an integer in the
       range 0 to 4) via the ESTCodec::encode() method.
    */
    inline void add(const char nucleotide) {
        freq[(int) ESTCodec::encode(nucleotide)]++;
    }

    /** Add frequency occurrences of all nucleotides from another
        object.

        This is a convenience method that can be used to add
        occurrence frequencies of various nucleotides from another
        NtDistr object to corresponding entries in this object.

        \param[in] other The other NtDistr object whose occurence
        frequencies are to be added to corresponding occurrence
        frequencies in this class.
    */
    void add(const NtDistr& other) {
        freq[0] += other.freq[0];
        freq[1] += other.freq[1];
        freq[2] += other.freq[2];
        freq[3] += other.freq[3];
        freq[4] += other.freq[4];            
    }

    /** Get the consensus nucleotide for this entry.

        This method is a convenience method that can be used to
        obtain the consensus nucleotide. The consensus nucleotide
        is the base that has the highest occurrence frequency.

        \return The nucleotide (aka base) that has the highest
        occurrence frequency. Valid return values is one of \c
        ATCG-.
    */
    inline char getConsensus() const {
        const int freqBase = std::max_element(freq, freq+4) - freq;
        ASSERT( freq[freqBase] > 0 );
        return ESTCodec::decode(freqBase);
    }

    /** Determine if all occurrence frequencies are zero.

        This method returns true if all occurrence frequencies for
        this entry are zero.

        \return Returns true if all occurrence frequencies are
        zero (indicating not a single nucleotide is present) or
        false otherwise.
    */
    inline bool empty() const {
        return ((freq[0] + freq[1] + freq[2] + freq[3] + freq[4]) == 0);
    }
		
    /**
       Array to track the occurrence frequency of each base
       pair. Index into this array are computed via the
       ESTCodec::encode() method.
    */
    int freq[5];
};

/** \typedef std::vector<NtDistr> NtDistrList
        
    \brief A vector that contains a list of NtDistr objects.
        
    This typedef provides a convenient short cut to refer to a
    vector that is used to hold a list of NtDistr objects.  This
    vector is used in more than one spot in this class
    code. Consequently, it has been defined here to serve as a
    common data structure.
*/
typedef std::vector<NtDistr> NtDistrList;

/** \fn std::ostream& operator<<(std::ostream&, const NtDistr&)
    
    Insertion operator to stream NtDistr information to a given output
    stream.  This method provides a convenient mechanism to dump the
    complete NtDistr information for debugging purposes.  It prints
    the occurrence frequencies and the final consensus nucleotide on
    one line in the form: <tt>A[A=5,T=2,C=0,G=1]</tt>
*/
extern std::ostream& operator<<(std::ostream&, const NtDistr&);

#endif
