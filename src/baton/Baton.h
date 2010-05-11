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

#include <iostream>

/** \file Baton.h

    \brief Class declaration for a Baton class, that is used to
    encapsulate the necessary data.

    This class is used to encapsulate the data associated with a
    single Baton
*/

/** A single Baton entry in a BatonList.

    This class is used to encapsulate the data associated with a
    single Baton.  A Baton is a contiguous part a given cDNA fragment
    with the same <i>n</i>-mers at the beginning and end.  Here are a
    few examples of Batons:

    <ul>

    <li><tt>AAAA<i><b>TCG</b>AGTACCTAGGCT<b>TCG</b></tt> (this is a
    3-mer Baton of length 18 bp.</li>

    <li><tt>TATA<i><b>CC</b>TTACATTCA<b>CC</b></tt> (this is a
    2-mer Baton of length 13 bp.</li>


    <li><tt>AAAA<i><b>TCG</b>AGTACCTAGGCT<b>TCG</b>AGTACT<b>TCG</b></tt>
    (this is an example of two, 3-mer Batons of length 18 bp and 12
    bp.</li>
    
    </ul>

    This class merely stores the meta data information associated with
    a given <i>n</i>-mer Baton.  Specifically, this class encapsulates
    the following information:

    <ol>

    <li>The starting or beginning \i index position of the Baton in
    the cDNA fragment.</li>

    <li>The length of the the baton in number of base pairs (bp).
    This includes the <i>n</i>-mers that constitute the Baton
    ends.</li>

    <li>Instance variable to track the index of the cDNA fragment (\c
    otherEST) towards which this baton has been counted (or
    accounted-for) when searching for identifical batons. The default
    initial value is false. Whenever a baton is matched with another
    baton, this instance variable is set to \c true, indicating this
    baton has been acounted.</li>
    
    </ol>

    Here are some of the assumptions/background on how the batons
    work:

    <ul>

    <li>Note that an Baton object does not actually contain the
    <i>n</i>-mer encoding associated with its its ends.  This
    information is logically maintained by the BatonList class when it
    builds a list of similar (batons with the same <i>n</i>-mer endds)
    baton. </li>

	<li>It does not contain the actual nucleotide sequences.  This
	information is obtained from the EST class.  The index of the EST
	with which this baton is associated is stored in the BatonList
	class.</li>
	
    <ul>
*/
class Baton {
    // Insertion operator to make dumping Batons easier.
    friend std::ostream& operator<<(std::ostream&, const Baton&);  
public:
    /**
       Default constructor.

       The default constructor initializes the start index and length
       to zero.  In addition, the countedFor instance variable is set
       to \c false to indicate that this Baton has not be counted
       towards another EST when searching for matching batons.

       \note The default constructor is handy to enable use of this
       class with standard STL classes.
    */
    Baton() {
        startIndex = 0;
        length     = 0;
        accounted  = false;
    }

    /**
       Convenience constructor.

       This constructor is a convenience constructor that can be used
       to set the various instance variables to appropriate values.

       \param[in] startIndex The starting or beginning \i index
       position of the Baton in the cDNA fragment.

       \param[in] length The length of the the baton in number of base
       pairs (bp)</li>. This includes the <i>n</i>-mers that
       constitute the Baton ends.

       \param[in] accounted If this value is \c false, that means
       this baton has not already been counted (or accounted-for) when
       counting identical batons.
    */
    Baton(const int startIndex, const int length,
          const bool accounted = false) {
        this->startIndex = startIndex;
        this->length     = length;
        this->accounted  = accounted;
    }

    /** Determine the zero-based index position from where this baton
        starts.

        This method must be used to obtain the zero-based index
        position within a given cDNA fragment from where this baton
        starts. The actual cDNA fragment is maintained by the
        BatonList class (which logically contains a set of batons).

        \return The zero-based index position from where this baton
        starts.
    */
    inline int getStartIndex() const { return startIndex; }

    /** Determine the length of this baton.

        This method must be used to obtain the length of the baton.
        The length includes the <i>n</i>-mers that constitute the two
        ends of the baton.  For example, in the cDNA fragment
        <li><tt>AAAA<i><b>TCG</b>AGTACCTAGGCT<b>TCG</b></tt>, the
        <i>3</i>-mer baton (with ends <tt>TCG</tt>) has a length of 18
        nt.  Note that the actual cDNA fragment is maintained by the
        BatonList class (which logically contains a set of batons).

        \return The length of this baton (in number of nucleotides).
    */
    inline int getLength() const { return length; }

    /** Determine if this baton has already been accounted in a search.

         This method is a convenience method that can be used to track
         if a baton has been accounted for when searching and
         comparing batons from two different cDNA fragments.  The
         setAccounted() method must be used to set/reset the \c
         accounted instance variable.  This method is primarily used
         by the BatonList class.

         \return A \true value indicates this baton has already been
         accounted in a given search and should be ignored.  If it
         returns \c false, that indicates the baton has <b>not</b>
         been accounted in the current search phase.
    */
    inline bool isAccounted() const { return accounted; }

    /** Set/reset if this baton has already been accounted in a
        search.
        
        This must be used to set/reset if a baton has been accounted
        for when searching and comparing batons from two different
        cDNA fragments.  The isAccounted() method must be used to
        determine if a baton has been accounted.  This method is
        primarily used by the BatonList class.
        
        \param[in] accounted A \true value indicates this baton has
        already been accounted in the current search and should be
        ignored for remainder of the current search. At the end of the
        search this value is typically reset back to \c false via this
        method by the BatonList class.
    */
    inline void setAccounted(const bool accounted) {
        this->accounted = accounted;
    }

    /** Compares two Batons using just their lengths.

        This overloaded operator< (less than) compares two batons
        based purely on their lengths. None of the other instance
        variables are involved in the comparison.

        \note This operator is used for sorting batons based on their
        lengths in the BatonList::buildBatons() method.
        
        \param[in] rhs The Baton object to the right-hand-side of the
        comparison.

        \return This method returns \c true if the length of \c this
        baton is less than the length of the other baton.
    */
	inline bool operator<(const Baton& rhs) const {
        return (length < rhs.length);
	}
	
private:
    /**
       The starting or beginning \i index position of the Baton in the
       cDNA fragment.  This value is typically set when a Baton is
       created and is never changed during the life time of the Baton.
       The actual fragments into whose nucleotide sequence this
       variable indexes is only present in the BatonList class, that
       logically owns this class.
    */
    int startIndex;

    /**
       The length of this baton (in bp) including the <i>n</i>-mers at
       the two ends of the baton.  This value is typically set when a
       Baton is created and is never changed during the life time of
       the Baton.  The actual EST into whose nucleotide sequence this
       variable indexes is only present in the BatonList class, that
       logically owns this class.
    */    
    int length;

    /**
       This instance variable is a convenience instance variable that
       is used to track if a baton has been accounted for when
       searching and comparing batons from two different cDNA
       fragments.  This instance variable is initialized and reset
       after a search to \c false indicating it is not accounted.  A
       \true value indicates this baton has already been accounted in
       a given search and should be ignored.  Its value is
       appropriately changed by the setAccounted() method in this
       class.
    */
    bool accounted;
};

/** \fn std::ostream& operator<<(std::ostream&, const Baton&)

    Insertion operator to stream Baton information to a given output
    stream.  This method provides a convenient mechanism to dump the
    complete Baton information for debugging purposes.

    \param[out] os The output stream to which the baton information is
    to be written.

    \param[in] baton The baton whose information is to be displayed.

    \return As the generalized API contract for insertion operators,
    this method returns os (first parameter)
*/
extern std::ostream& operator<<(std::ostream& os, const Baton& baton);

#endif
