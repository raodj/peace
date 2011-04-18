#ifndef SEQUENCE_ALIGNER_H
#define SEQUENCE_ALIGNER_H

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

#include "Component.h"
#include <cstdlib>

// Forward declarations to keep compiler fast and happy
class AlignmentInfo;
class BatonListCache;
class EST;

/** \file SequenceAligner.h

	\brief Declearation for a base class that can be used to pairwise
	align a given nucleotide sequence with another reference sequence.

	This class contains a generic base class definition for a helper
	class that is used by baton assemblers to pairwise-align a given
	nucleotide sequence with another reference sequence.
*/

/** \brief The base class for different types of nucleotide sequence
    aligners.

    This class serves as the base class for various types of sequence
    aligners that are used in the system.  A sequence aligner
    essentially performs the task of aligning a given nucleotide
    sequence with another reference sequence.  The base class does not
    stipulate if the SequenceAligner performs \e local (such as:
    Smith-Waterman), \e global (such as: Needleman-Wunsch) or a \e
    hybrid alignment (aka "glocal" methods).

    \note This class extended Component to merely reuse the command
    line argument processing features.  However, the run() method in
    this class must not be used (and if run() method is invoked it
    reports an error and aborts execution).    
*/
class SequenceAligner : public Component {
public:
    /** The destructor.

        The destructor is merely present to adhere to coding
        conventions.  Currently, this does not have any special tasks
        to perform.
    */
    virtual ~SequenceAligner();

    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are aligned
        via a call to the align(EST *) method.  This method currently
        saves the index in the \c refESTidx instance variable for
        further look up.  Next, it builds the normal BatonList for the
        reference EST.  The BatonList is stored in the blCache (read
        as baton list cache) instance variable for future reference.

        \note This method must be called only after the \c
        initialize() method is called.

		\param[in] est A pointer to an immutable EST object to be used
		as the reference.  The reference EST will be subsequently
		aligned with a number of other ESTs via the align() method.

        \return This method returns zero if the estIdx was within the
        given range of values.  Otherwise this method returns a
        non-zero value as the error code.
    */
    virtual int setReferenceEST(const EST* est) = 0;
    
    /** Main method to align a given EST with the reference sequence.

        This method provides the primary API into this analyzer for
        attempting to align a given EST with the reference sequence.
        Note that the reference sequence is set via a call to one of
        the setReferenceEST methods in this class.  Typically (this is
        not a stringent requirement), this method is implemented to
        determine the alignment (if relevant), in the following
        manner:

        <ol>

        <li>It first (via various helper methods) computes possible
        alignment using the normal (not reverse complement) nucleotide
        sequence for \c otherEST.</li>

        <li>If a successful alignment is not found in the previous
        step, then this method compute possible alignment using the
        reverse complementary nucleotide sequence for \c
        otherEST.</li>

        <li>If a valid alignment is found in either of the two steps,
        this method updates the supplied output parameters and returns
        \c true. If a valid alignment is not found, then the method
        returns \c false.</li>
        
        </ol>

        \param[in] otherEST A pointer to an immutable EST object with
		which the reference EST (or sequence) is to be analyzed and if
		valid/possible aligned.

        \param[out] info The information regarding the alignment (if
        one is computed) is populated into this object.  The alignment
        information includes the following information:

        <ul>

        <li>The index of the reference EST and the otherEST that were
        analyzed, found to be sufficiently similar, and aligned.</li>

        <li>The relative nucleotide index positions that indicates the
        centroid of the best possible alignment was found.</li>

        <li>The alignment score (for the best alignment) associated
        with the alignment between \c otherEST and the reference
        sequence.</li>

        <li>A normal or reverse complement (RC) flag to indicate if
        the alignment was accomplised with the regular or the reverse
        complement version of the \c otherEST.  This information is \c
        true if the alignment was accomplished with the reverse
        complement representation of the \c otherEST.  If the
        alignment is to be done with the normal sequence, then this
        parameter is set to \c false.

        </ul>
        
        \return This method returns \c true if a valid alignment was
        found between the reference sequence and \c otherEST.  If a
        valid alignment was not found (that means the pairs were
        determined to be unrelated), then this method returns \c
        false.
    */
    virtual bool align(const EST* otherEST, AlignmentInfo& info) = 0;

    /** Set the baton list cache to be used by this class.

        This method must be used to set the BatonListCache to be used
        by this object.  This method is typically used by an Assembler
        to set a shared cache for use by this class.

        \param[in] cache The cache to be used by this object.  This
        class saves and uses this pointer internally.  Consequently,
        the object must not be deleted.
    */
    void setBLCache(BatonListCache *cache);
    
protected:
    /** The constructor.
        
        This is the only constructor for this class. The constructor
        is not public to ensure that this class is not instantiated
        directly.  Instead one of the derived classes must be
        instantiated and used.  Currently the constructor does not
        have any special tasks to perform and merely initializes
        instance variables to their default initial value.

		\param[in] name A name to be associated with the derived
		class.  This passed to the base class that stores the name.
		The name is just a debugging/troubleshooting aid.
    */
    SequenceAligner(const std::string& name);

    /** The cache of baton lists to be used by this sequence aligner.

        This instance variable blCache (read as: baton list cache) is
        used used to maintain a cache of baton list objects that have
        been pre-computed via an on-demand approach.  The entries in
        this list consist of baton lists for both normal and reverse
        complement (RC) representation of a given cDNA fragment.  A
        valid pointer is set via the setBLCache() method.
    */
    BatonListCache* blCache;
    
private:
    /** Typically used to permit the component to perform its core
        tasks <b>but it is not used in a similar fashion here</b>.

        The SequenceAligner class hierarchy does not use the run()
        method as it is a helper class to be used with an Assembler.
        Consequently, the run() method is not meanigful and should not
        be used in this specific hierarchy.  Consequently, this method
        is used to override the default implementation in the base
        class and generate an error if it is ever invoked (even
        accidentally).
        
        \return This method returns zero on success.  If errors occur
        during initialization then this method must return an non-zero
        error code.  Currently this method does not return as it
        aborts execution after reporting an error.
    */
    int run();
};

#endif
