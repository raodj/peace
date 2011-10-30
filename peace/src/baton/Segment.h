#ifndef SEGMENT_H
#define SEGMENT_H

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

#include <ostream>
#include <vector>

/** \file Segment.h

    \brief Class to encapsulate sub-alignment information within an
    BatonAlignmentInfo object.

    <p>This class essentially contains a sub-alignment between a given
    reference-cDNA fragment and another cDNA fragment.  A
    sub-alignment information essentially identifies sub-fragments
    from two given cDNA fragments that are sufficiently identical to
    each other indicating that they were formed from the same
    gene. This class is not meant to be directly instantiated and
    used. Instead the BatonAlignmentInfo class provides the primary
    interface for this class.</p> */

/** Class to encapsulate sub-alignment information.

	<p>This class essentially contains a sub-alignment between a given
    reference-cDNA fragment and another cDNA fragment.  A
    sub-alignment information essentially identifies sub-fragments
    from two given cDNA fragments that are sufficiently identical to
    each other indicating that they were formed from the same gene.</p>

    <p>For example consider the following two cDNA fragments (assume
    that the first one is the reference EST and the second one is the
    \i other EST). The nucleotides in capital letters are
    intentionally used to highlight two regions of alignment between
    the two fragments. These two regsions of alignment would be
    represented using two different Segment objects (in one
    BatonAlignmentInfo object) similar to:

    <ol>

    <li>Segment(9, 5, 7)</li>

    <li>Segment(19, 15, 6)</li>

    </ol>

    </p>

    \code

    atttcatagTATTTCGataGGATGGatta
        gtaagTATTTCGcatGGATGGcccg

    \endcode

	<p>This class is primarily used to encapsulate a part of the
    alignment information pertaining to a given pair of cDNA
    fragments.  Furthermore this class does not contain the nucleotide
    sequence associated with the actual fragment (the nucleotide
    sequence can be obtained via the various API methods associated
    with the EST class).  Consequently, this class is relatively
    straightforward.</p>

	<p>This class is not really meant to be used as a standalone
    class.  Instead it is meant to be indirectly used via the
    BatonAlignmentInfo class. Objects of this type are typically
    created by Sequence Aligner class and added to the
    BatonAlignmentInfo class for further use.</p>
	
    \note This object is sent/received over the wire via
    MPI_SEND/MPI_RECV calls. Consequently, this class cannot have
    pointers!

    \see BatonAlignmentInfo
*/
class Segment {
    // Insertion operator to make dumping Segments for debugging easier.
    friend std::ostream& operator<<(std::ostream&, const Segment&);    
public:
    /** The primary constructor for creating a sub-alignment
        segment.

        <p>This constructor is primarily used by the
        BatonAlignmentInfo class to create a Segment object whenever
        new sub-alignments are added.  This constructor is relatively
        straightforward in that it merely initializes the instance
        variables with the values supplied as parameters.</p>
        
        <p>Ideally the instance variables in this class would all be
        constant because this class is meant to contain immutable data
        --- that is, once the data set we never want to change them.
        However, objects containing constant instance variables, don't
        play well with many of the STL container classes.
        Consequently, this class is designed to permit setting values
        via the constructor and once an object is created, its values
        cannot be changed as there are (intentionally) no setter
        methods in this class.</p>

        \param[in] refOffset The offset (zero based) into the
        reference EST where this segment of alignment begins in the
        reference EST.  The actual index of the reference EST is
        present in the BatonAlignmentInfo class that contains this
        Segment.

        \param[in] otherOffset The offset (zero based) into the other
        EST (which is being aligned with the reference EST) from where
        this segment of alignment begins.  The actual index of the
        other EST is present in the BatonAlignmentInfo class that
        contains this Segment.

        \param[in] len The length (in number of nucleotides) that are
        aligned as a part of this segment. This value remains the same
        for both the reference and other EST.
    */
    inline Segment(const int refOffset, const int otherOffset, const int len) :
        refESTOffset(refOffset), othESTOffset(otherOffset),
        alignmentLen(len) {}
    
	/** Copy constructor.

		This is a convenience copy constructor that is used to copy
		the information from a given source alignment information
		object to a newly created one.

		\param[in] src The source alignment information object from
		where the data for this object is to be copied.
	 */
    inline Segment(const Segment& src) :
        refESTOffset(src.refESTOffset), othESTOffset(src.othESTOffset),
        alignmentLen(src.alignmentLen) {}
    
    /** A default constructor.

        The default constructor essentially initializes all the
        instance variables to some initial invalid values. The default
        constructor is needed so that this class can be conveniently
        used with various STL containers (such as: std::vector).
    */
    Segment() : refESTOffset(-1), othESTOffset(-1), alignmentLen(-1) {}

    /** Obtain the offset of the first nucleotide in a cDNA fragment
        from where this sub-alignment commences.
        
        This method can be used to obtain the zero-based offset of the
        first nucleotide in the cDNA fragment from where this
        sub-alignment commences.  This method essentially returns the
        value that was set when this object was created.

        \return This method returns the zero-based offset into the
        nucleotide sequence (stored as a string).
    */
    inline int getOthESTOffset() const { return othESTOffset; }

    /** Obtain the offset of the first nucleotide in the reference
        cDNA fragment from where this sub-alignment commences.
        
        This method can be used to obtain the zero-based offset of the
        first nucleotide in the reference fragment from where this
        sub-alignment commences.  This method essentially returns the
        value that was set when this object was created.

        \return This method returns the zero-based offset into the
        reference nucleotide sequence (stored as a string).
    */
    inline int getRefESTOffset() const { return refESTOffset; }
    
    /** Obtain the number of nucleotides constituting this
        sub-alignment.

        This method can be used to determine the number of nucleotides
        constituting this segment of alignment.  The length of the
        sub-alignment is the same for both the cDNA fragments being
        represented by this Segment.

        \return The number of nucleotides constituting this
        sub-alignment.
    */
    inline int getLength() const { return alignmentLen; }

	/** Obtain offset the last nucleotide in the reference cDNA
		fragment constituting this sub-alignment

		This is a convenience method that can be used to obtain the
		zero-based offset of the last nculeotide in the reference cDNA
		fragment where this sub-alignment ends. This method
		essentially returns the value of the expression:
		getRefESTOffset() + getLength() - 1.

		\return The zero-based index/offset (into the cDNA string) of
		the last nucleotide constituting this sub-alignment.
	*/
	inline int getRefESTEndOffset() const
	{ return refESTOffset + alignmentLen - 1; }

	/** Obtain offset the last nucleotide in the other cDNA fragment
		constituting this sub-alignment

		This is a convenience method that can be used to obtain the
		zero-based offset of the last nculeotide in the other cDNA
		fragment where this sub-alignment ends. This method
		essentially returns the value of the expression:
		getOthESTOffset() + getLength() - 1.

		\return The zero-based index/offset (into the cDNA string) of
		the last nucleotide constituting this sub-alignment.
	*/
	inline int getOthESTEndOffset() const
	{ return othESTOffset + alignmentLen - 1; }

	/** Obtain \i relative offset difference between the first
		nucleotide of the reference and other cDNA fragment
		constituting this sub-alignment.
		
		This is a convenience method that can be used to obtain the
		relative offset difference between first nculeotide of the
		reference and other cDNA fragment. This method essentially
		returns the value of the expression: getESTOffset() -
		getOthESTOffset().

		\return The difference between the index/offset (into the cDNA
		string) of the first last nucleotide of the reference and
		other cDNA fragment. The returned value can be negative, zero,
		or positive.
	*/
	inline int getRelativeOffset() const
	{ return refESTOffset - othESTOffset; }

    /** Convenience method to update the offsets by a given detla value.

        This method is a convenience method that is used to shift the
        offsets for this segment by a given constant value.  Shifting
        segments is used when sub-contigs are merged and the position
        of this Segment in the final contig needs to be updated

        \param[in] delta The constant delta value to be added to
        <i>only the reference EST offset values</i> stored in this
        class. The otherESTOffset remains unchanged by this method.
    */
	inline void adjustRefOffset(const int delta) {
        refESTOffset += delta;
	}
	
protected:
    // Currently this class does not have any protected members
        
private:
    /** The zero-based offset of the first nucleotide in the reference
        cDNA fragment that is part of this sub-alignment.
        
        This instance variable tracks the zero-based offset of the
        first nucleotide in the reference cDNA fragment (with which
        another cDNA is being aligned).  This value is in the range 0
        (zero) to the length of the reference cDNA fragment.  This
        value is set in the constructor and is never changed during
        the life time of this class.
    */
    int refESTOffset;

    /** The zero-based offset of the first nucleotide in the other
        cDNA fragment that is part of this sub-alignment.
        
        This instance variable tracks the zero-based offset of the
        first nucleotide in the other cDNA fragment (which is being
        aligned with the reference cDNA).  This value is in the range
        0 (zero) to the length of the other cDNA fragment.  This value
        is set in the constructor and is never changed during the life
        time of this class.
    */    
    int othESTOffset;
    
    /** The number of nucleotides constituting the segment of
        alignment.

        This instance variable is used to maintain the length (in
        number of nucleotides) that constitute this segment of
        sub-alignment.
    */
    int alignmentLen;
};


/** \typedef std::vector<Segment> SegmentList

    \brief A vector that contains a list of SegmentList objects.

    This typedef provides a convenient short cut to refer to a vector
    that is used to hold a list of Segment objects.  This vector is
    used in more than one spot in the code. Consequently, it has been
    defined here to serve as a common data structure.
*/
typedef std::vector<Segment> SegmentList;

/** \fn std::ostream& operator<<(std::ostream&, const Segment&)

    Insertion operator to stream Segment information to a given output
    stream.  This method provides a convenient mechanism to dump the
    complete Segment information for debugging purposes.

    \param[out] os The output stream to which the segment information
    is to be written.

    \param[in] segment The segment whose information is to be
    displayed.
    
    \return As the generalized API contract for insertion operators,
    this method returns os (first parameter)
*/
extern std::ostream& operator<<(std::ostream& os, const Segment& baton);

#endif
