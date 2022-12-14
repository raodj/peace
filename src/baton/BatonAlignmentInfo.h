#ifndef BATON_ALIGNMENT_INFO_H
#define BATON_ALIGNMENT_INFO_H

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

#include "Segment.h"

/** \file BatonAlignmentInfo.h

    \brief Class to encapsulate alignment information that is
    generated, maintained, and used by various classes constituting
    the baton analyzer.

    This class essentially contains the set of information that is
    mainained for the assembly generated by the baton assembler.
    Refer to the class and method documentation for details on the
    data encapsulated by this class.  This file also provides a
    convenience type definition for the BatonAlignmentInfoList data type.
*/

/** Class to encapsulate alignment information.

    This class is primarily used to encapsulate alignment information
    pertaining to a given pair of cDNA fragments.  Note that this
    class contains only the alignment information and not the
    nucleotide sequence associated with the actual fragment (the
    nucleotide sequence can be obtained via the various API methods
    associated with the EST class).  Consequently, this class is
    relatively straightforward.  This class is used by several other
    classes constituting the baton assembler.  Therefore, it has been
    refactored into its own independent class.

    \note This object is <b>not meant to be</b> sent/received over the
    wire via MPI_SEND/MPI_RECV calls (as it contains vectors that are
    meant to be used on a local machine).

    \see BatonAnalyzer
    \see ConsensusMaker
    \see BatonAssembler
*/
class BatonAlignmentInfo {
public:
    /** The primary default constructor that constitutes the main API
        for this class.

        The constructor merely initializes all the instance variables
        to invalid values.  Various setter and getter methods can then
        be used to suitably update the values in this object.
    */
    BatonAlignmentInfo() : estIndex(-1), refESTindex(-1), score(-1),
                      rcAlignment(false) {}


    /** Convenience method to set most (but not all) the alignment
        information in one shot.

        \param[in] refESTIdx The index (zero based) of the reference
        EST with which the fragment at estIdx index position has been
        aligned.  for which an alignment information object is being
        created.  This value must be in the range 0 \f$ \le \f$
        estIdx &lt EST::getSequenceCount().
        
        \param[in] othESTIdx The index (zero based) of the EST for
        which an alignment information (with-respect-to the reference
        cDNA) object is being created.  This value must be in the
        range 0 \f$ \le \f$ estIdx &lt EST::getSequenceCount().

        \param[in] revCompFlag This parameter is set to \c true if
        the alignment was accomplished with the reverse complement
        representation of the estIdx.  If the alignment is to be done
        with the normal sequence, then this parameter is set to \c
        false.
        
        \param[in] alignScore The alignment score (for the best
        alignment) associated with the alignment between \c estIdx and
        the reference sequence at \c refESTidx respectively. This
        value is typically the sum of the lengths of the Segments
        constituting the alignment.
    */
    void setInfo(const int estIdx, const int refESTIdx,
                 const int alignScore, const bool revCompFlag);
    
    /** Obtain the index (or ID) of the fragment for which this object
        contains alignment information.

        This method can be used to obtain the index of the fragment
        for which this object contains the necessary alignment data.
        This value can be used to obtain the actual nucleotide
        sequence via a call to EST::getEST(ai.getESTIndex()).  This
        method essentially returns the value that was set when this
        object was created.

        \return This method returns the index into the list of
        fragments maintained by the EST class.  For valid alignment
        objects, this value is in the range 0 \f$ \le \f$ \e retVal
        &lt EST::getSequenceCount().
    */
    inline int getESTIndex() const { return estIndex; }

    /** Obtain the index of the parent or reference fragment with
        which the alignment was determined.

        The alignment information for fragments is determined with
        respect to another reference fragment.  For the first (or
        root) fragment, this value is -1.  All other fragments have a
        valid reference fragment index.  This value can be used to
        obtain the actual reference nucleotide sequence via a call to
        EST::getEST(ai.getRefESTIndex()).  This method essentially
        returns the value that was set when this object was created.

        \return This method returns the index into the list of
        fragments maintained by the EST class.  For valid alignment
        objects, this value is in the range 0 \f$ \le \f$ \e retVal
        &lt EST::getSequenceCount().
    */
    inline int getRefESTIndex() const { return refESTindex; }

    /** Obtain a mutable reference to the vector that contains the
        various Segment object that contain sub-alignment information.

        This method can be used to obtain the list of Segment objects
        that constitute the sub-alignment regions. Each Segment
        delineates the sub-fragments that have batching
        nucleotides. These regions of alignment are identified using
        Batons.

        \note The Segment objects are sorted in a logical
        left-to-right manner.
        
        \return List of Segment objects that identify regions of
        alignment between the reference and other cDNA sequences.
    */
    inline SegmentList& getSegmentList() { return segmentList; }

    /** Obtain a immutable (read-only) reference to the vector that
        contains the various Segment object that contain sub-alignment
        information.

        This method can be used to obtain the list of Segment objects
        that constitute the sub-alignment regions. Each Segment
        delineates the sub-fragments that have batching
        nucleotides. These regions of alignment are identified using
        Batons.

        \note The Segment objects are sorted in a logical
        left-to-right manner.
        
        \return List of Segment objects that identify regions of
        alignment between the reference and other cDNA sequences. The
        returned list cannot be modified.
    */
    inline const SegmentList& getSegmentList() const { return segmentList; }
    
    /** Obtain the score or strength of the alignment between this
        fragment and its reference.

        This method can be used to determine the score between the
        fragment (whose index is returned by getESTIndex() method) and
        the reference fragment (whose index is returned by the
        getRefESTIndex() method) in this alignment info.

        \return The score between the given fragment and its reference
        fragment.
    */
    inline int getScore() const { return score; }

    /** Determine if the fragment is best aligned via its reverse
        complementary sequence.

        This method can be used to determine if the regular or the
        reverse-complementary representation of the nucleotide
        sequence associated with fragment (whose index is returned by
        getESTIndex() method) must be used for building the consensus
        sequence.  For example, if the nucleotide sequence is \c
        AATGCC and this method returns \c false, then the reverse
        complementary sequence \c GGCATT must be used to construct the
        final consensus sequence.
        
        \return This method returns \c true if the resulting consensus
        sequence must use the reverse complementary representation of
        the fragment.  Otherwise this method returns \c false.
     */
    inline bool isRCAlignment() const { return rcAlignment; }

    /** Convenience method to update the offsets in the segments by a
		given detla value.

        This method is a convenience method that is used to shift the
        offsets for all the segments in this class by a given constant
        value.  Shifting segments is used when sub-contigs are merged
        and the position of this Segment in the final contig needs to
        be updated. This method essentially calls
        Segment::adjustOffsets() with the delta value on all objects
        in the segmentList.

        \param[in] delta The constant delta value to be added to \i
        only the reference EST offset values stored in all objects in
        the segmentList associated this object.
    */
	void adjustRefOffsets(const int delta);

    /** Generates a formatted contents of this alignment for
        convenient analysis.

        This method provides a textual view of the alignment
        information for a given cDNA fragment.  The display is handy
        for analyzing alignments and troubleshooting any issues.
        
        \param[out] os The output stream to which the contents of this
        contig are to be written.

		\param[in] seq The nucleotide sequence of the cDNA associated
		with this alignment.  This sequence should be the same
		sequence associated with EST with index estIndex.
		
		\param[in] refESTpos The relative starting position of the
		reference fragment in the contig to which this object
		belongs. This information is used to appropriately
		(horizontally) place the fragment in the contig data being
		printed.

		@see ContigMaker::prettyPrint()
    */
    void prettyPrint(std::ostream& os, const std::string& seq,
					 const int refESTpos) const;

    /** Helper method to get a SAM CIGAR string for the segments
        represented by this alignment.

        The CIGAR string that provides a concise representation of the
        alignment for the given cDNA fragment.  The alignment
        information is with respect to the contig to which this
        fragment belongs.  The CIGAR string representation used by SAM
        file format (see: http://samtools.sourceforge.net/SAM1.pdf) is
        used as the standard here.  Here are a few examples of CIGAR
        strings: \c 8M2I4M1D3M, \c 3S6M1P1I4M, \c 9M.  The default
        value is empty string (\c "").

		\param[in] seq The nucleotide sequence of the cDNA associated
		with this alignment.  This sequence should be the same
		sequence associated with EST with index estIndex.

        \return The SAM CIGAR string that provides a concise summary
        of the alignment.
    */
    std::string getCIGAR(const std::string& seq) const;

    /** Obtain the starting offset in the reference contig for this
        alignment.

        This is a simple convenience method that can be used to obtain
        the starting offset in the contig, for the first-aligned
        nucleotide.  This method essentially returns the reference EST
        offset stored in the first aligned segment (if one is
        available) as shown below:

        \code

        return (segmentList.size() > 0) ?
        segmentList[0].getRefESTOffset() : -1;
            
        \endcode

        \return The starting offset in the contig, for the
        first-aligned nucleotide.  If no alignment is present, then
        this method returns -1.
    */
    inline int getStartOffset() const {
        return (segmentList.size() > 0) ?
            segmentList[0].getRefESTOffset() : -1;
    }
    
protected:
    // Currently this class does not have any protected members
    
private:
    /** The index (or ID) of the fragment whose alignment information
        is contained in this object.
        
        This instance variable tracks the index of the fragment for
        which this object contains the necessary alignment data.  This
        value can be used to obtain the actual nucleotide sequence via
        a call to EST::getEST(ai.getESTIndex()).  Typically, this
        value is set via the setInfo() method in this class.
    */
    int estIndex;

    /** The reference fragment with respect to which the alignment
        information was computed.

        The alignment information for fragments is determined with
        respect to another reference fragment.  This instance variable
        maintains the index of the reference cDNA fragment. For an
        artificially generated seed cDNA (that is obtained from an
        intermediate consensus sequence), this value is -1.  All other
        fragments have a valid reference fragment index.  This value
        can be used to obtain the actual reference nucleotide sequence
        via a call to EST::getEST(ai.getRefESTIndex()).  This method
        essentially returns the value that was set via the setInfo
        method in this class.
    */
    int refESTindex;

    /** The list of sub-alignment Segment objects that constitute the
        overall alignment.

        This instance variable is used to maintain the list of Segment
        objects that constitute the alignment between two given cDNA
        fragments.  This object is currently populated by
        DefaultSequenceAligner. Maybe these operations could be
        migrated into this class.
    */
    SegmentList segmentList;

    /** Tthe score or strength of the alignment between this fragment
        and its reference

        This instance variable holds the score between the fragment
        (whose index is returned by getESTIndex() method) and the
        reference fragment (whose index is returned by the
        getRefESTIndex() method) in this alignment info. This value is
        typically the sum of the length of the Segment objects in
        segmentList.
    */
    int score;

    /** Flag to indicate of regular or reverse-complement sequence
        must be used for building consensus sequence.

        This flag contains \c true if the resulting consensus sequence
        must use the reverse complementary representation of the
        fragment obtained via a call to EST::getEST(ai.getESTIndex()).
        If the regular nucleotide sequence is to be used, then this
        flag is set to \c false.
    */
    bool rcAlignment;
};


/** \typedef std::vector<BatonAlignmentInfo> BatonAlignmentInfoList

    \brief A vector that contains a list of BatonAlignmentInfo objects.

    This typedef provides a convenient short cut to refer to a vector
    that is used to hold a list of BatonAlignmentInfo objects.  This
    vector is used in more than one spot in the code. Consequently, it
    has been defined here to serve as a common data structure.
*/
typedef std::vector<BatonAlignmentInfo> BatonAlignmentInfoList;

#endif
