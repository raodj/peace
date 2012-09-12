#ifndef ALIGNMENT_INFO_H
#define ALIGNMENT_INFO_H

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

#include <string>

/** \file AlignmentInfo.h

    \brief Class to encapsulate generic (that is -- genome-assembler
    neutral) alignment information about a specific cDNA
    fragment. This information is primarily used to write assembly
    results to in appropriate file formats (such as: FASTA or SAM).

    <p>This class contains only the essential alignment information
    about one single cDNA fragment.  This information is shared
    between multiple sub-systems constituting PEACE.  The assembly
    sub-system typically creates the data while the output sub-system
    consumes the information.  This class is primarily used to
    efficiently write alignemnt information to a file (in different
    formats such as: FASTA or SAM). Currently, most of the alignment
    information in this file is geared towards writing SAM file
    formats (see SAMWriter).</p>

	<p>Various genomic-assemblers directly or indirectly utilize this
    class to hold alignment information for each cDNA
    fragment. Typically, genomic-assemblers will require additional
    information (than what is encapsulated by this class) for their
    operations.  In such cases and when appropriate, this class can be
    used as a base class. Once a contig has been formed, the
    additional intermediate data can be dropped retaining only the
    core information in this class.  This minimizes peak memory
    footprint of assemblers (enabling assembly of larger data
    files).</p>

    <p>Refer to the class and method documentation for details on the
    data encapsulated by this class.</p>
*/

/** Class to encapsulate generic (genomic-assembler neutral) alignment
	information for a single cDNA fragment.

    This class is primarily used to encapsulate alignment information
    pertaining to a given cDNA fragments.  Note that this class
    contains only the alignment information and not the nucleotide
    sequence associated with the actual fragment (the nucleotide
    sequence can be obtained via the various API methods associated
    with the EST class).  Consequently, this class is relatively
    straightforward.

    \note This object is <b>not meant to be \i directly </b>
    sent/received over the wire via MPI_SEND/MPI_RECV calls (as it
    contains string that are meant to be used on a local machine).

    \see Contig
    \see SAMWriter
*/
class AlignmentInfo {
public:
    /** The primary default constructor that constitutes the main API
        for this class.

        The constructor merely initializes all the instance variables
        to invalid values.  Various setter and getter methods can then
        be used to suitably update the values in this object.  A
        default constructor is needed to ease the use of this class
        with various STL containers (such as: std::vector).
    */
    inline AlignmentInfo() : estIndex(-1), rcAlignment(false), position(-1),
                             mapQuality(255), cigar(""), templateLen(-1),
                             score(-1) {}

    /** Convenience constructor to set most all the alignment
        information in one shot.

        \param[in] estIdx The index (zero based) of the cDNA for which
        an alignment information object is being created.  This value
        must be in the range 0 \f$ \le \f$ estIdx &lt
        EST::getSequenceCount().

        \param[in] revCompFlag This parameter is set to \c true if
        the alignment was accomplished with the reverse complement
        representation of the estIdx.  If the alignment is to be done
        with the normal sequence, then this parameter is set to \c
        false.

        \param[in] position The relative starting position on the
        contig for the cDNA fragment. This value is zero-based (that
        is the first position starts at zero). Consequently, this
        value cannot be negative.

        \param[in] cigarString The CIGAR string that contains the
        concise representation of the alignment between a given cDNA
        fragment (estIdx) and the contig to which this fragment
        belongs.  The CIGAR string representation used by SAM file
        format (see: http://samtools.sourceforge.net/SAM1.pdf) is used
        as the standard here.  Here are a few examples of CIGAR
        strings: \c 8M2I4M1D3M, \c 3S6M1P1I4M, \c 9M.
        
        \param[in] mapQuality The mapping quality value for the
        alignment associated with the given cDNA fragment.  It equals
        <em>-10 log<sub>10</sub>Pr{mapping position is wrong}</em>,
        rounded to the nearest integer. A value 255 indicates that the
        mapping quality is not available.

        \param[in] templateLen The signed observed Template LENgth for
        a \b read-pair. If all segments are mapped to the same
        reference, the unsigned observed template length equals the
        number of bases from the leftmost mapped base to the rightmost
        mapped base. The leftmost segment has a plus sign and the
        rightmost has a minus sign. The sign of segments in the middle
        is undefined. It is set as 0 for single-segment template or
        when read-pair information is unavailable.
        
        \param[in] alignScore An optional alignment score (if any) for
        this alignment.  The alignment score is completely dependent
        on the genomic-assembler being used as different approaches
        may be used to generate the score.
    */
    inline AlignmentInfo(const int estIdx,  const bool revCompFlag,
                         const int pos, const std::string& cigarStr,
                         const int mapQual, const int templLen,
                         const int alignScore) :
        estIndex(estIdx), rcAlignment(revCompFlag), position(pos),
        mapQuality(mapQual), cigar(cigarStr), templateLen(templLen),
        score(alignScore) {
        // Nothing else to be done for now in the constructor.
    }
    
    /** Convenience method to set the alignment information in one
        shot.

        \param[in] estIdx The index (zero based) of the cDNA for which
        an alignment information object is being created.  This value
        must be in the range 0 \f$ \le \f$ estIdx &lt
        EST::getSequenceCount().

        \param[in] revCompFlag This parameter is set to \c true if
        the alignment was accomplished with the reverse complement
        representation of the estIdx.  If the alignment is to be done
        with the normal sequence, then this parameter is set to \c
        false.

        \param[in] position The relative starting position on the
        contig for the cDNA fragment. This value is zero-based (that
        is the first position starts at zero). Consequently, this
        value cannot be negative.

        \param[in] cigarString The CIGAR string that contains the
        concise representation of the alignment between a given cDNA
        fragment (estIdx) and the contig to which this fragment
        belongs.  The CIGAR string representation used by SAM file
        format (see: http://samtools.sourceforge.net/SAM1.pdf) is used
        as the standard here.  Here are a few examples of CIGAR
        strings: \c 8M2I4M1D3M, \c 3S6M1P1I4M, \c 9M.
        
        \param[in] mapQuality The mapping quality value for the
        alignment associated with the given cDNA fragment.  It equals
        <em>-10 log<sub>10</sub>Pr{mapping position is wrong}</em>,
        rounded to the nearest integer. A value 255 indicates that the
        mapping quality is not available.

        \param[in] templateLen The signed observed Template LENgth for
        a \b read-pair. If all segments are mapped to the same
        reference, the unsigned observed template length equals the
        number of bases from the leftmost mapped base to the rightmost
        mapped base. The leftmost segment has a plus sign and the
        rightmost has a minus sign. The sign of segments in the middle
        is undened. It is set as 0 for single-segment template or when
        read-pair information is unavailable.
        
        \param[in] alignScore An optional alignment score (if any) for
        this alignment.  The alignment score is completely dependent
        on the genomic-assembler being used as different approaches
        may be used to generate the score.
    */
    void setInfo(const int estIdx,  const bool revCompFlag,
                 const int position, const std::string& cigar,
                 const int mapQuality = 255, const int templateLen = 0,
                 const int alignScore = -1);

    /** Obtain the index (or ID) of the fragment for which this object
        contains alignment information.

        This method can be used to obtain the index of the fragment
        for which this object contains the necessary alignment data.
        This value can be used to obtain the actual nucleotide
        sequence via a call to ESTList::getEST(ai.getESTIndex()).  This
        method essentially returns the value that was set when this
        object was created.

        \return This method returns the index into the list of
        fragments maintained by the EST class.  For valid alignment
        objects, this value is in the range 0 \f$ \le \f$ \e retVal
        &lt EST::getSequenceCount().
    */
    inline int getESTIndex() const { return estIndex; }
    
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

    /** Obtain the relative starting position on the contig for the
        cDNA fragment.

        This method returns the relative starting position for the
        first aligned nucleotide for this cDNA fragment.  This
        position does not include any left-most mismatched or clipped
        nucleotide.  This value is zero-based (that is the first
        position starts at zero) position on the contig. Consequently,
        this value cannot be negative.

        \return The non-negative, zero-based starting position of the
        first aligned nucleotide.
    */
    inline int getPosition() const { return position; }

    /** Obtain the mapping quality value for the alignment associated
        with the given cDNA fragment.

        The mapping qulity provides the overall quality of the
        mapping.  It equals <em>-10 log<sub>10</sub>Pr{mapping
        position is wrong}</em>, rounded to the nearest integer. A
        value of 255 (which is the default value) indicates that the
        mapping quality is not available.

        \return The overall mapping quality value set for this
        alignment.
    */    
    inline int getMapQuality() const { return mapQuality; }
    
    /** Obtain the CIGAR string that provides a concise representation
        of the alignment for the given cDNA fragment.

        This method returns the CIGAR string associated with the
        alignment information.  The alignment information is with
        respect to the contig to which this fragment belongs.  The
        CIGAR string representation used by SAM file format (see:
        http://samtools.sourceforge.net/SAM1.pdf) is used as the
        standard here.  Here are a few examples of CIGAR strings: \c
        8M2I4M1D3M, \c 3S6M1P1I4M, \c 9M.

        \return The CIGAR string for this alignment information. If a
        valid CIGAR string is not set then this method returns an
        empty string (\c "").
    */
    inline const std::string& getCigar() const { return cigar; }

    /** Obtain the signed observed Template LENgth for a \b read-pair
        alignment.

        This value is meangiful only when this alignment corresponds
        to a read-pair (and assuming the genomic-assembler explicitly
        detected and uses read-pairs for assembly).  In all other
        cases this value is zero.

        \return This method returns the signed observed template
        length. The default return value is zero.
    */    
    inline int getTemplateLen() const { return templateLen; }
    
    /** Obtain the score or strength of the alignment between this
        fragment and its reference.

        This method can be used to an optional determine the
        alignment-score for this fragment (whose index is returned by
        getESTIndex() method). Note that alignment scores are
        optional. The default score is -1.

        \return The score between the given fragment and its reference
        fragment.
    */
    inline int getScore() const { return score; }

protected:
    /** The index (or ID) of the fragment whose alignment information
        is contained in this object.
        
        This instance variable tracks the index of the fragment for
        which this object contains the necessary alignment data.  This
        value can be used to obtain the actual nucleotide sequence via
        a call to EST::getEST(ai.getESTIndex()).  Typically, this
        value is set via the setInfo() method in this class.
    */
    int estIndex;

    /** Flag to indicate of regular or reverse-complement sequence
        must be used for building consensus sequence.

        This flag contains \c true if the resulting consensus sequence
        must use the reverse complementary representation of the
        fragment obtained via a call to EST::getEST(ai.getESTIndex()).
        If the regular nucleotide sequence is to be used, then this
        flag is set to \c false.
    */
    bool rcAlignment;

    /** The relative starting position on the contig for the cDNA
        fragment.

        This value is zero-based (that is the first position starts at
        zero) position of the <em>first aligned</em>
        nucleotide. Consequently, this value cannot be negative.  This
        value corresponds to the \c POS field in a SAM file.
    */
    int position;

    /** The mapping quality value for the alignment associated with
        the given cDNA fragment.

        The mapping qulity provides the overall quality of the
        mapping.  It equals <em>-10 log<sub>10</sub>Pr{mapping
        position is wrong}</em>, rounded to the nearest integer. A
        value of 255 (which is the default value) indicates that the
        mapping quality is not available.
    */
    int mapQuality;

    /** The CIGAR string that provides a concise representation of the
        alignment for the given cDNA fragment.

        This CIGAR string associated with this alignment information.
        The alignment information is with respect to the contig to
        which this fragment belongs.  The CIGAR string representation
        used by SAM file format (see:
        http://samtools.sourceforge.net/SAM1.pdf) is used as the
        standard here.  Here are a few examples of CIGAR strings: \c
        8M2I4M1D3M, \c 3S6M1P1I4M, \c 9M.  The default value is empty
        string (\c "").
    */
    std::string cigar;
	
    /** The signed observed Template LENgth for a \b read-pair.

        This value is meangiful only when read-pairs are explicitly
        detected and used for assembly.  In all other cases this value
        must be set to zero.  In the case of explicitly processing
        read-pairs, the following applies:

        <dl>
        <dt></dt>
        <dd>If all segments are mapped to the same reference, the
        unsigned observed template length equals the number of bases
        from the leftmost mapped base to the rightmost mapped
        base. The leftmost segment has a plus sign and the rightmost
        has a minus sign. The sign of segments in the middle is
        undened. It is set as 0 for single-segment template or when
        read-pair information is unavailable.</dd>
        </dl>
    */
    int templateLen;
    
    /** The score or strength of the alignment between this fragment
        and its reference

        This instance variable holds the score between the fragment
        (whose index is returned by getESTIndex() method) and the
        reference fragment (whose index is returned by the
        getRefESTIndex() method) in this alignment info. This value is
        typically the sum of the length of the Segment objects in
        segmentList.
    */
    int score;

private:
    // Currently this class does not have any private members.
};

#endif
