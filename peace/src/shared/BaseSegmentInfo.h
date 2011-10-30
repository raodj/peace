#ifndef BASE_SEGMENT_INFO_H
#define BASE_SEGMENT_INFO_H

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

/** \file BaseSegmentInfo.h

    \brief Class to encapsulate generic (that is -- genome-assembler
    neutral) base segment information associated a given Contig.  This
    information is primarily used to write assembly results to in
    appropriate file formats (such as: ACE).

    <p>A base segment is a sub-fragment of a cDNA sequence whose
	nucleotides have been primarily used to build a given region of a
	Contig. A contig can have multiple base segments.  For example
	consider the contig built from the following four fragments:</p>

	<table>
	  <tr><td>cDNA#0</td><td>\c   ATTCGATTAGATTAGATA            </td></tr>
	  <tr><td>cDNA#1</td><td>\c     TCGATTAGATT                 </td></tr>
	  <tr><td>cDNA#2</td><td>\c            GATTAGATATACTTACT    </td></tr>
	  <tr><td>cDNA#3</td><td>\c                      ACTTACTGGGA</td></tr>
	  <tr><td>Contig</td><td>\c\b ATTCGATTAGATTAGATATACTTACTGGGA</td></tr>
	</table>

	<p>In the above example, the Contig has the following three base
	segment regions:

	<table>
	  <th><td>B.Seg</td><td>cDNA</td><td>Start</td><td>End</td></th>
	  <tr><td>1</td><td>cDNA#0</td><td>1</td> <td>18</td></tr>
	  <tr><td>2</td><td>cDNA#2</td><td>19</td><td>26</td></tr>
	  <tr><td>3</td><td>cDNA#3</td><td>27</td><td>30</td></tr>	  
	</table>
	</p>
	
	This class encapsulates the essential information to identify a
	single base segment region.  This information is shared between
	multiple sub-systems constituting PEACE.  The assembly sub-system
	typically creates the data while the output sub-system consumes
	the information.  This class is primarily used to efficiently
	write alignemnt information to a file (in different formats such
	as: ACE). Currently, most of the alignment information in this
	file is geared towards writing ACE file formats (see
	ACEWriter).</p>

	<p>Various genomic-assemblers directly or indirectly utilize this
    class to hold alignment information for each contig. Typically,
    genomic-assemblers will require additional information (than what
    is encapsulated by this class) for their operations.  In such
    cases and when appropriate, this class can be used as a base
    class. Once a contig has been formed, the additional intermediate
    data can be dropped retaining only the core information in this
    class.  This minimizes peak memory footprint of assemblers
    (enabling assembly of larger data files).</p>

    <p>Refer to the class and method documentation for details on the
    data encapsulated by this class.</p>
*/

/** Class to encapsulate generic (genomic-assembler neutral) base
    segment information for a given contig.

    This class is primarily used to encapsulate alignment information
    pertaining to a given base segment.  Note that this class contains
    only the necessary information and not the nucleotide sequence
    associated with the actual fragment (the nucleotide sequence can
    be obtained via the various API methods associated with the EST
    class).  Consequently, this class is relatively straightforward.

    \see Contig
    \see ACEWriter
*/
class BaseSegmentInfo {
public:
    /** The primary default constructor that constitutes the main API
        for this class.

        The constructor merely initializes all the instance variables
        to invalid values.  Various setter and getter methods can then
        be used to suitably update the values in this object.  A
        default constructor is needed to ease the use of this class
        with various STL containers (such as: std::vector).
    */
    inline BaseSegmentInfo() : estIndex(-1), startPos(-1), endPos(-1) {}

    /** Convenience constructor to set most all the alignment
        information in one shot.

        \param[in] estIdx The index (zero based) of the cDNA for which
        a base segmeent object is being created.  This value must be
        in the range 0 \f$ \le \f$ estIdx &lt EST::getSequenceCount().

        \param[in] start The starting position on the contig for
        the base segment. This value is zero-based (that is the first
        position starts at zero). Consequently, this value cannot be
        negative.

        \param[in] end The ending position on the contig for the base
        segment. This value is zero-based (that is the first position
        starts at zero). Consequently, this value cannot be negative.
    */
    inline BaseSegmentInfo(const int estIdx,  const int start,
                           const int end) :
        estIndex(estIdx), startPos(start), endPos(end) {
        // Nothing else to be done for now in the constructor.
    }
    
    /** Convenience method to set the base segment information in one
        shot.

        \param[in] estIdx The index (zero based) of the cDNA for which
        an alignment information object is being created.  This value
        must be in the range 0 \f$ \le \f$ estIdx &lt
        EST::getSequenceCount().

        \param[in] start The starting position on the contig for
        the base segment. This value is zero-based (that is the first
        position starts at zero). Consequently, this value cannot be
        negative.

        \param[in] end The ending position on the contig for the base
        segment. This value is zero-based (that is the first position
        starts at zero). Consequently, this value cannot be negative.
    */
    void setInfo(const int estIdx,  const int start,
                 const int end);

    /** Obtain the index (or ID) of the fragment that contributes its
        sub-fragment to the contig.

        This method can be used to obtain the index of the cDNA entry
        which constitutes the base segment for a given region of the
        contig.  This value can be used to obtain the actual
        nucleotide sequence or its information via a call to
        EST::getEST(ai.getESTIndex()).  This method essentially
        returns the value that was set when this object was created.

        \return This method returns the index into the list of
        fragments maintained by the EST class.  For valid alignment
        objects, this value is in the range 0 \f$ \le \f$ \e retVal
        &lt EST::getSequenceCount().
    */
    inline int getESTIndex() const { return estIndex; }

    /** Obtain the relative starting position on the contig for this
        base segment entry.

        This method returns the starting position for the first
        nucleotide position on the base segment. Consequently, this
        value cannot be negative.

        \return The non-negative, zero-based starting position of the
        first nucleotide that contributes to the base segment.
    */
    inline int getStartPos() const { return startPos; }

    /** Obtain the relative ending position on the contig for this
        base segment entry.

        This method returns the ending position for the last
        nucleotide position on the base segment. Consequently, this
        value cannot be negative.

        \return The non-negative, zero-based starting position of the
        last nucleotide in this base segment that contributes to the
        contig.
    */
    inline int getEndPos() const { return endPos; }

protected:
    /** The index (or ID) of the fragment whose base segment
        information is contained in this object.
        
        This instance variable tracks the index of the fragment for
        which this object contains the necessary alignment data.  This
        value can be used to obtain the actual nucleotide sequence via
        a call to EST::getEST(ai.getESTIndex()).  Typically, this
        value is set via the setInfo() method in this class.
    */
    int estIndex;

    /** The relative starting position on the contig for this base
        segment.

        This value is zero-based (that is the first position starts at
        zero) position of the <em>first</em> nucleotide on the
        contig. Consequently, this value cannot be negative.
    */
    int startPos;

    /** The relative ending position on the contig for this base
        segment.

        This value is zero-based (that is the first position starts at
        zero) position of the <em>last</em> nucleotide on the
        contig. Consequently, this value cannot be negative.
    */
    int endPos;

private:
    // Currently this class does not have any private members.
};

#endif
