#ifndef CONTING_MAKER_H
#define CONTIG_MAKER_H

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

#include "AlignmentInfo.h"
#include <vector>

/** \file ContigMaker.h

    \brief Class to encapsulate a set of related alignment information
    and finally construct a contig.

    This class essentially contains the set of related fragments that
    are finally combined together to form a contig. This class is
    primarily used by the BatonAssemblerManager to build contigs.
    Refer to the class and method documentation for details on this
    class.
*/

/** \typedef std::vector< std::vector<int> > TwoDArray

    \brief A vector-of-vector-of-integers that serves as a classical
    two-dimensional array.

    This typedef provides a convenient short cut to refer to a
    vector-of-vectors that is used to hold a list of integers.
    Currently this list is used to build the base distributions, \i
    i.e., the number of times a given nucleotide occurs at a given
    index position.  Using this construct facilitates the creation and
    management of 2-D arrays.
*/
typedef std::vector< std::vector<int> > TwoDArray;

/** Class to encapsulate a set of related alignment information
    and finally construct a contig.

    This class is used by the BatonAssemblerManager to gather
    alignment information about the set of sequences being assembed
    and finally to build a consensus sequence from the subset of
    sequences added to this contig maker object.  Note that for each
    set of related (and aligned) sequences, a new contig maker object
    is expected to be created.  The contig maker is not responsible
    for analyzing and computing the necessary assembly information.
    It only performs the final stages of contig formation by utilizing
    the assembly information provided by other components of the baton
    assembler.

    @see BatonAssemblerManager
*/
class ContigMaker {
public:
    /** The only constructor.

        The constructor is relatively straightforward and adds the
        given root fragment as the first entry in the set of fragments
        to be aligned together.

        \param[in] rootESTidx The index of the root (aka first)
        fragment constituting the contig being assembled.  This value
        must be int the range: 0 \f$ \lte \f$ rootESTidx &lt
        EST::getSequenceCount().
    */
    ContigMaker(const int rootESTidx);

    /** The destructor.

        The destructor does not have any specific tasks to perform (as
        this class currently does not \i directly use any dynamic
        memory) and is present merely to adhere to coding conventions.
    */
    ~ContigMaker();

    /** Obtain the next reference fragment to be used for growing this
        contig.

        This method must be used to determine the next fragment (that
        has already been added to this contig) that must be used as
        the reference to further grow the contig.  This method
        utilizes the information that is setup by the constructor
        and/or that tracks the left-most and right-most fragments in
        this contig.  If (to the left or right).  This method
        alternates between left and right extensions, assuming that
        those fragments have not already been used as references.
        This method uses the nextReference() helper method to
        streamline its operations.

        \return The index of the next fragment in this contig to be
        used as the reference sequence to grow the contig (this value
        is in the range 0 \f$ \lte \f$ retVal &lt
        EST::getSequenceCount()). If both the left-most and right-most
        fragments have already been used as a reference this method
        returns -1.
    */
    int nextReference();

    /** Add a new fragment and its alignment information to this
        contig.

        This method must be used to add a new fragment to be assembled
        as a part of this contig.  The alignment information
        associated with the contig is the only parameter to this
        method.  This method performs the following tasks:

        <ol>

        <li>Checks are made to ensure that that the reference fragment
        in the alignment information is the same as the reference
        fragment set in this contig maker object.</li>

        <li>The newly provided alignment information is added to the
        list of fragments to be assembled.</li>

        <li>The left-most fragment information is checked and updated.</li>

        <li>The right-most fragment information is checked and updated.</li>
        
        </ol>

        \param[in] info The alignment information about the new
        fragment to be added to the list of fragments that are to be
        fused together to form the final contig.
    */
    void add(const AlignmentInfo& info);

    /** Method to trigger contig formation and reporting.

        This method must be used to initiate the final phases of
        contig formation once all the pertinent fragments have been
        added to this contig maker.  This method performs the
        following tasks:

        <ol>

        <li>It computes the base distributions via the
        getBaseDistributions() method.  Refer to the documentation on
        that method for further details.</li>

        <li>It calls the formConsensus() helper method to form the
        actual nucleotide sequences constituting the final
        contig.</li>
        
        </ol>

        \param[in] verbose If this flag is \c true, then additional
        details regarding the information being processed (and
        generated) by this method are displayed on standard output.
        This is useful for debugging and troubleshooting the program.
    */
    void formContig(const bool verbose = false);
	
protected:
    /** Helper method to build the reverse complement nucleotide sequence.

        This is a helper method that is used to build the reverse
        complementary representation of a given nucleotide sequence.
        In other words, given the sequence \c AATCGG this method
        returns \c CCGATT.

        \param[in] seq The nucleotide sequence for which the reverse
        complementary representation is to be computed and returned.

        \return The reverse complementary nucleotide sequence for the
        given seq.  The length of the returned sequence will be
        exactly the same as the given sequence.
    */
    std::string makeRevComp(const std::string& seq) const;

    /** Helper method to compute the base distributions that is
        required to form the contig.

        This method is a helper method that is invoked from the
        formContig() method to compute the base distributions.  Base
        distribution is set of the frequencies with which a given
        nucleotide occurs at a specific column in the final consensus
        sequence.  Here is an example consisting of four aligned
        fragments and the corresponding base distribution:

        <table cellpadding="3" cellspacing="0">
        <th><td>Fragment</td><td>Nucleotides</td></th>

        <tr><td>EST 1: </td><td>\c AACTAGGTACG</td></tr>
        <tr><td>EST 2: </td><td>\c ---TAGGTA--</td></tr>
        <tr><td>EST 3: </td><td>\c AACTCAGG----</td></tr>
        <tr><td>EST 4: </td><td>\c ----TGGTACG</td></tr>
        
        <tr><td colspan="2">Base distributions</td></tr>
        
        <tr><td>Nt-A : </td><td>\c 22002140300</td></tr>
        <tr><td>Nt-T : </td><td>\c 00031003000</td></tr>
        <tr><td>Nt-C : </td><td>\c 00201000020</td></tr>
        <tr><td>Nt-G : </td><td>\c 00000341002</td></tr>
        </table>

        \param[out] distributions A two-dimensional array of integers
        that is populated with the frequency of occurrence of each
        base character.  The base characters are encoded to integers
        in the range 0 to 3 using the ESTCodec::encode() method.

        \param[in] verbose If this flag is \c true, then additional
        details regarding the information being processed (and
        generated) by this method are displayed on standard output.
        This is useful for debugging and troubleshooting the program.
    */
	void getBaseDistributions(TwoDArray& distributions,
							  const bool verbose = false) const;

    /** Helper method to determine the next reference fragment to be
        used.

        This method is used from the main nextReference() method to
        try and explore/expand the contig to the left or to the right.

        \param[in] entry The index (in the alignedESTs vector) of the
        entry suggested as the next reference fragment to be used to
        grow the contig.  Typially either the leftMostEntry or
        rightMostEntry instance variables are passed in as the
        arguments to this parameter.

        \param[out] processed This flag is set to \c true to indicate
        that the entry was used as the next reference entry and it
        should not be reconsidered (as the reference) in the future.
        Typially either the leftMostProcessed or rightMostProcessed
        instance variables are passed in as the arguments to this
        parameter.
    */
	int nextReference(int& entry, bool& processed);

    /** Helper method to form the final consensus sequence.

        This method utilizes the base distributions generated by the
        getBaseDistributions() method to create the final consensus
        sequence for the contig being formed.  The operation of this
        method is relatively straightforward.  For each nucleotide
        position in the base distribution, this method determines the
        base that has the highest frequency and uses that nucleotide
        character as the consensus character.

        \param[in] distributions The base distributions for the
        consensus sequence.  This value is typically obtained via a
        call to the getBaseDistributions() method.

        \return This method returns the final consensus sequence that
        represents the contig assembled.
    */
	std::string formConsensus(const TwoDArray& distributions) const;
	
private:
    /** The list of fragments that have been added to form the contig.

        This array maintains the list of alignment information
        associated with the fragments to be fused together to form the
        final consensus sequence.  The first entry to this list is
        added in the constructor.  Subsequent entries are added each
        time the add() method in this class is called.
    */
    AlignmentInfoList alignedESTs;

    /** The index of the entry in the alignedESTs that is serving as
        the current reference fragment.

        This instance variable tracks the entry in the alignedESTs
        vector that is currently serving as the reference fragment.
        This value is set and updated in the nextReference() method in
        this class.
    */
    int currRefEntry;

    /** The left-most fragment that has been added to this contig
        maker.

        This instance variable tracks the left-most fragment that has
        been added (thus far) to this contig maker.  This value is
        initialized in the constructor (to refer to the root fragment)
        and is later updated in the add() method as new fragements are
        added to this contig maker.  This instance variable always
        tracks the left-most entry in this contig maker.  This
        information is used to grow the contig to the left and to form
        the final consensus sequence for the contig.
    */
    int leftMostEntry;
    
    /** Flag to indicate if the left-most fragment has already been
        used to grow the contig.

        This instance variable is used to determine if the fragment
        indicated by the leftMostEntry has been used as the reference
        sequence to try and grow the contig to the left.  This value
        is initialized in the constructor to \c true (so that the root
        element is used to grow the contig), reset to \c false in the
        nextReference() method, and suitably set to \c true in the
        add() method if a fragment that is more-left than the current
        leftMostEntry is added.
    */
	bool leftMostProcessed;

    /** The right-most fragment that has been added to this contig
        maker.

        This instance variable tracks the right-most fragment that has
        been added (thus far) to this contig maker.  The right-most
        fragment is defined as the fragment whose right-most
        nucleotide occupies the right-most position in the
        contig. This value is initialized in the constructor (to refer
        to the root fragment) and is later updated in the add() method
        as new fragements are added to this contig maker.  This
        instance variable always tracks the right-most entry in this
        contig maker.  This information is used to grow the contig
        to the right and to form the final consensus sequence for the
        contig.
    */    
    int rightMostEntry;

    /** Flag to indicate if the right-most fragment has already been
        used to grow the contig.

        This instance variable is used to determine if the fragment
        indicated by the rightMostEntry has been used as the reference
        sequence to try and grow the contig to the right.  This value
        is initialized in the constructor to \c false (so that the
        root element is not reused to the right as it is already used
        to grow as the left fragment), reset to \c false in the
        nextReference() method, and suitably set to \c true in the
        add() method if a fragment that is more-right than the current
        rightMostEntry is added.
    */    
	bool rightMostProcessed;

    /** Flag to indicate if preference is to be given to the left or
        the right most entry to grow the contig.

        This instance variable is toggled back-and-forth (each time
        the nextReference method is called) from \c true and \c false
        to indicate if the contig is to be grown to the left or to the
        right.  This flag is essentially used to give some preference
        to grow in a specific direction.
    */
	bool leftOrRight;
};

#endif
