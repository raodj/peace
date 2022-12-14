#ifndef DEFAULT_SEQUENCE_ALIGNER_H
#define DEFAULT_SEQUENCE_ALIGNER_H

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

#include "SequenceAligner.h"
#include "BatonList.h"
#include "Segment.h"

// Forward declarations (if any) to keep compiler happy and fast


/** \def SEED_SEARCH_DIST

	\brief Number of coded n-mers to be checked to the left and right
	to find matching seeds.

	This named constant is used to determine the maximum number of
	coded n-mers to be checked to the left and right of a Segment to
	determine matching seeds. The matching seeds are used to start
	construction of the next segment of alignment.  This constant
	is used by DefaultSequenceAlignner::findSeedToLeft() and
	DefaultSequenceAlignner::findSeedToRight() methods.
*/
#define SEED_SEARCH_DIST 9

/** A simple pair-wise sequence aligner.

    This class provides a simple implementation of a SequenceAligner.
    This 
*/
class DefaultSequenceAligner : public SequenceAligner {
public:
    /** The default (and only) constructor for this class.

        The constructor does not have any special tasks to perform.
    */
    DefaultSequenceAligner();
    
    /** The destructor.
        
        The destructor is merely present to adhere to coding
        conventions.  Currently, this does not have any special tasks
        to perform.
    */
    virtual ~DefaultSequenceAligner();
    
    /** Add the set of command line parameters for this component.
        
        This method is invoked by the SubSystem (that logically owns
        the component) when it is being requested for command line
        arguments.  This method is expected to add the various command
        line arguments that can be used to further customize the
        operations of this component.

        \note The default implementation does not add any arguments.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);

    /** Method to begin EST analysis.
        
        This method is invoked just before commencement of alignment
        operations.  This method first invokes the base class's \c
        initialize method.  It ensures a valid BatonListCache has been
        set via call to setBLCache() method.  In addition, it
        validates command-line arguments.
                
		\return If initialization was successful then this method
        returns \c true. Otherwise this method returns with \c false
        signalling an error.
    */
    virtual bool initialize();

    /** Main method to align a given EST with the reference sequence.

        This method provides the primary API into this aligner for
        attempting to align a given EST with the reference sequence.
        Note that the reference sequence is set via a call to the the
        setReferenceEST() method in this class.  This method computes
        the alignment using an internal overloaded align() helper
        method, in the following manner:

        <ol>

        <li>It uses the internal, overloaded align() helper method to
        compute possible alignment using the normal (not reverse
        complement) nucleotide sequence for \c otherEST.</li>

        <li>If a successful alignment is not found in the previous
        step, then this method once again uses the internal,
        overloaded align() helper method to compute possible alignment
        using the reverse complementary nucleotide sequence for \c
        otherEST.</li>

        <li>If a valid alignment is found in either of the two steps,
        this method updates the supplied output parameters and returns
        \c true. If a valid alignment is not found, then the method
        returns \c false.</li>
        
        </ol>

        \param[in] otherEST A pointer to an immutable EST object with
		which the reference EST (or sequence) is to be analyzed and if
		valid/possible aligned.

        \param[out] alignInfo The information regarding the
        baton-based alignment (if one is computed) is populated into
        this object.  The alignment information includes the following
        information:

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

        \note The setReferenceEST() method must be called (to set the
        reference cDNA fragment to align with) prior to invoking this
        method.
        
        \return This method returns \c true if a valid alignment was
        found between the reference sequence and \c otherEST.  If a
        valid alignment was not found (that means the pairs were
        determined to be unrelated), then this method returns \c
        false.
    */
    virtual bool align(const EST* otherEST, BatonAlignmentInfo& alignInfo);

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
    virtual int setReferenceEST(const EST* est);

    /** Set the reference nucleotide sequence for analysis.

        This method is invoked just before a batch of ESTs are aligned
        via a call to the align(EST *) method.  This method saves the
        reference sequence passed in for further analysis in refSeq
        instance variable.  Next, it builds the normal BatonList for
        the given reference sequence.  The generated BatonList is not
        cached but is stored in refBatonList for further reference.

        \note This method must be called only after the \c
        initialize() method is called.

		\param[in] refSeq The nucleotide sequence to be used as the
		reference sequence for further alignment.

        \return This method always returns zero.
    */
    virtual int setReferenceEST(const std::string& refSeq);
    
protected:
    /** Helper method to align a given EST with the reference
		sequence.

        This method helps the primary API method for this analyzer for
        attempting to align a given EST with the reference sequence.
        Note that the reference sequence is set via a call to one of
        the setReferenceEST methods in this class.  This method
        computes the alignment (if relevant) in the following manner:

        <ol>

        <li>It first computes the normal or reverse complement baton
        list for the \c otherEST (depending on the \c doRC parameter),
        (if it is not already in the cache) via a call to the
        getBatonList method in this class.</li>

        <li>It then determines possible windows that could yield good
        alignment via a call to the BatonList::getWindowPairs()
        method.  If no such window-pair were found, then this method
        immediately returns with the result \c false.</li>

        <li>For each window-pair returned by
        BatonList::getWindowPairs() method, this method (in an
        iterative manner) uses the helper method
        getBestAlignmentInWindow() to determine the best possible
        alignment (and associated Segment objects) using the identical
        batons within the given window-pair.</li>

        <li>If a good alignment is found by the helper then method
        immediately returns with \c true indicating an alignment was
        found.</li>
        
        </ol>

        \param[in] otherEST A pointer to an immutable EST object with
        which the reference EST (or sequence) is to be analyzed and if
        valid/possible aligned.

        \param[out] doRC This parameter is set to \c true then a
        reverse complement nucleotide sequence of \c otherEST is used
        to compute candidate alignment.  Otherwise the normal
        nucleotide sequence of \c otherEST is used to compute
        alignment.
        
        \param[out] alignInfo The baton-based alignment information
        object to be populated by this method if a good alignment was
        found.  Note that the value in this object is meaningful only
        if this method returns \c true.

        \return This method returns \c true if a valid alignment was
        found between the reference sequence and \c otherEST.  If a
        valid alignment was not found (that means the pairs were
        determined to be unrelated), then this method returns \c
        false.
    */
    virtual bool align(const EST* otherEST, const bool doRC,
					   BatonAlignmentInfo& alignInfo);

    /** Helper method that determines best alignment using batons in a
        given window-pair.

        This method is a helper method that is used to by the
        getAlignment() method to determine the best possible alignment
        using the batons within a given a window for the two cDNA
        fragments being aligned.  This method operates in the
        following manner:

        <ol>

        <li>The entries in the given BatonList objects are primarily
        organized based on the <i>n</i>-mers that constitute the baton
        heads.  Accordingly, this method first locates identical (same
        length and identical <i>n</i>-mer heads) batons, with the same
        <i>n</i>-mer heads that fit within the specified windows.</li>

        <li>For each pair of batons that satisfy the aforementioned
        condition (that is, identical batons whose starting index
        position lies within the given, respective windows) this
        method uses the makeSegment() method to compute an initial
        segment of alignment.</li>

        <li>If the initial segment is sufficiently long (longer than
        goodAlignmentScore user-configurable parameter), then this
        method gathers additional Segment objects of alignment into
        the alignInfo parameter and then returns with the score.</li>

        </ol>

        \note The instance variables refBatonList and codedRefNmers
        are used by this method to determine the necessary information
        regarding the reference cDNA fragment with which other
        fragment is being aligned.

        \param[in] refBatonList The baton list associated with the
        reference sequence. The value is typically a reference to the
        same baton list pointed by the refBatonList instance variable.
        This is being passed-in just a convenience (and possibly for
        an extension via an independent class hierarchy at a later
        date).

        \param[in] codedRefNmers A vector containing the encoded list
        of <i>n</i>-mers from the reference cDNA fragment.  The value
        is typically a reference to the same baton list pointed by the
        codedRefNmers list.  This is being passed-in just a
        convenience (and possibly for an extension via an independent
        class hierarchy at a later date).
        
        \param[in] refWindow The zero-based logical window index
        number in the reference baton list.  This value must be in the
        range \f$0 \le refWindow \lt refBatonList.getWindowCount()
        \f$.

        \param[in] othBatonList The baton list containing the various
        batons (in all windows) for the other cDNA fragment which is
        being aligned with the reference fragment.
        
        \param[in] codedOthNmers A vector containing the encoded list
        of <i>n</i>-mers from the other cDNA fragment.  Typically, the
        reference fragment is a given sequence against which
        comparisons for sufficiently similar cDNA fragments is being
        performed.

        \param[in] othWindow The zero-based logical window index
        number in the other baton list.  This value must be in the
        range \f$0 \le othWindow \lt othBatonList.getWindowCount()
        \f$.

        \param[in] numPermittedErrs The number of permitted or
        acceptable errors/differences that are to be tolerated when
        computing the alignments.
        
        \param[in] goodAlignmentScore This method exhaustively
        compares the alignment score between all pairs of identical
        batons within the given windows to determine the best
        alignment.  However, exhaustive analysis may not be needed if
        a sufficiently good alignment is encountered.  This value
        specifies a sufficiently good alignment score which when
        encountered, this method immediate returns that value/position
        as the best alignment.  Setting this value of \c INT_MAX will
        cause this method to exhaustively search for best alignment.

        \param[out] segments The alignment information that is
        computed by this method.  The alignment information is placed
        in this list only if a good-enough alignment is found.
        Consequently, the data in this object is meanigful only if ths
        method returns positive score.
        
        \return This method returns the best alignment score that was
        found based on the set of identical batons in the given pair
        of windows. If a good-enough alignment was not found, then
        this method returns -1.
    */
    int getBestAlignmentInWindow(const std::string& refSeq,
								 const BatonList& refBatonList,
                                 const IntVector& codedRefNmers,
                                 const int refWindow,
								 const std::string& otherSeq,
                                 const BatonList& othBatonList,
                                 const IntVector& codedOthNmers,
                                 const int othWindow,
                                 const int numPermittedErrs,
                                 const int goodAlignmentScore,
                                 SegmentList& segments) const;
        
    /** Determine alignment score using the starting positions of a
        pair of batons in a given pair of encoded cDNA fragments.

        This method is a utility method that is used to determine a
        numerical score indicating the potential quality of aligning
        two cDNA fragments using the starting index positions of a
        pair of \i matching batons (in the two cDNA sequences).  This
        method computes the number of matching nucleotides (in the two
        sets of ccDNA sequences) to the left and right of the given
        baton index positions while permitting a given number of
        differences.  The number of matching nucleotides is returned
        as the alignment score (or quality of alignment) that would be
        obtained using the given pair of batons as the centroid for
        alignment.  This method uses the getLeftExtension and
        getRightExtension methods to explore matching <i>n</i>-mers to
        the left and right of the given baton positions to determine
        the alignment score.

        \param[in] refSeq The reference nucleotide sequence being used
        to determine alignment.  This method does not use this
        parameter but is present for completeness (and for possible
        use by derived class).
        
        \param[in] codedRefNmers A vector containing the encoded list
        of <i>n</i>-mers from a reference cDNA fragment.  Typically,
        the reference fragment is a given sequence against which
        comparisons for sufficiently similar cDNA fragments is being
        performed.  The list of encoded <i>n</i>-mers maintained as an
        instance variable in this class is passed-in as the argument
        to this parameter.

        \param[in] baton1Pos The starting index position of a baton in
        the reference cDNA fragment.  This is the same as the index
        position into the codedRefNmers vector.

        \param[in] otherSeq The other nucleotide sequence that is
        being aligned with the reference sequence.  This method does
        not use this parameter but is present for completeness (and
        for possible use by derived class).
        
        \param[in] codedOthNmers A vector containing the encoded list
        of <i>n</i>-mers from another cDNA fragment that is being
        aligned with the reference cDNA fragment.  The list of encoded
        <i>n</i>-mers is obtained via a call to
        BatonList::getCodedNMers() method.
        
        \param[in] baton2Pos The starting index position of a baton in
        the other cDNA fragment.  This is the same as the index
        position into the codedOthNmers vector.

        \param[in] numPermittedErrs The number of permitted or
        acceptable errors/differences that are to be tolerated when
        computing the alignment score.

        \return The alignment score between the two cDNA fragments,
        given the baton positions.  The score is the length of the
        matching nucleotide sequences around (to the left and right)
        the given baton positions such that there are no more than
        numPermittedErrs differing nucleotides.
    */
    int getAlignmentScore(const std::string& refSeq,
                          const IntVector& codedRefNmers, int baton1Pos,
                          const std::string& otherSeq,
                          const IntVector& codedOthNmers, int baton2Pos,
                          int numPermittedErrs) const;
    
    /** Helper method to explore the alignment quality to the left of
        the given starting index positions.

        This helper method is invoked from the getAlignmentScore()
        method to determine how well two given ESTs match each other
        from a given starting position.  This method essentially
        searches for a position where the encoded <i>n</i>-mers
        different.  The search starts at the given index positions in
        the two cDNA fragments.  The search proceeds as long as there
        are at least one more <i>n</i>-mer is present (in both ESTs)
        to the left of the given starting index position and the coded
        <i>n</i>-mers are identical.

        \note This method first moves one step to the left (if
        possible) and then performs the comparison.
        
        \param[in] codedRefNmers A vector containing the encoded list
        of <i>n</i>-mers from a reference cDNA fragment.  Typically,
        the reference fragment is a given sequence against which
        comparisons for sufficiently similar cDNA fragments is being
        performed.

        \param[in] refIndexPos The starting index position in the
        reference codedRefNmers vector from where the search (to the
        left) must commence.

        \param[in] codedOthNmers A vector containing the encoded list
        of <i>n</i>-mers from another cDNA fragment that is being
        aligned with the reference cDNA fragment.

        \param[in] refIndexPos The starting index position in the
        codedOthNmers vector from where the search (to the left) must
        commence.

        \return If one of the refIndexPos is already at the left-most
        <i>n</i>-mer entry, then this method returns zero.  Otherwise,
        this method returns the number of nucleotides that are similar
        to the left of the starting positions in both the
        codedRefNmers and codedOthNmers array.  In this case, the
        minimum return value is one.
    */
	int
    tryLeftExtension(const IntVector& codedRefNmers, int refIndexPos,
                     const IntVector& codedOthNmers, int othIndexPos) const;

    /** Helper method to explore the alignment quality to the right of
        the given starting index positions.

        This helper method is invoked from the getAlignmentScore()
        method to determine how well two given ESTs match each other
        from a given starting position.  This method essentially
        searches for a position where the encoded <i>n</i>-mers
        different.  The search starts at the given index positions in
        the two cDNA fragments.  The search proceeds as long as there
        are at least one more <i>n</i>-mer is present (in both ESTs)
        to the right of the given starting index positions and the
        coded <i>n</i>-mers are identical.

        \note This method first moves one step to the right (if
        possible) and then performs the comparison.
        
        \param[in] codedRefNmers A vector containing the encoded list
        of <i>n</i>-mers from a reference cDNA fragment.  Typically,
        the reference fragment is a given sequence against which
        comparisons for sufficiently similar cDNA fragments is being
        performed.

        \param[in] refIndexPos The starting index position in the
        reference codedRefNmers vector from where the search (to the
        right) must commence.

        \param[in] codedOthNmers A vector containing the encoded list
        of <i>n</i>-mers from another cDNA fragment that is being
        aligned with the reference cDNA fragment.

        \param[in] refIndexPos The starting index position in the
        codedOthNmers vector from where the search (to the right) must
        commence.

        \return If one of the refIndexPos is already at the right-most
        <i>n</i>-mer entry, then this method returns zero.  Otherwise,
        this method returns the number of nucleotides that are similar
        to the right of the starting positions in both the
        codedRefNmers and codedOthNmers array.  In this case, the
        minimum return value is one.
    */    
	int
    tryRightExtension(const IntVector& codedRefNmers, int refIndexPos,
                      const IntVector& codedOthNmers, int othIndexPos) const;

    /** Pointer to the reference EST to align against.

        This member object is used to hold a pointer to a given
        reference EST.  This member is initialized to NULL in the
        constructor and is changed by the setReferenceEST() methods.
        This pointer can be NULL if a direct reference sequence has
        been set (instead of a reference EST).
    */
    const EST* refEST;

    /** The reference sequence for align against.

        This member object is used to hold the reference nucleotide
        sequence against which ESTs are to be aligned. This member
        always as a reference sequence immaterial of which one of the
        setReferenceEST() method is used. It either points to the
        nucleotide sequence of refEST (if refEST is not NULL) or the
        reference sequence that has been directly set.
    */
    std::string refSeq;
    
    /** The baton list for the reference cDNA sequence.

        This instance variable is used to maintain the baton list
        associated with the reference cDNA sequence with which other
        fragments are to be analyzed, for clustering and assembly.
        This list is set whenever one of the overloaded
        setReferenceEST() methods in this class are invoked.
    */
    BatonList* refBatonList;

    /** The vector containing the encoded list of <i>n</i>-mers from
        the reference cDNA.

        This vector is used to hold the encoded list of <i>n</i>-mers
        from the reference cDNA fragment.  The reference cDNA fragment
        is a given sequence with which other fragments are analyzed,
        for clustering and assembly.  The list of encoded
        <i>n</i>-mers is obtained via a call to
        BatonList::getCodedNMers() method.  This list is maintained as
        an instance variable to minimize variables passed-in to
        various methods.  The list is not cached in BatonList class
        because its usage lifetime is rather short -- it is needed
        only during alignment, once we have determined two cDNA
        fragments must be aligned.  It is not needed for general
        searching.
    */
    IntVector codedRefNmers;
    
    /** Instance variable to track the minimum number of identical
        batons in a window for two cDNA fragments to be deemed
        related.

        This member is used to hold the threshold value.  The
        threshold value is a command line argument that can be set by
        the user via the \c --batonThresh command line argument.
    */
    int threshold;

    /** The number of permitted or acceptable errors/differences that
        are to be tolerated when computing the alignments.

        This member is used to hold the number of errors/differences
        (in number of nucleotides) that are to be tolerated when
        exploring for candidate alignments.  This value is ultimately
        used by the getAlignmentScore() method while exploring
        candidate alignment lengths by exploring matches to the left
        and right of a given possible alignment position pair.  the
        positions are determined by the starting index value of
        identical batons that are being used to compute candiate
        alignment.  The default value is 10.  This value is a command
        line argument that can be set by the user via the \c
        --permErrs command line argument.
    */
    int numPermittedErrs;

    /** An alignment score, obtained using a single pair of batons, to
        be considered as a good enough alignment score to short
        circuit exhaustive search of candidate alignments.

        <p>This analyzer can be configured to exhaustively compare the
        alignment score between all pairs of identical batons within
        the given windows to determine the best possible alignment.
        However, exhaustive analysis may not be needed if a
        sufficiently good alignment is encountered.  This value
        specifies a sufficiently good alignment score which when
        encountered, this method immediate returns that value/position
        as the best alignment.  The default value is 25.  Setting this
        value of \c -1 will cause this analyzer to exhaustively search
        for best alignment.  This value is a command line argument
        that can be set by the user via the \c --goodScore command
        line argument.</p>

		<p>Note the difference between this parameter and
		minAlignScore. This parameter is the score (currently, number
		of matching nucleotides) between a single-pair of contiguous
		fragments of cDNA from the reference EST and other EST being
		aligned.  On the other hand, the minAlignScore is the sum of
		the scores from all pairs of aligned cDNA fragments from the
		reference EST and other EST.</p>
   
    */
    int goodAlignmentScore;

	/** An alignment score, obtained from full alignment (using all
		segments) to be considered as good enough alignment (thereby
		short circuiting further/exhaustive search).

		<p>This parameter determines when a complete alignment is
		considered sufficiently good to be short circuit further
		search for alignments. The default value is 15.  This value
		can become larger as the average length of cDNA fragments to
		be aligned increases.  However, care must be taken when
		setting it to large values as it may cause good alignments to
		be rejected resulting in too many singletons.  This value is
		command line argument that can be set by the user via the \c
		--minAlignScore command line argument.</p>
		
		<p>Note the difference between this parameter and
		goodAlignmentScore. This parameter is the score (currently,
		number of matching nucleotides) between all pairs of aligned
		sub-fragments. In contrast, the goodAlignmentScore is just
		between one pair of sub-fragments and is used to decide if
		constructing a full alignment is worthwhile.</p>
	*/
	int minAlignScore;
	
private:
	/** This is a custom helper method that is used to search to the
        left and right of an identical baton-pair starting positions
        to form the longest Segment possible.

        <p>This helper method is used to form a Segment that
        identifies a sequence of identical nucleotides starting from a
        given pair of batons. This method searches to both the left
        and the right of the starting position. However, the span of
        the Segment will not exceed the specified bounds. This feature
        is used by this class to repeatedly form Segments after an
        inital, sufficiently long (at least as long as
        DefaultSequenceAligner#goodAlignmentScore) Segment is
        formed.</p>

        <p>This method essentially uses the countMatchingBasesToLeft()
        and countMatchingBasesToRight() method to form the
        Segment. Consequently this method by itself is relatively
        straightforward.</p>

        \attention We use a cusomtimzed methods for left and right
        exploration because these methods are frequently used method
        and we would like to keep the code as small and customized as
        possible to get maximum performance.
        
        \param[in] refSeq The nucleotide sequence for the reference
        cDNA fragment against which we are trying to locate other
        candidate sequences for matching.

        \param[in] baton1Pos The starting position of the baton in the
        refSeq from where the search to the left in refSeq must
        commence.

        \param[in] refSeqStartPos The starting position in refSeq
        where search must terminate. This value must not be less than
        0 (zero).

        \param[in] refSeqEndPos The ending position in refSeq where
        search must terminate. This value cannot exceed the length of
        the nucleotide sequence (i.e., refSeq.size()).
        
        \param[in] otherSeq The nucleotide sequence for the other cDNA
        fragment which is being tested for candidacy in current
        consensus sequence.

        \param[in] baton2Pos The starting position of an
        identical-baton in the otherSeq from where the search to the
        left in otherSeq must commence.

        \param[in] othSeqStartPos The starting position in otherSeq
        where search must terminate. This value is typically 0 (zero).

        \param[in] othSeqEndPos The ending position in otherSeq where
        search must terminate. This value cannot exceed the length of
        the nucleotide sequence (i.e., refSeq.size()).
        
        \return This method returns the left-distance from baton1Pos
        where the first, two consecutively differing nucleotides were
        found.  For example, the following method call (nucleotide
        sequences intentionally camel-cased and indented to help
        illustration) should return 8.
    */
    Segment makeSegment(const std::string& refSeq, const int baton1Pos,
                        const int refSeqStartPos, const int refSeqEndPos,
                        const std::string& otherSeq, const int baton2Pos,
                        const int othSeqStartPos, const int othSeqEndPos) const;

	/** This is a custom helper method that is used to continue to
        construct segments to the left-and-right of one segment with
        sufficiently long alignment.

        <p>This helper method is used to form multiple segments that
        identify the series of segments with good alignment between
        the reference and other sequence. This method assumes that an
        initial reference segment (with maximum contiguous alignment)
        has already been constructed.  This method then searches to
        the left and then to the right of the reference segment to
        identify other sub-alignment segments. </p>

        \attention We use a cusomtimzed methods for left and right
        exploration because these methods are frequently used method
        and we would like to keep the code as small and customized as
        possible to get maximum performance.
        
        \param[in] refSeq The nucleotide sequence for the reference
        cDNA fragment against which we are trying to locate other
        candidate sequences for matching.

        \param[in] codedRefNmers A vector containing the encoded list
        of <i>n</i>-mers from another cDNA fragment that is being
        aligned with the reference cDNA fragment.

		\param[in] otherSeq The nucleotide sequence for the other cDNA
        fragment which is being aligned with the refSeq.
		
        \param[in] codedOthNmers A vector containing the encoded list
        of <i>n</i>-mers from another cDNA fragment that is being
        aligned with the reference cDNA fragment.

		\param[in] refSeg The reference segment that has already been
		formed indicating the region of maximum alignment between the
		reference and other sequences.  This segment is appropriately
		added to the segments (outgoing vector/list of segments)
		parameter.

		\param[out] segements This list is populated with segments (as
		they are formed) by this method. Any existing entries are
		initially cleared out prior to adding segments.  When this
		method returns this vector will containin the list of sgements
		formed by this method.  This vector/list is sorted with
		segments appearing in a left-to-right manner. The refSeg is
		appropriately included in this list.  Consequently, at a bare
		minimum, this list has at least one entry when this method
		returns.
		
        \return This method returns the count of total number of
        nucleotides that were aligned into segments. This value is
        essentially the sum of the lengths of Segments constituting
        the alignment between the reference and other cDNA sequences.
    */
    int makeMultipleSegments(const std::string& refSeq,
                             const IntVector& codedRefNmers,
                             const std::string& otherSeq,
                             const IntVector& codedOthNmers,
                             const Segment& refSeg,
                             SegmentList& segments) const;
    
	/** This is a custom helper method that is used to search to the
        left of an identical baton-pair starting positions to locate
        first position where two consecutive nucleotide differences
        occur.

        This helper method is used to determine the maximum left-most
        initial position of similarity between the reference sequence
        and another sequence for determining the best possible
        alignment using a given pair of identical batons from the two
        sequences.  This method searches for two consecutively
        different nucletodies (one difference may indicate a SNP and
        we intentionally skip over them) in the given pair of
        sequences.  This method proceeds towards the left of the two
        sequences (until the beginning of either of the sequence is
        hit).

        \note We use a cusomtimzed methods for left and right
        exploration because these methods are frequently used method
        and we would like to keep the code as small and customized as
        possible to get maximum performance.
        
        \param[in] refSeq The nucleotide sequence for the reference
        cDNA fragment against which we are trying to locate other
        candidate sequences for matching.

        \param[in] baton1Pos The starting position of the baton in the
        refSeq from where the search to the left in refSeq must
        commence.

        \param[in] refSeqStartPos The starting position in refSeq
        where search must terminate. This value is typically 0
        (zero).

        \param[in] otherSeq The nucleotide sequence for the other cDNA
        fragment which is being tested for candidacy in current
        consensus sequence.

        \param[in] baton2Pos The starting position of an
        identical-baton in the otherSeq from where the search to the
        left in otherSeq must commence.

        \param[in] othSeqStartPos The starting position in otherSeq
        where search must terminate. This value is typically 0 (zero).
        
        \return This method returns the left-distance from baton1Pos
        where the first, two consecutively differing nucleotides were
        found.  For example, the following method call (nucleotide
        sequences intentionally camel-cased and indented to help
        illustration) should return 8.

        \code
        
        countMatchingBasesToLeft("atcgatCGatcGaTcgatcgatcg",
                                 16, 0, 
                                 "atACatcTaAcgatcgatcgatcg",
                                 12, 0);


        \endcode

        Note that if even a single step cannot be made, then this
        method returns 0 (zero).
    */
	int countMatchingBasesToLeft(const std::string& refSeq,
								 const int baton1Pos,
								 const int refSeqStartPos,
								 const std::string& otherSeq,
								 const int baton2Pos,
								 const int othSeqStartPos) const;

	/** This is a custom helper method that is used to search to the
        right of an identical baton-pair starting positions to locate
        first position where two consecutive nucleotide differences
        occur.
        
        This helper method is used to determine the maximum right-most
        initial position of similarity between the reference sequence
        and another sequence for determining the best possible
        alignment using a given pair of identical batons from the two
        sequences.  This method searches for two consecutively
        different nucletodies (one difference may indicate a SNP and
        we intentionally skip over them) in the given pair of
        sequences.  This method proceeds towards the right of the two
        sequences (until the end of either of the sequence is hit).

        \note We use a cusomtimzed methods for left and right
        exploration because these methods are frequently used method
        and we would like to keep the code as small and customized as
        possible to get maximum performance.
        
        \param[in] refSeq The nucleotide sequence for the reference
        cDNA fragment against which we are trying to locate other
        candidate sequences for matching.

        \param[in] baton1Pos The starting position of the baton in the
        refSeq from where the search to the right in refSeq must
        commence.

        \param[in] refSeqEndPos The ending position in refSeq where
        search must terminate. This value is typically the length of
        the nucleotide sequence (i.e., refSeq.size()).

        \param[in] otherSeq The nucleotide sequence for the other cDNA
        fragment which is being tested for candidacy in current
        consensus sequence.

        \param[in] baton2Pos The starting position of an
        identical-baton in the otherSeq from where the search to the
        right in otherSeq must commence.

        \param[in] othSeqEndPos The ending position in otherSeq where
        search must terminate. This value is typically the length of
        the nucleotide sequence (i.e., refSeq.size()).
        
        \return This method returns the right-distance from baton1Pos
        where the first, two consecutively differing nucleotides were
        found.  For example, the following method call (nucleotide
        sequences intentionally camel-cased and indented to help
        illustration) should return 6.

        \code
        
        countMatchingBasesToRight("nnnnatCgaTCgatcgatcgatcg",
                                  4, 24, 
                                      "atGgaCGgatcgatcgatc",
                                  0, 20);


        \endcode

        Note that if even a single step cannot be made, then this
        method returns 0 (zero).
    */
	int countMatchingBasesToRight(const std::string& refSeq,
                                  const int baton1Pos,
                                  const int refSeqStartPos,
                                  const std::string& otherSeq,
                                  const int baton2Pos,
                                  const int othSeqStartPos) const;

    /** Helper method to determine start (or seed) for a segment to
        the left of a given pair of positions in two cDNA fragments.

        <p>This is a helper method that is used to locate the start of
        a new sub-alignment segment to the left of a given pair of
        starting positions on the reference and other cDNA fragments.
        This method is repeatedly invoked from the
        makeMultipleSegments() method in this class to try and form
        multiple sub-alignment segments.</p>
        
        This method essentially searches for the first position where
        the encoded <i>n</i>-mers are \b not different.  The two
        identical coded <i>n</i>-mers are called \i seeds in this
        context. The seed positions are used to form the next
        sub-alignment Segment. The search starts at the given pair of
        positions in the two cDNA fragments and proceeds to the left.
        The search does not exceed 10 bases from the given starting
        position.

        \note This method first moves one step to the left (if
        possible) and then performs the comparison.
        
        \param[in] codedRefNmers A vector containing the encoded list
        of <i>n</i>-mers from a reference cDNA fragment.  Typically,
        the reference fragment is a given sequence against which
        comparisons for sufficiently similar cDNA fragments is being
        performed.

        \param[in] refStartPos The starting index position in the
        reference codedRefNmers vector from where the search (to the
        left) must commence.

        \param[in] codedOthNmers A vector containing the encoded list
        of <i>n</i>-mers from another cDNA fragment that is being
        aligned with the reference cDNA fragment.

        \param[in] othStartPos The starting index position in the
        codedOthNmers vector from where the search (to the left) must
        commence.

        \param[out] refSeedPos The seed position on the reference cDNA
        fragment from where the next Segment can be formed.  This
        value is meaningful only if this method returns \c true to
        indicate that a valid seed was found.

        \param[out] othSeedPos The seed position on the other cDNA
        fragment from where the next Segment can be formed.  This
        value is meaningful only if this method returns \c true to
        indicate that a valid seed was found.

        \return This method returns true if a valid pair of seed were
        found to the left of the starting positions.  If a valid seed
        was not found, then this method returns false.
    */
	bool findSeedToLeft(const IntVector& codedRefNmers, const int refStartPos,
                        const IntVector& codedOthNmers, const int othIndexPos,
                        int& refSeedPos, int& othSeedPos) const;

    /** Helper method to determine start (or seed) for a segment to
        the right of a given pair of positions in two cDNA fragments.

        <p>This is a helper method that is used to locate the start of
        a new sub-alignment segment to the right of a given pair of
        starting positions on the reference and other cDNA fragments.
        This method is repeatedly invoked from the
        makeMultipleSegments() method in this class to try and form
        multiple sub-alignment segments.</p>
        
        This method essentially searches for the first position where
        the encoded <i>n</i>-mers are \b not different.  The two
        identical coded <i>n</i>-mers are called \i seeds in this
        context. The seed positions are used to form the next
        sub-alignment Segment. The search starts at the given pair of
        positions in the two cDNA fragments and proceeds to the right.
        The search does not exceed 10 bases from the given starting
        position.

        \note This method first moves one step to the right (if
        possible) and then performs the comparisons.
        
        \param[in] codedRefNmers A vector containing the encoded list
        of <i>n</i>-mers from a reference cDNA fragment.  Typically,
        the reference fragment is a given sequence against which
        comparisons for sufficiently similar cDNA fragments is being
        performed.

        \param[in] refStartPos The starting index position in the
        reference codedRefNmers vector from where the search (to the
        right) must commence.

        \param[in] codedOthNmers A vector containing the encoded list
        of <i>n</i>-mers from another cDNA fragment that is being
        aligned with the reference cDNA fragment.

        \param[in] othStartPos The starting index position in the
        codedOthNmers vector from where the search (to the right) must
        commence.

        \param[out] refSeedPos The seed position on the reference cDNA
        fragment from where the next Segment can be formed.  This
        value is meaningful only if this method returns \c true to
        indicate that a valid seed was found.

        \param[out] othSeedPos The seed position on the other cDNA
        fragment from where the next Segment can be formed.  This
        value is meaningful only if this method returns \c true to
        indicate that a valid seed was found.

        \return This method returns true if a valid pair of seed were
        found to the right of the given pair of starting positions.
        If a valid seed was not found, then this method returns false.
    */
	bool findSeedToRight(const IntVector& codedRefNmers, const int refStartPos,
                         const IntVector& codedOthNmers, const int othIndexPos,
                         int& refSeedPos, int& othSeedPos) const;

    bool findNearestSeed(const IntVector& codedRefNmers, const int refStartPos,
                         const int refEndPos,
                         const IntVector& codedOthNmers, const int othStartPos,
                         const int othEndPos,
                         int& refSeedPos, int& othSeedPos) const;

    /** \typedef std::pair<int, int> SeedInfo

        \brief A typedef for convenient reference to seed information
        in the findNearestSeed() method.

        This is just a convenient typedef that is used to make the
        code more readable. The pair of integers (SeedInfo.first,
        SeedInfo.second) hold the coded n-mer value and its index
        position respectively. This information is used to create a
        sorted std::vector and locate matching seeded that are nearest
        to each other quickly.
    */
    typedef std::pair<int, int> SeedInfo;

    /** Comparison structure for SeedInfo
        
        The following structure essentially provides the comparison
        operator needed for sorting and other operations that require
        comparing values. The comparison is done using only the coded
        n-mer value (while ignoring its index position).
    */
    struct LessSeedInfo {
        inline bool operator() (const SeedInfo& si1,
                                const SeedInfo& si2) const {
            // Compare (and sort) based on n-mer values.
            return si1.first < si2.first;
        }
    };
    
};

#endif
