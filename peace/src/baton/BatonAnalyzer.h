#ifndef BATON_ANALYZER_H
#define BATON_ANALYZER_H

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

#include "ESTAnalyzer.h"
#include "BatonList.h"
#include "BatonListCache.h"

// Forward declarations to keep compiler fast and happy
class BatonAlignmentInfo;

/** \brief EST Analyzer that uses Baton approach to compute similarity
    between two ESTs.

    <p>This analyzer provides the mechanism to use Batons to compute
    the similarity values between a pair of ESTs. The baton approach
    essentially uses the maximum number of identical batons between
    two given ESTs to determine if a pair of ESTs are related and can
    be clustered together.</p>

	\note This class handles two different code paths running through
	it.  The first case is when this class is used as an analyzer for
	clustering.  The second case is when this class is used with the
	baton assembler.
	
    <p>This class essentially implements the interfaces required by
    the ESTAnalyzer.  This permits the BatonAnalyzer to be used in
    conjunction with a ClusterMaker.</p>

*/
class BatonAnalyzer : public ESTAnalyzer {
    friend class ESTAnalyzerFactory;
	friend class BatonAssembler;
public:
    /** The destructor.
        
        Currently the destructor does not have any specific task to
        perform.  It is merely present to adhere to coding conventions
        and to facilitate continued polymorphism.
    */
    virtual ~BatonAnalyzer();

    /** Add valid command line arguments for this analyzer.

        This method must be used to add all valid command line options
        that are supported by this analyzer.  Note that derived
        classes may override this method to add additional command
        line options that are applicable to it.  This method is
        invoked when the clustering sub-system is initialized.
        
        \note Derived EST analyzer classes may override this method to
        display help for their custom command line arguments.  When
        this method is overridden don't forget to call the
        corresponding base class implementation to add common options.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);

    /** Method to begin EST analysis.
        
        This method is invoked just before commencement of EST
        analysis.  This method first invokes the base class's \c
        initialize method (if the analyzer was created for clustering
        and not for assembly) that initializes the heuristic chain, if
        a chain has been set for this analyzer.

		\return If initialization was successful then this method
        returns \c true. Otherwise this method returns with \c false
        signalling an error.
    */
    virtual bool initialize();

    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        currently saves the index in the \c refESTidx instance
        variable for further look up.  Next, it builds the normal
        BatonList for the reference EST.  The BatonList is stored in
        the blCache (read as baton list cache) instance variable for
        future reference.

        \note This method must be called only after the \c
        initialize() method is called.

		\param[in] est A pointer to an immutable EST object to be used
		as the reference.  The reference EST will be subsequently
		analyzed with a number of other ESTs via the analyze() method.

        \return This method returns zero if the estIdx was within the
        given range of values.  Otherwise this method returns a
        non-zero value as the error code.
    */
    virtual int setReferenceEST(const EST* est);

    /** Set the reference cDNA sequence for analysis and assembly.

        <p>This method is invoked just before a batch of ESTs are
        analyzed and assembled via a call to the assemble() method in
        this class.  This version of the method is typically used to
        align ESTs using a computed consensus sequence for which an
        actual EST entry does not exist (in the list of ESTs loaded
        from a data file).</p>

        <p>This method currently sets the \c refESTidx instance
        variable to -1.  Next, it builds the normal BatonList for the
        supplied reference sequence.  The BatonList is stored in the
        refBatonList instance variable directly for future reference.
        Finally, the coded <i>n</i>-mers for the given cDNA sequence
        is constructed and saved in the \c codedRefNmers instance
        variable for use during analysis, clustering, and assembly.
        
        \note This method must be called only after the \c
        initialize() method is called.
        
        \return This method returns zero if the reference sequence was
        successfully processed.  Otherwise this method returns a
        non-zero value as the error code.
    */
    virtual int setReferenceEST(const char* refSeq);
    
    /** Method to perform exhaustive EST analysis.
        
        This method is used to perform the core tasks of comparing all
        ESTs to one another for full analysis of ESTs.  This is an
        additional feature of PEACE that is \em not used for
        clustering but just doing an offline analysis.  Currently,
        this method merely calls the corresponding base class
        implementation that performs all the necessary operations.

        \note This method was introduced to avoid the unnecessary
        warning about partial overloading of virtual methods (warning
        #654) in ICC.
	
        \return This method returns zero if all the processing
        proceeded successfully. On errors this method returns a
        non-zero value.
    */
    virtual int analyze();

    /** Main method to align a given EST with the reference sequence.

        This method provides the primary API into this analyzer for
        attempting to align a given EST with the reference sequence.
        Note that the reference sequence is set via a call to one of
        the setReferenceEST methods in this class.  This method
        computes the alignment (if relevant) using an internal
        overloaded align() helper method, in the following manner:

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
    virtual bool align(const EST* otherEST, BatonAlignmentInfo& info);
    
protected:
    /** Helper method to align a given EST with the reference sequence.

        This method helps the primary API into this analyzer for
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
        getBestAlignmentInWindow to determine the best possible
        alignment using the identical batons within the given
        window-pair.</li>

        <li>This method tracks the best alignment information using
        the supplied output parameters and returns \c true indicating
        a possible alignment was found.</li>
        
        </ol>

        \param[in] otherEST A pointer to an immutable EST object with
        which the reference EST (or sequence) is to be analyzed and if
        valid/possible aligned.

        \param[out] refAlignPos The index position within the
        reference sequence where the best possible alignment with the
        \c otherEST was found.

        \param[out] othAlignPos The index position within the \c
        otherEST where the best possible alignment with the reference
        EST (or sequence) was found.

        \param[out] score The alignment score (for the best alignment)
        associated with the alignment between \c otherEST and the
        reference sequence at \c othAlignPos and \c refAlignPos
        respectively.

        \param[out] doRC This parameter is set to \c true then a
        reverse complement nucleotide sequence of \c otherEST is used
        to compute candidate alignment.  Otherwise the normal
        nucleotide sequence of \c otherEST is used to compute
        alignment.
        
        \return This method returns \c true if a valid alignment was
        found between the reference sequence and \c otherEST.  If a
        valid alignment was not found (that means the pairs were
        determined to be unrelated), then this method returns \c
        false.
    */
    virtual bool align(const EST* otherEST, int& refAlignPos, int& othAlignPos,
                       int& score, const bool doRC);
    
    /** Analyze and obtain a similarity metric.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.  This method uses the overloaded getMetric() helper
        method and operates as follows:

        <ol>

        <li>It obtains the reference EST's baton list from the cache.</li>

        <li>It uses the overloaded getMetric() helper method to
        compute a similarity metric between the reference EST's normal
        baton list and the other EST's normal baton list.</li>

        <li>It uses the overloaded getMetric() helper method to
        compute a similarity metric between the reference EST's normal
        baton list and the other EST's reverse-complementary baton
        list.</li>

        <li>It returns the maximum of the above two metrics.</li>

        </ol>

        \param[in] otherEST A pointer to an immutable EST object with
        which the reference EST is to be compared.

        \return This method returns the similarity value between the
        referenceEST and otherEST.
    */
    virtual float getMetric(const EST* otherEST);

    /** Helper method to compute a metric by comparing two given baton
        lists.

        This helper method is invoked from the other overloaded
        getMetric() method (typically twice to compute a single
        similarity metric).  This method computes similarity metric in
        the following manner:

        <ol>

        <li>The pairs of sections that have more than \c threshold
        number of identical batons is obtained via the
        BatonList::getWindowPairs() method.</li>

        <li>This method currently returns the number of pairs of
        sections as the similarity metric.</li>

        </ol>

        \note Maybe this method could be modified to use the number of
        identical batons in the candidate windows rather than just the
        count of windows.

        \param[in] list1 The first baton list to be used to determine
        windows that have more than \c threshold number of identical
        batons.  This is typically the "normal" baton list for the
        reference EST.

        \param[in] list1 The first baton list to be used to determine
        windows that have more than \c threshold number of identical
        batons.  This is usually either the normal or the reverse
        complement baton list for the other EST being analyzed.

        \return This method currently returns the number of windows in
        list1 and list2 that have more than \c threshold number of
        identical batons.
    */
    float getMetric(const BatonList* const list1,
                    const BatonList* const list2) const;
    
    /** Method to compare two metrics generated by this class.

        This method provides the interface for comparing metrics
        generated by this ESTAnalyzer when comparing two different
        ESTs.  This method returns \c true if \c metric1 is
        comparatively better than or equal to \c metric2.

        \note As per the ESTAnalyzer API requirements, only EST
        analyzers that are based on distance measures (such as this D2
        analyzer) need to override this method.  However, for
        completeness of the API this method is implemented here to
        perform the same operation as the base class method.
        
        \param[in] metric1 The first metric to be compared against.

        \param[in] metric2 The second metric to be compared against.
        
        \return This method returns \c true if \c metric1 is
        comparatively better than \c metric2.
    */
    inline bool compareMetrics(const float metric1, const float metric2) const
    { return (metric1 > metric2); }
    
    /** Obtain an invalid (or the worst) metric generated by this
        analyzer.
        
        This method can be used to obtain an invalid metric value for
        this analyzer.  This value can be used to initialize metric
        values.
        
        \return This method returns an invalid (or the worst) metric
        of 0 for this EST analyzer.
    */
    inline float getInvalidMetric() const { return 0; }

    /** Get alignment data for the previous call to analyze method.

        This method can be used to obtain alignment data that was
        obtained typically as an byproduct of the previous call to the
        analyze() method.  This method essentially returns the
        difference between the two windows that contained the maximum
        number of identical batons.

        \param[out] alignmentData The parameter is updated to the
        alignment information generated as a part of the the
        immediately preceding analyze(const int) method call is
        returned in the parameter.

        \return This method always returns \c true to indicate that
        alignment data is computed by this ESTAnalyzer.
    */
    virtual bool getAlignmentData(int &alignmentData);
    
    /** Determine if this EST analyzer provides distance metrics or
        similarity metrics.
        
        This method can be used to determine if this EST analyzer
        provides distance metrics or similarity metrics.  If this
        method returns \c true, then this EST analyzer returns
        distance metrics (smaller is better).  On the other hand, if
        this method returns \c false, then this EST analyzer returns
        similarity metrics (bigger is better).

        \return This method returns \c false to indicate that this EST
        analyzer operates using similarity metrics.
    */
    bool isDistanceMetric() const { return false; }

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
        method uses the getAlignmentScore method to compute the
        alignment score for the given pair of batons. It tracks the
        best alignment found and returns the corresponding
        information.  If an alignment score exceeds the
        goodAlignmentScore value, then it immediately returns that
        alignment without further search.</li>

        <li>This method returns the best alignment score (and sets the
        refAlignPos and othAlignPos) determined within the given pair
        of windows.</li>
        
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

        \param[out] refAlignPos The index position in the reference
        cDNA fragment from where the best alignment score was
        determined.  This value corresponds to the starting index
        position of the baton in the reference cDNA which yielded the
        best alignment score.

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

        \param[out] othAlignPos The index position in the reference
        cDNA fragment from where the best alignment score was
        determined.  This value corresponds to the starting index
        position of the baton in the other cDNA which yielded the best
        alignment score.

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

        \return This method returns the best alignment score that was
        found based on the set of identical batons in the given pair
        of windows.
    */
    int getBestAlignmentInWindow(const BatonList& refBatonList,
                                 const IntVector& codedRefNmers,
                                 const int refWindow,
                                 int&  refAlignPos,
                                 const BatonList& othBatonList,
                                 const IntVector& codedOthNmers,
                                 const int othWindow,
                                 int&  othAlignPos,
                                 const int numPermittedErrs,
                                 const int goodAlignmentScore) const;
        
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
    int getAlignmentScore(const IntVector& codedRefNmers, int baton1Pos,
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
	inline int
    tryLeftExtension(const IntVector& codedRefNmers, int refIndexPos,
                     const IntVector& codedOthNmers, int othIndexPos) const {
        // Save starting point to compute delta-moves for return value below.
        const int startRefIndexPos = refIndexPos; 
        while ((refIndexPos > 0) && (othIndexPos > 0)) {
            refIndexPos--; // Move one step left in reference sequence
            othIndexPos--; // Move one step left in other
            if (codedRefNmers[refIndexPos] != codedOthNmers[othIndexPos]) {
                break;
            }
        }
        return startRefIndexPos - refIndexPos;
    }

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
    inline int
    tryRightExtension(const IntVector& codedRefNmers, size_t refIndexPos,
                      const IntVector& codedOthNmers, size_t othIndexPos) const {
        if ((refIndexPos >= codedRefNmers.size()) ||
            (othIndexPos >= codedOthNmers.size())) {
            // One of the indexes is too large. Cannot scan to right.
            return 0;
        }
        // Save starting point to compute delta-moves for return value below.
        const size_t startRefIndexPos = refIndexPos;
        do {
            refIndexPos++; // Move one step right in reference sequence
            othIndexPos++; // Move one step right in other
        } while ((refIndexPos < codedRefNmers.size()) &&
                 (othIndexPos < codedOthNmers.size()) &&
                 (codedRefNmers[refIndexPos] != codedOthNmers[othIndexPos]));
        // Return the delta moved to the right. It will at least be one.
        return refIndexPos - startRefIndexPos;
    }
	
private:
	/** The local cache of Baton Lists.

		This instance variable blCache (read as: baton list cache) is
		used used to maintain a cache of baton list objects that have
		been pre-computed via an on-demand approach.  The entries in
		this list consist of baton lists for both normal and reverse
		complement (RC) representation of a given cDNA fragment.
	*/
	BatonListCache blCache;
	
    /* The constructor for this class.
       
       The constructor is made private so that this class cannot be
       directly instantiated.  However, since the ESTAnalyzerFactory
       is a friend of this class; therefore it can instantiate the
       BatonAnalyzer.  Accordingly, the ESTAnalyzerFactory::create()
       method must be used to instantiate this class.
       
       \param[in] usedForAssembly This parameter is used by the
       BatonAssembler to indicate that this analyzer is being
       exclusively used for assembly.  This class performs slightly
       different initialization tasks depending on the value of this
       parameter.
    */
    BatonAnalyzer(const bool usedForAssembly = false);

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

    /** An alignment score to be considered as a good enough alignment
        score to short circuit exhaustive search of candidate
        alignments.

        This analyzer can be configured to exhaustively compare the
        alignment score between all pairs of identical batons within
        the given windows to determine the best possible alignment.
        However, exhaustive analysis may not be needed if a
        sufficiently good alignment is encountered.  This value
        specifies a sufficiently good alignment score which when
        encountered, this method immediate returns that value/position
        as the best alignment.  The default value is 15.  Setting this
        value of \c -1 will cause this analyzer to exhaustively search
        for best alignment.  This value is a command line argument
        that can be set by the user via the \c --goodScore command
        line argument.
    */
    int goodAlignmentScore;

    /** Flag to indicate if this analyzer is being exclusively used by
        the baton assembler.

        This flag is set in the constructor to indicate if this
        analyzer is being exclusively used by the BatonAssembler to
        aid in assembly process.  When this class is created purely
        for clustering, then this flag is set to \c false.
    */
    const bool usedByAssembler;
};

#endif
