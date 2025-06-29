#ifndef D2_H
#define D2_H

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

#include "FWAnalyzer.h"
#include "ESTList.h"

#include <string>
#include <vector>

// Forward declaration to keep compiler happy
class EST;
class ResultLog;

/** \brief EST Analyzer that uses the D2 algorithm to compute
    distances between two ESTs.

    <p>This analyzer provides the mechanism to use D2 algorithm to
    compute the distance values between a pair of ESTs. The D2
    implementation has been adapted purely from the implementations of
    WCD, Zimmerman, and CLU.<p>

    \note This D2 analyzer uses a word size of 8 base pairs.  This is
    hard coded into the algorithm in the form of various assumptions
    (data types of variables etc.). A similar approach is used by
    other D2 implementations as well.  This is a compromise between
    performance and flexibility of the implementation.
*/
class D2 : public FWAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~D2();

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
        initialize method that initializes the heuristic chain, if a
        chain has been set for this analyzer. It then initializes
        memory for the word tables used to store hash values for the
        pairs of ESTs to be analyzed.

        \return Currently, this method always returns 0 (zero) to
        indicate initialization was successfully completed.
    */
    virtual bool initialize();

    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        currently saves the index in the \c refESTidx instance
        variable for further look up.  Next it creates a "word table"
        (via a call to the \c buildWordTable method) mapping integer
        indices to integer hashes of words, in effect translating the
        sequence from a sequence of characters to a sequence of
        n-words (where n = wordSize).  This word table is kept until
        the reference EST is changed, which reduces overhead.

        \note This method must be called only after the \c
        initialize() method is called.

		\param[in] est A pointer to an immutable EST object to be used
		as the reference.  The reference EST will be subsequently
		analyzed with a number of other ESTs via the analyze() method.

        \return This method returns \c true if the estIdx was within
        the given range of values.  Otherwise this method returns a
        non-zero value as the error code.
    */
    virtual int setReferenceEST(const EST* est);
    
    /** Method to perform exhaustive EST analysis.
        
        This method is used to perform the core tasks of comparing all
        ESTs to one another for full analysis of ESTs.  This is an
        additional feature of PEACE that is \em not used for clustering
        but just doing an offline analysis.  Currently, this method
        merely calls the corresponding base class implementation that
        performs all the necessary operations.

        \note This method was introduced to avoid the unnecessary
        warning about partial overloading of virtual methods (warning
        #654) in ICC.
	
        \return This method returns zero if all the processing
        proceeded successfully. On errors this method returns a
        non-zero value.
    */
    virtual int analyze() { return FWAnalyzer::analyze(); }
    
protected:
    /** Analyze and obtain a distance metric.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \param[in] otherEST Pointer to an immutable EST object with
        which the reference EST is to be compared.

        \return This method returns the distance value reported by the
        D2 algorithm.
    */
    virtual float getMetric(const EST* otherEST);

    /** Method to compare two metrics generated by this class.

        This method provides the interface for comparing metrics
        generated by this ESTAnalyzer when comparing two different
        ESTs.  This method returns \c true if \c metric1 is
        comparatively better than or equal to \c metric2.

        \note As per the ESTAnalyzer API requirements, EST analyzers
        that are based on distance measures (such as this D2 analyzer)
        \b must override this method.
        
        \param[in] metric1 The first metric to be compared against.

        \param[in] metric2 The second metric to be compared against.
        
        \return This method returns \c true if \c metric1 is
        comparatively better than \c metric2.
    */
    inline bool compareMetrics(const float metric1, const float metric2) const
    { return (metric1 < metric2); }
    
    /** Obtain an invalid (or the worst) metric generated by this
        analyzer.
        
        This method can be used to obtain an invalid metric value for
        this analyzer.  This value can be used to initialize metric
        values.
        
        \note Derived distance-based metric classes (such as this D2
        analyzer) \b must override this method to provide a suitable
        value.
        
        \return This method returns an invalid (or the worst) metric
        of 1e7 for this EST analyzer.
    */
    inline float getInvalidMetric() const { return 4.0f * frameSize; }

    /** Get alignment data for the previous call to analyze method.

        This method can be used to obtain alignment data that was
        obtained typically as an byproduct of the previous call to the
        analyze() method.  This method essentially returns the
        difference between the two windows that provided the minimum
        d2 distance value.

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

        \return This method returns \c true to indicate that this EST
        analyzer operates using distance metrics.
    */
    bool isDistanceMetric() const { return true; }

    /** Helper method to create a word table.
        
        Creates a "word table" mapping integer indices to integer
        hashes of words, in effect translating the sequence from a
        sequence of characters to a sequence of n-words (where n =
        wordSize).

        \param encoder The encoder class to be used for building the
        word table.
		
        \param[out] wordTable The word table to be populated with with
        hash values.

        \param[in] estSeq The sequence of base pairs associated with
        this EST that must be converted to hash values.
    */
    template <typename Encoder>
    void buildWordTable(std::vector<int>& wordTable, const char* estSeq,
                        Encoder encoder) {
        ASSERT( ESTAnalyzer::estList != NULL );
        // First clear out old data in word table.
        wordTable.clear();
        const size_t wordTableSize = estList->getMaxESTLen() + frameSize;
        // Reserve space to avoid repeated reallocations as entries
        // are added to the wordTable vector.
        wordTable.reserve(wordTableSize);
        // Compute the has for the first word using the encoder. The
        // encoder may do normal or reverse-complement encoding for us.
        int hash       = 0;
        int ignoreHash = 0;
        for(register int i = 0; ((*estSeq != 0) && (i < wordSize - 1));
            i++, estSeq++) {
            hash = encoder(hash, *estSeq, ignoreHash);
        }
        // Now compute the hash for each word
        for(int entry = 0; (*estSeq != 0); entry++, estSeq++) {
            hash = encoder(hash, *estSeq, ignoreHash);
            if (!ignoreHash) {
                // This hash does not have a 'n' in it. So use it.
                wordTable.push_back(hash);
            }
        }
    }

    /** Instance variable that maps an index in the reference sequence
        (sequence s1) to the hash of a word.  This hash can then be
        used as an index in the delta array to get the frequency
        differential for that word.

        The word table is created in the initialize() method and
        filled in using the buildWordTable() method.  For the reference
        EST, buildWordTable() is called in the setReferenceEST() method,
        meaning it does not need to be rebuilt every time we analyze
        a new comparison sequence.
    */
    std::vector<int> s1WordTable;

    /** Instance variable that maps an index in the comparison sequence
        (sequence s2) to the hash of a word.

        This hash is used as an index in the delta array to get the
        frequency differential for that word.

        The word table is created in the initialize() method and
        filled in using the buildWordTable() method.  For the comparison
        EST, buildWordTable() must be called in the analyze() method because
        a new comparison sequence is given every time analyze() is called.
    */
    std::vector<int> s2WordTable;

    /** Array to track the delta values in the core D2 algorithm.

        This array is initialized to point to an array of size
        4<sup>wordSize</sup>.  This array is used to track the delta
        values generated in the core D2 algorithm.
    */
    int *delta;
    
    /** Parameter to define number of characters to shift the frame
        on the reference sequence, when computing D2.
		
        This parameter is used to enable D2-asymmetric behavior.
        The default value is 1, which means D2 symmetric: all frames
        in both sequences will be compared.  Higher values mean that
        the algorithm will shift by more than one character when
        shifting the frame on the reference sequence, resulting in
        fewer computations but a possible loss of accuracy from not
        comparing every frame in both sequences.
		
        \note Currently this value is not used.
    */
    int frameShift;
    
    /** Instance variable to store the number of bits to be shifted to
        create hash values.

        <p>This instance variable is set to the value of 2 * (\em
        wordSize - 1) (in the \c initialize method) to reflect the
        number of bits that need to be shifted in order to build the
        hash values for common words (including the values stored in
        \c s1WordMap and \c s1RCWordMap).</p>

        <p>This instance variable is actually passed on to the
        ESTCodec::NormalEncoder or ESTCodec::RevCompEncoder when
        computing hash values.  Since this is value is passed in a
        template parameter, it is defined to be static (to ensure that
        it has external linkage as per the ISO/ANSI standards
        requirement).</p>
    */    
    static int bitShift;

    /** The threshold score below which two ESTs are considered
        sufficiently similar to be clustered.

        This instance variable tracks the threshold value to be used
        to break out of the core D2 loop.  This value essentially
        represents the D2 score below which two ESTs are considered
        sufficiently similar. Currently, the default threshold value
        is set to 40.  This value can be overridden by the user with
        the \c --threshold command line argument.
    */
    int threshold;
    
private:
    /* The default constructor for this class.
       
       The default constructor for this class.  The constructor is
       made private so that this class cannot be directly
       instantiated.  However, since the ESTAnalyzerFactory is a
       friend of this class; therefore it can instantiate the D2
       analyzer.  Accordingly, the ESTAnalyzerFactory::create()
       method must be used to instantiate this class.
    */
    D2();

    /** Instance variable to track alignment metric computed by the
        analyze() method.
        
        This instance variable is used to hold the alignment metric
        that was computed in the previous \c analyze method call. By
        default this value is set to zero.  The alignment metric is
        computed as the difference in the window positions (on the two
        ESTs being analyzed) with the minimum d2 distance.
    */
    int alignmentMetric;

    /** The bitmask to be used when build hashing values.
            
        This bitmask is used to drop higher order bits when building
        hash values for \c wordSize (defined in base class) base pairs
        in a given EST. This instance variable is actually passed on
        to the ESTCodec::NormalEncoder or ESTCodec::RevCompEncoder
        when computing hash values.  Since this is value is passed in
        a template parameter, it is defined to be static (to ensure
        that it has external linkage as per the ISO/ANSI standards
        requirement).
    */
    static int BitMask;

    /** Instance variable to track the number of words (of \c
        wordSize) that can fit into a window (of \c frameSize).

        This instance variable is set in the initialize() method and
        its value is used by various methods in this class.  Rather
        than computing this value repeatedly, it is computed once and
        used throughout.
    */
    int numWordsInWindow;

    /** Helper method to update the scores based on a sliding window.

        This method is invoked from several different spots from the
        runD2() method to update the d2 scores as the window slides
        across the two sequences being analyzed. The hash values of
        the words moving into and out of the window are used to update
        the scores.

        \param[in] wordIn The hash value of the word moving into the
        window.

        \param[in] wordOut The hash value of the word moving out of
        the window.

        \param[in,out] score The current running score for this
        window. This value is updated using the delta array.

        \param[in,out] minScore The current minimum score. This value
        is updated after the score is updated to reflect the minimum
        of score and minScore.
    */
    inline void updateWindow(const int wordIn, const int wordOut,
                             int& score, int& minScore) {
        // Update score and delta for word moving in
        score -= (delta[wordIn] << 1) - 1;
        delta[wordIn]--;
        // Update score and delta for word moving out
        score += (delta[wordOut] << 1) + 1;
        delta[wordOut]++;
        // Track the minimum score.
        minScore = std::min(score, minScore);
    }

    /** Helper method to update the scores based on a sliding window
        and track the index positions of the minimum scores.

        This method is invoked from several different spots from the
        runD2() method to update the d2 scores as the window slides
        across the two sequences being analyzed. The hash values of
        the words moving into and out of the window are used to update
        the scores.

        \param[in] wordIn The hash value of the word moving into the
        window.

        \param[in] wordOut The hash value of the word moving out of
        the window.

        \param[in,out] score The current running score for this
        window. This value is updated using the delta array.

        \param[in,out] minScore The current minimum score. This value
        is updated after the score is updated to reflect the minimum
        of score and minScore.

        \param[in] s1Pos The current index position in the 1st strain.

        \param[in] s2Pos The current index position in the 2nd strain.
        
        \param[out] s1MinPos The index position (in 1st strain) of the
        minimum score thus far.  This value is updated only when a new
        minimum score is encountered.

        \param[out] s1MinPos The index position (in 2nd strain) of the
        minimum score thus far.  This value is updated only when a new
        minimum score is encountered.
    */
    inline void updateWindow(const int wordIn, const int wordOut,
                             int& score, int& minScore, int s1Pos, int s2Pos,
                             int& s1MinPos,
                             int& s2MinPos) {
        // Update score and delta for word moving in
        score -= (delta[wordIn] << 1) - 1;
        delta[wordIn]--;
        // Update score and delta for word moving out
        score += (delta[wordOut] << 1) + 1;
        delta[wordOut]++;
        // Track the minimum score.
        if (score < minScore) {
            minScore = score;
            s1MinPos = s1Pos;
            s2MinPos = s2Pos;
            std::cout << "D2(" << s1MinPos << ", " << s2MinPos
                      << ") = " << minScore << std::endl;
        }
    }
    
    /** The core method that run's the D2 algorithm.

        This method performs the core D2 analysis. This method is
        invoked from the getMetric() method to run the core D2
        algorithm. This method operates as follows:

        <ol>

        <li>If a heuristic chain has been set, then this method
        obtains a hint from the chain to decide if the normal or
        reverse complement analysis must be performed. The default is
        normal analysis.</li>

        <li>Next, a hash table of words is built for the otherEST
        sequence via a call to the buildWordTable method.</li>

        <li>Next, it computes the score using the first two windows.</li>

        <li>Finally, for each window in the reference EST it iterates
        over the windows in the other EST (sliding windows in each
        inner iteration) to compute the minimum d2 score by calling
        the updateWindow() method.</li>

        <li>It finally returns the minimum d2 score recorded.</li>

        </ol>

        \param[in] otherEST A pointer to an immutable EST object to be
        analyzed by this method.
        
        \return The d2 score between the reference EST (set via call
        to setReferenceEST()) and the otherEST (parameter).
    */
    float runD2(const EST* otherEST);

    /** The hint key that is used to add hint for normal or
        reverse-complement D2 computation.
		
        This hint key is used to set a hint in the \c hints hash
        map. This string is defined as a constant to save compute time
        in the core \c runHeuristics method.
    */
    const std::string hintKey;
};

#endif
