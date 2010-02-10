#ifndef TWOPASS_D2_H
#define TWOPASS_D2_H

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
#include <string>
#include <vector>

// Forward declaration to keep compiler happy
class EST;
class ResultLog;

/** \brief EST Analyzer that uses the D2 algorithm in both its
	symmetric and asymmetric variations to compute distance between
	two ESTs.

    <p>This analyzer provides the mechanism to use D2 algorithm to
    compute the distance values between a pair of ESTs. For improved
	performance, the analyzer uses D2 asymmetric and symmetric together.
	It runs a fast "first pass" using D2 asymmetric and finds the minimum
	D2 score, as well as the two windows between which the score was
	computed.  The analyzer then computes bounds and runs a bounded
	version of D2 symmetric.  Ideally this analyzer will provide
	faster runtime because of the asymmetric pass, and not sacrifice any
	precision because of the bounded symmetric pass.<p>

    \note This D2 analyzer uses a word size of 8 base pairs.  This is
    hard coded into the algorithm in the form of various assumptions
    (data types of variables etc.). A similar approach is used by
    other D2 implementations as well.  This is a compromise between
    performance and flexibility of the implementation.
*/
class TwoPassD2 : public FWAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~TwoPassD2();
    
    /** Display valid command line arguments for this analyzer.
	
        This method must be used to display all valid command line
        options that are supported by this analyzer.  Currently, this
        analyzer does not require any special command line parameters.

        \note The ESTAnalyzer base class requires that derived EST
        analyzer classes <b>must</b> override this method to display
        help for their custom command line arguments.  When this
        method is overridden don't forget to call the corresponding
        base class implementation to display common options.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.

        <p> This method is used to process command line arguments
        specific to this EST analyzer.  This method is typically used
        from the \c main method just after the EST analyzer has been
        instantiated.  This method consumes all valid command line
        arguments.  If the command line arguments were valid and
        successfully processed, then this method returns \c true. </p>

        <p>Currently, this EST analyzer accepts an optional \c
        frameShift value to control the stride used for D2
        algorithm. However, the default value of 1 for \c frameShift
        is the preferred value for this parameter.</p>

        \note The \c ESTAnalyzer base class requires that derived EST
        analyzer classes <b>must</b> override this method to process
        any command line arguments that are custom to their operation.
        However, as per API requirement, this methodcalls the
        corresponding base class implementation to display additional
        options.
        
        \param[in,out] argc The number of command line arguments to be
        processed.

        \param[in,out] argv The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.  This method checks to ensure that a valid
        frame size and a valid word size have been specified.
    */
    virtual bool parseArguments(int& argc, char **argv);
    
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
    int initialize();

    /** Method to obtain human-readable name for this EST analyzer

        This method provides a human-readable string identifying the
        EST analyzer.  This string is typically used for
        display/debugging purposes (particularly via the PEACE
        Interactive Console).

        \return This method returns the string "layeredD2" identifiying 
        this analyzer.
    */    
    virtual std::string getName() const { return "layeredD2"; }
    
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

        \return This method returns \c true if the estIdx was within
        the given range of values.  Otherwise this method returns a
        non-zero value as the error code.
    */
    virtual int setReferenceEST(const int estIdx);

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

        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.

        \return This method returns the distance value reported by the
        D2 algorithm.
    */
    virtual float getMetric(const int otherEST);
    
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
		
        \param[out] wordTable The word table to be populated with with
        hash values.
		
        \param[in] estSeq The sequence of base pairs associated with
        this EST that must be converted to hash values.

		\param encoder The encoder object to be used to generate
		encoding for the words added to the generated word table.
    */
    template <typename Encoder>
    void buildWordTable(int* wordTable, const char* estSeq, Encoder encoder) {
        // Compute the has for the first word using the encoder. The
        // encoder may do normal or reverse-complement encoding for us.
        register int hash = 0;
        for(register int i = 0; ((*estSeq != 0) && (i < wordSize - 1));
            i++, estSeq++) {
            hash = encoder(hash, *estSeq);
        }
        // Now compute the hash for each word
        for(int entry = 0; (*estSeq != 0); entry++, estSeq++) {
            hash = encoder(hash, *estSeq);
            wordTable[entry] = hash;
        }
    }

    /** Instance variable that maps an index in the reference sequence
        (sequence s1) to the hash of a word.  This hash can then be used
        as an index in the fdHashMap to get the frequency differential
        for that word.

        The word table is created in the initialize() method and
        filled in using the buildWordTable() method.  For the reference
        EST, buildWordTable() is called in the setReferenceEST() method,
        meaning it does not need to be rebuilt every time we analyze
        a new comparison sequence.
    */
    int* s1WordTable;

    /** Instance variable that maps an index in the comparison sequence
        (sequence s2) to the hash of a word.

        This hash is used as an index in the fdHashMap to get the
        frequency differential for that word.

        The word table is created in the initialize() method and
        filled in using the buildWordTable() method.  For the comparison
        EST, buildWordTable() must be called in the analyze() method because
        a new comparison sequence is given every time analyze() is called.
    */
    int* s2WordTable;

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
    static int frameShift;
    
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
        the "--threshold" command line argument.
    */
    static int threshold;
	
    /** The threshold score above which two ESTs are considered
        sufficiently dissimilar that a more accurate symmetric
        D2 pass would still not cluster the ESTs.
		
        This instance variable tracks the threshold value to be used
        to call the runD2Bounded method.  If the calculated D2 score
        from the runD2Asymmetric method is above this threshold,
        bounded D2 will not be run.  Currently, the default value is
        130.  This value can be overridden by the user with the
        "--maxThreshold" command line argument.
    */
    static int maxThreshold;
    
private:
    /** The set of arguments specific to the D2 algorithm

        This instance variable contains a static list of arguments
        that are specific only to the D2 analyzer class.  This
        argument list is statically defined and shared by all
        instances of this class.

        \note Use of static arguments and parameters makes D2 class
        hierarchy not MT-safe.
    */
    static arg_parser::arg_record argsList[];
	
    /* The default constructor for this class.
       
       The default constructor for this class.  The constructor is
       made private so that this class cannot be directly
       instantiated.  However, since the ESTAnalyzerFactory is a
       friend of this class; therefore it can instantiate the D2
       analyzer.  Accordingly, the ESTAnalyzerFactory::create()
       method must be used to instantiate this class.

       \param[in] refESTidx The reference EST index value to be used
       when performing EST analysis.  This parameter should be >= 0.
       This value is simply passed onto the base class.
        
       \param[in] outputFile The name of the output file to which the
       EST analysis data is to be written.  This parameter is ignored
       if this analyzer is used for clustering.  If this parameter is
       the empty string then output is written to standard output.
       This value is simply passed onto the base class.
    */
    TwoPassD2(const int refESTidx, const std::string& outputFileName);

    /** Instance variable to track alignment metric computed by the
        analyze() method.
        
        This instance variable is used to hold the alignment metric
        that was computed in the previous \c analyze method call. By
        default this value is set to zero.  The alignment metric is
        computed as the difference in the window positions (on the two
        ESTs being analyzed) with the minimum d2 distance.
    */
    int alignmentMetric;

    /** The bitmak to be used when build hash values.
            
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

    /** Helper method to update the delta table as well as the
        minimum score (if it has changed) when the window is shifted.
		
        For d2 asymmetric, the updateWindow method also has to keep
        track of the location on each sequence where the minimum
        score was found, which is the function of the variables
        s1MinScoreIdx and s2MinScoreIdx.  The variables s1CurWindowIdx
        and s2CurWindowIdx mark the start of the current window on each 
        sequence.  All of these variables are numbers, corresponding to
        indices of the word tables.
    */
    inline void updateWindowAsym(const int wordIn, const int wordOut,
                                 int& score, int& minScore,
                                 const int s1CurWindowIdx,
                                 const int s2CurWindowIdx,
                                 int* s1MinScoreIdx, int* s2MinScoreIdx) {
        // Update score and delta for word moving in
        score -= (delta[wordIn] << 1) - 1;
        delta[wordIn]--;
        // Update score and delta for word moving out
        score += (delta[wordOut] << 1) + 1;
        delta[wordOut]++;
        // Track the minimum score.
        if (score < minScore) {
            minScore = score;
            // Update the indices for the new min score
            *s1MinScoreIdx = s1CurWindowIdx;
            *s2MinScoreIdx = s2CurWindowIdx;
        }
    }
    
    /** Helper method to update the delta table as well as the
        minimum score (if it has changed) when the window is shifted.
        This method is identical to the updateWindow method found in
        D2.h (Peace's base D2 implementation).
    */
    inline void updateWindow(const int wordIn, const int wordOut,
                             int& score, int& minScore,
                             const int windowDistance) {
        // Update score and delta for word moving in
        score -= (delta[wordIn] << 1) - 1;
        delta[wordIn]--;
        // Update score and delta for word moving out
        score += (delta[wordOut] << 1) + 1;
        delta[wordOut]++;
        // Track the minimum score.
        if (score < minScore) {
            minScore = score;
            alignmentMetric = windowDistance;
        }
    }

    /** Method to run the D2 algorithm, asymmetric version.
        This method runs D2 asymmetric on the two EST sequences
        and finds the minimum D2 score.  Since this is D2 asymmetric,
        the score is not guaranteed to be the true minimum.  Thus,
        this method also finds the indices of the windows on each
        sequence where the minimum D2 score was found.  The analyzer
        will use these to compute appropriate bounds for the bounded
        symmetric D2 function.
    */
    float runD2Asymmetric(const int otherEST, int* s1MinScoreIdx,
                          int* s2MinScoreIdx); 	

    /** Method to run the D2 algorithm, symmetric version.
        This method is essentially identical to the runD2 method
        in the D2.cpp class (Peace's base D2 implementation) with
        the exception that this method places bounds on each sequence.
        It will not search beyond these bounds.  This saves a great
        deal of computation time (though the asymptotic runtime
        is unchanged) by avoiding unnecessary computations and only
        searching in the area where the best match might be found.
    */
    float runD2Bounded(const int otherEST, int sq1Start, int sq1End, 
                       int sq2Start, int sq2End);

    /** The hint key that is used to add hint for normal or
        reverse-complement D2 computation.

        This hint key is used to set a hint in the \c hints hash
        map. This string is defined as a constant to save compute time
        in the core \c runHeuristics method.
    */
    const std::string hintKey;
};

#endif
