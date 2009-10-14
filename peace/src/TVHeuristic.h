#ifndef TV_HEURISTIC_H
#define TV_HEURISTIC_H

//---------------------------------------------------------------------------
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//          James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "NewUVHeuristic.h"
#include "EST.h"
#include <utility>

/** Heuristic based upon the T/V heuristic used in WCD, a type of
    common word heuristic.

    The idea of the <i>t/v</i>-heuristic is to require the common
    words in a pair of ESTs to be found reasonably close to each other
    but not too close. The rule of this heuristic is as follows:

    <ol>

    <li>Assume we are given two ESTs to analyze, say e<sub>i</sub> and
    e<sub>j</sub> and a threshold \i t.</li>

    <li>Consider all \i v words appearing in e<sub>i</sub>.</li>

    <li>At least \i t of these \i v words must appear in \i j so that
    they do not overlap (their starting positions must be at least \i
    v base pairs different) and are at least 100 base pairs of each
    other.</li>

    <li>If there are at least \i t \i v words the huristic passes. If
    not, the pair need not be considered further.</li>
    
    </ol>
*/
class TVHeuristic : public NewUVHeuristic {
    friend class HeuristicFactory;
public:
    /** Display valid command line arguments for this heuristic.
        
        This method is used to display all valid command line options
        that are supported by this heuristic.  Note that this method
        invokes the corresponding method in the base class to display
        any options supported by the base class.  This method is
        typically used in the \c main() method when displaying usage
        information.

        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.
        
        This method is used to process command line arguments specific
        to this heuristic.  This method is typically used from the \c
        main method just after the heuristic has been instantiated.
        This method consumes all valid command line arguments.  If the
        command line arguments were valid and successfully processed,
        then this method returns \c true.

        \note This method calls the corresponding base class
        implementation to display common options.
        
        \param[inout] argc The number of command line arguments to be
        processed.  When arguments are processed/consumed this
        parameter is changed to reflect the number of arguments
        processed.
        
        \param[inout] argc The array of command line arguments to be
        processed. When arguments are processed/consumed, the consumed
        arguments are removed from this array.
        
        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.
    */
    virtual bool parseArguments(int& argc, char **argv);
    
    /** Method to begin heuristic analysis (if any).
        
        This method is invoked just before commencement of EST
        analysis.  This method essentially passes control to the base
        class that merely creates the arrays for building hash maps.
        
        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns with a
        non-zero error code.
    */
    virtual int initialize();
    
    /** Set the reference EST id for analysis.
        
        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  Setting the
        reference EST provides this heuristic an opportunity to
        pre-compute the normal and reverse-complement hash tabes for
        words of varying sizes. The hash table enables rapid searching
        for words in the \c runHeuristic method.

        \note This method must be called only after the initialize()
        method is called.

        \return If the hash table creation process was sucessful, then
        this method returns 0.  Otherwise this method returns an error
        code.
    */
    virtual int setReferenceEST(const int estIdx);

    /** The destructor.

        The destructor frees memory allocated for holding any dynamic
        data.  in the base class.
    */
    virtual ~TVHeuristic();

    /** Obtain the window length used for <i>t/v</i> heuristic.

        The window length defines the length of the window within
        which common words are tracked and reported by this heuristic.
        Typically, this window length must match the window length
        used for D2 analysis for the heuristic to be meanigful. The
        default value is 100.  This value can be overridden by the
        user via suitable command line arguments.

        \return The current window (or frame) size set for <i>t/v</i>
        heuristic.
    */
    static int getWindowLen() { return windowLen; }
    
protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead it should be
        created via a suitable call to the HeuristicFactory API
        method(s).

        \param[in] outputFileName The output file to which any
        heuristic data is to be written. Currently, this value is
        ignored.
    */
    TVHeuristic(const std::string& outputFileName);
    
    /** Determine whether the analyzer should analyze, according to
	this heuristic.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method. This method operates as follows:

        <ol>

        <li>It invokes the corresponding method in the base class to
        first run the UV-sample heuristic on the pair of ESTs. If the
        pair fails UV-sample heuristic this method returns immediately
        with \c false (indicating further analysis is not needed).</li>

        <li>If the pair passes UV-sample heuristic then this method
        invokes the overloaded \c runHeuristic method with a suitable
        encoder (normal or reverse-complement encoder depending on the
        value of NewUVHeuristic::bestMatchIsRC flag) to analyze
        the pair of ESTs using the TV heuristic.</li>
        
        </ol>
        
        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.
        
        \return This method returns \c true if the heuristic says the
	EST pair should be analyzed, and \c false if it should not.
    */
    virtual bool runHeuristic(const int otherEST);

    template <typename Encoder>
    int countCommonWords(const int otherEST, Encoder encoder,
                         const char* refWordMap) {
        // First compute the hash for the first word using the supplied
        // encoder.
        const char *otherSeq   = EST::getEST(otherEST)->getSequence();
        const int  otherESTLen = strlen(otherSeq);
        register int hash      = 0;
        for(int i = 0; (i < NewUVHeuristic::v - 1); i++) {
            hash = encoder(hash, otherSeq[i]);
        }
        // Set first window length entries to zero.
        memset(matchTable, 0, sizeof(char) *( windowLen + v));
        // Skip first windowLen entries to simplify logic in loop below.
        char *matchTicker = matchTable + windowLen;
        // Now see how many common words exist in the two ESTs
        int numMatch     = 0, maxMatch = 0;
	int oldWindowPos = -windowLen;
        for(int i = NewUVHeuristic::v - 1; (i < otherESTLen); i++) {
            hash = encoder(hash, otherSeq[i]);
            matchTicker[i] = refWordMap[hash];
            numMatch += refWordMap[hash];
            numMatch -= matchTicker[oldWindowPos++];
            maxMatch  = std::max(maxMatch, numMatch);
        }
        // Return number of matches encountered.
        return maxMatch;
    }

    /** Method to display statistics regarding operation of this
        heuristic.
            
        This method can be used to obtain a dump of the statistics
        gathered regarding the operation of this heuristic. This
        method calls the base class method first which prints some
        common statistics.  It then displays the number of times the
        <i>u/v</i> sample heuristic (the base class) returned success
        causing the <i>t/v</i> heuristic to be run.
        
        \param[out] os The output stream to which the statistics
        regarding the heuristic is to be dumped.
    */
    virtual void printStats(std::ostream& os) const;
    
private:
    /** The set of arguments specific to the TV heuristic.

        This instance variable contains a static list of arguments
        that are specific only to this analyzer class.  This argument
        list is statically defined and shared by all instances of this
        class.

        \note Use of static arguments and parameters renders this UV
        sample heuristic class not to be MT-safe.
    */
    static arg_parser::arg_record argsList[];

    /** The number of minumum number of common words.

        This instance variable contains the minimum number of words
        (that are close but not too close) that have matching values
        in pairs of ESTs.  The default is 65. However, this value can
        be overridden by a command line argument.
    */
    static int t;

    /** The window length to be used for <i>t/v</i> heuristic.

        The window length defines the length of the window within
        which common words are tracked and reported by this heuristic.
        Typically, this window length must match the window length
        used for D2 analysis for the heuristic to be meanigful. The
        default value is 100.  This value can be overridden by the
        user via suitable command line arguments.
    */
    static int windowLen;
    
    /** A large table to track matches.

        This instance variable contains a large table that tracks
        matches encountered as this heuristic tracks matching words.
    */
    char *matchTable;

    /** Instance variable to track the number of times UV sample
        heuristic passed.

        This instance variable is used to track the number of times
        the UV sample heuristic passed.  This value indicates the
        number of times the TV heuristic was actually run.  This value
        is incremented in the runHeuristic method and is displayed by
        the printStats() method.
    */
    int uvSuccessCount;
};

#endif
