#ifndef TV_HEURISTIC_H
#define TV_HEURISTIC_H

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

#include "NewUVHeuristic.h"
#include "ESTList.h"
#include <utility>
#include <vector>

/** Heuristic based upon the T/V heuristic, a type of common word
    heuristic.

    The idea of the <i>t/v</i>-heuristic is to require the common
    words in a pair of ESTs to be found reasonably close to each other
    but not too close. The rule of this heuristic is as follows:

    <ol>

    <li>Assume we are given two ESTs to analyze, say e<sub>i</sub> and
    e<sub>j</sub> and a threshold \em t.</li>

    <li>Consider all \em v words appearing in e<sub>i</sub>.</li>

    <li>At least \em t of these \em v words must appear in \em j so
    that they do not overlap (their starting positions must be at
    least \em v base pairs different) and are at least 100 base pairs
    of each other.</li>

    <li>If there are at least \em t \em v words the huristic
    passes. If not, the pair need not be considered further.</li>
    
    </ol>
*/
class TVHeuristic : public NewUVHeuristic {
    friend class HeuristicFactory;
public:
   /** Add valid command line arguments for this heuristic.
        
        This method is invoked by the SubSystem (that logically owns
        the component) when it is being requested for command line
        arguments.  This method must be used to add all valid command
        line options that are supported by this heuristic.  Note that
        derived classes may override this method to add additional
        command line options that are applicable to it.

        \note Derived heuristic classes may override this method to
        add their custom command line arguments.  When this method is
        overridden <b>don't forget to call the corresponding base
        class implementation</b> to add common options.

        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);
    
    /** Method to begin heuristic analysis (if any).
        
        This method is invoked just before commencement of EST
        analysis.  This method essentially passes control to the base
        class that merely creates the arrays for building hash maps.
        
        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns with a
        non-zero error code.
    */
    virtual bool initialize();
    
    /** Set the reference EST id for analysis.
        
        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  Setting the
        reference EST provides this heuristic an opportunity to
        pre-compute the normal and reverse-complement hash tabes for
        words of varying sizes. The hash table enables rapid searching
        for words in the \c runHeuristic method.

        \note This method must be called only after the initialize()
        method is called.

		\param[in] refEST The reference EST against which a large
		number of other EST entries are going to be compared via
		subsequent calls to the shoudAnalyze(const EST*) method.
        
        \return If the hash table creation process was sucessful, then
        this method returns 0.  Otherwise this method returns an error
        code.
    */
    virtual int setReferenceEST(const EST* refEST);

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
    int getWindowLen() { return windowLen; }
    
protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead it should be
        created via a suitable call to the HeuristicFactory API
        method(s).

        \param[in] chain The heuristic chain that logically contains
        this heuristic.  This value is used by derived classes to
        update hints (used by heavy weight analyzers) in the heuristic
        chain class.  This value is saved in a pointer in the base
        class (Heuristic::heuristicChain) and is never changed or
        deleted by this class (or its children).
    */
    TVHeuristic(HeuristicChain *chain);
    
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
        
        \param[in] otherEST A pointer to the EST with which the
        reference EST is to be compared.
        
        \return This method returns \c true if the heuristic says the
        EST pair should be analyzed, and \c false if it should not.
    */
    virtual bool runHeuristic(const EST* refEST);

    /** Method to obtain and update the parameters for the heuristic
		based on the parameter set manager.

        \param[in] otherEST A pointer to the other EST with which the
        reference EST is to be compared.
		
		\return true if we should analyze these ESTs, false otherwise
    */
    bool updateParameters(const EST* otherEST);

    /** Templatized-method for counting common woards between two ESTs.

        This method is a helper method that is invoked from the
        runHeuristic method to count the number of common words
        between the reference EST (set via call to setReferenceEST)
        and otherEST (parameter). This method operates as follows:

        <ol>

        <li>First the matchTable (instance variable) is cleared to all
        zeros.</li>

        <li>Next the initial word of length NewUVHeuristic::v is
        constructed while ignoring bases marked as 'n' (this may
        require processing of more than the first NewUVHeuristic::v
        bases if one of them is a 'n'.</li>
        
        </ol>

     */
    template <typename Encoder>
    int countCommonWords(const EST* otherEST, Encoder encoder,
                         const char* refWordMap) {
        // First compute the hash for the first word using the supplied
        // encoder.
        const char *otherSeq   = otherEST->getSequence();
		const int   otherESTLen= otherEST->getSequenceLength();
        register int hash      = 0;
        int ignoreMask         = 0;
        // First, reset the array that tracks matching counts as the
        // window slides across the otherEST.
        memset(&matchTable[0], 0, sizeof(char) * matchTable.size());
        // Compute hash for initial word while skipping over bases
        // makred 'n'. This may require processing of more then v-1
        // bases
        for(int i = 0; (i < NewUVHeuristic::v - 1); i++) {
            hash = encoder(hash, otherSeq[i], ignoreMask);
        }
        // Skip first windowLen entries to simplify logic in loop below.
        char *matchTicker = &matchTable[0] + windowLen;
        // Now see how many common words exist in the two ESTs
        int numMatch     = 0, maxMatch = 0;
		int oldWindowPos = -windowLen;
        for(int i = NewUVHeuristic::v - 1; (i < otherESTLen); i++) {
            hash = encoder(hash, otherSeq[i], ignoreMask);
	    numMatch -= matchTicker[oldWindowPos++];
            if (!ignoreMask) {
                // If ignoreMask is zero, it tells us that the hash is
                // not tainted due to 'n' bp and it is good to be
                // used for the operations below.
                matchTicker[i] = refWordMap[hash];
                numMatch += refWordMap[hash];
                maxMatch  = std::max(maxMatch, numMatch);
            }
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
    /** The number of minumum number of common words.

        This instance variable contains the minimum number of words
        (that are close but not too close) that have matching values
        in pairs of ESTs.  The default is 65. However, this value can
        be overridden by a command line argument.
    */
    int t;

    /** The window length to be used for <i>t/v</i> heuristic.

        The window length defines the length of the window within
        which common words are tracked and reported by this heuristic.
        Typically, this window length must match the window length
        used for D2 analysis for the heuristic to be meanigful. The
        default value is 100.  This value can be overridden by the
        user via suitable command line arguments.
    */
    int windowLen;
    
    /** A large table to track matches.

        This instance variable contains a large table that tracks
        number matches encountered for each word-hash, as this
        heuristic tracks matching words.  Earlier this was an array of
        characters which made it a bit cumbersome to
        track/troubleshoot issues.  Hence, it has been convered to a
        vector<char>.
    */
    std::vector<char> matchTable;

    /** Instance variable to track the number of times UV sample
        heuristic passed.

        This instance variable is used to track the number of times
        the UV sample heuristic passed.  This value indicates the
        number of times the TV heuristic was actually run.  This value
        is incremented in the runHeuristic method and is displayed by
        the printStats() method.
    */
    int uvSuccessCount;

	/** Index of the current parameter set being used for analysis.

		This instance variable is used to track the index of the
		current parameter set that is being used for analysis.  This
		instance variable is initialized to -1.  The
		updateParameters() method sets this value to reflect the
		current parameter set being used.  The parameters are updated
		only when a new parameter set needs to be used for analysis.
	*/
	int currParamSetIndex;
};

#endif
