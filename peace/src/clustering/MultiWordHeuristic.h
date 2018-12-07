#ifndef MULTI_WORD_HEURISTIC_H
#define MULTI_WORD_HEURISTIC_H

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

#include "Heuristic.h"
#include "HeuristicFactory.h"
#include "Utilities.h"
#include "HashMap.h"

#include <string>
#include <vector>

/** \typedef HashMap<int, bool> ESTMap;

    \brief A hash map to rapidly determine if a given EST (identified
    by its integer ID) is present in a given Map.

    This typedef provides a convenient short cut to refer to a hash
    map that is used to track if a given EST has a given word
    (identified via a hash value).  The zero-based index of the EST is
    used as the key into this hash map.
*/
typedef HashMap<int, bool> ESTMap;

/** \typedef HashMap<int, ESTMap> WordMap;

    \brief A hash map to rapidly add a given EST to the list of
    entries associated with a given hash word.

    This typedef provides a convenient short cut to refer to a hash
    map that is used to track if a given EST has a given word
    (identified via a hash value).  The hash of a given word is used
    as the key into this hash map.
*/
typedef HashMap<int, ESTMap> WordMap;

/** A Multi-Word (MW) Heuristic based upon the number of matching
    h-words between two sequences.
*/
class MultiWordHeuristic : public Heuristic {
    friend class HeuristicFactory;
public:
    /** Add valid command line arguments for this heuristic.
        
        This method is invoked by the SubSystem (that logically owns
        the component) when it is being requested for command line
        arguments.  As per API, this method is used to add the
        following valid command line options that are supported by
        this heuristic:

	<ul>

	<li>\c --mw-wordLen: This command-line argument can be used to
	change the length of the word used for finding matches. The
	default value is 8 nucleotides. </li>

	<li>\c --mw-wordCount: The number of words (that are \c
	wordLen each) that must match in order to report two ESTs are
	sufficiently similar.  The default value is 4.</li>
	
	</ul>
	
        \note Derived heuristic classes may override this method to
        add their custom command line arguments.  When this method is
        overridden <b>don't forget to call the corresponding base
        class implementation</b> to add common options.
	
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);
    
    /** Method to begin heuristic analysis (if any).

        This method is invoked when the SubSystem (that logically owns
        this component) is being initialized.  This method is invoked
        just before commencement of EST analysis.  This method
        typically loads additional information that may be necessary
        for a given heuristic from data files.  In addition, it may
        perform any pre-processing as the case may be.

        \return If the initialization process was sucessful, then this
        method returns \c true.  Otherwise this method returns \c
        false to signify error (errors are reported on std::cerr).
    */
    virtual bool initialize();
    
    /** Set the reference EST for analysis.
        
        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  Setting the
        reference EST provides heuristics an opportunity to optimize
        certain operations, if possible.

        \note This method must be called only after the initialize()
        method is called.

	\param[in] est A pointer to an immutable EST object to be used
	as the reference.  The reference EST will be subsequently
	analyzed with a number of other ESTs via the analyze() method.

        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns an error code.
    */
    virtual int setReferenceEST(const EST* est);

    /** The destructor.

        The destructor frees memory allocated for holding any dynamic
        data in this class.
    */
    virtual ~MultiWordHeuristic();

    /** The core method that computes matching words for all pairs of
	reads in a distributed manner.

	This method is a convinence method that has been introduced
	here to facilitate the process of computing matching words for
	all pairs of reads without comparing all-pairs.  This method
	is run on all the parallel processes being used and operates
	in a distributed manner as follows:


	<ol>

	<li>It assumes that the MultiWordHeuristic::initialize() has
	been called as per the PEACE framework work flow.</li>
	
	<li>Next it determines the sub-set of reads that this process
	is expected to operate.  For each read, it computes the words
	and adds the read-number to the list of reads associated with
	each word building the reverse-look-up data.</li>
	
	<li>It participates in interative broadcast in which it
	receives filter data from other processes (if any) and
	broadcasts its own data to others. This process ensures that
	all processes have a consistent view of the entries that have
	been filtered out on other processes.</li>

	<li>If all the operations were successfully completed, then
	this method returns 0 (zero).</li>
		
	</ol>

        \note This method assumes that the runtime context has been
        setup for the ClusteringSubSystem (as per the normal runtime
        operations).  The context is used by this heuristic class to
        obtain pointers to the ESTAnalyzer and ClusterMaker objects
        for its use.
        
	\return This method returns 0 (zero) on success. On errors it
	returns a non-zero error code.
    */
    virtual int run();
    
protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived Heuristic classes must be instantiated via the
        HeuristicFactory API methods.

        \param[in] chain The heuristic chain that logically contains
        this heuristic.  This value is used by derived classes to
        update hints (used by heavy weight analyzers) in the heuristic
        chain class.  This value is saved in a pointer in the base
        class (Heuristic::heuristicChain) and is never changed or
        deleted by this class (or its children).
		
        \param[in] name The name (specified by final derived child
        class) to be set for this heuristic.  This value is typically
        used for reporting various types of messages.
    */
    MultiWordHeuristic(HeuristicChain* chain,
		       const std::string& name = "mw");
    
    /** Determine whether the analyzer should analyze, according to
        this heuristic.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

	\param[in] otherEST A pointer to an immutable EST with which
        the reference EST is to be compared.
        
        \return This method returns true if the heuristic says the
        EST pair should be analyzed, and false if it should not.
    */
    virtual bool runHeuristic(const EST* otherEST);

    /** Method to compute the <i>u/v</i> hash values for a given EST.

        This method is a utility method that was introduced to
        streamline the process of computing and caching <i>u/v</i>
        hash values for a given EST.  This method uses the variables
        \c u, \c v, and \c wordShift (all of them user configurable)
        along with a ESTCodec::NormalEncoder object to compute hash
        values into a std::vector. The vector is added to the \c
        uvCache hash map for future reference.

        \param[in] est A pointer to an immutable EST object whose hash
        values is to be computed and cached.
    */
    void computeHash(const EST* est);
    
    /** Flag to indicate if normal or reverse-complement version
        provided best match.

        This flag is set at the end of the runHeuristic method in this
        class to indicate if the normal or the reverse-complement
        check yielded the best possible match. If this flag is \c
        false, then the normal check yielded the best match. If this
        value is \c true, then the reverse-complement check yielded
        the best match.
    */
    bool bestMatchIsRC;

    /** A convenience instance variable that is initialized once and
        used to call templatized helpers.

        When using templatized-helpers it seems to be convenient to
        have an instance variable to pass to them.  Since this is
        value is passed in a template parameter, it is defined to be
        static (to ensure that it has external linkage as per the
        ISO/ANSI standards requirement).  This instance variable is
        set in the initialize() method and used when running the
        heuristic.
    */
    static int BitMask;

    /** Instance variable to store the number of bits to be shifted to
        create hash values.

        <p>This instance variable is set to the value of 2 * (\em v -
        1) (in the \c initialize method) to reflect the number of bits
        that need to be shifted in order to build the hash values for
        common words (including the values stored in \c s1WordMap and
        \c s1RCWordMap).</p>

        <p>This instance variable is actually passed on to the
        ESTCodec::NormalEncoder or ESTCodec::RevCompEncoder when
        computing hash values.  Since this is value is passed in a
        template parameter, it is defined to be static (to ensure that
        it has external linkage as per the ISO/ANSI standards
        requirement).</p>
    */
    static int bitsToShift;

    /** Update the word map with the words occurring in the given EST.
        
        This method is invoked from the MultiWordHeuristic::run()
        method to update the occurrences of the words in the wordMap.

        \note This method must be called only after the initialize()
        method is called.

	\param[in] est A pointer to an immutable EST object for which
	the words occurring in its sequence are to be updated in the
	wordMap.
    */
    void updateWordMap(const EST* est);
    
private:
    /** Instance variable to maintain the \em wordLen parameter for
	the <i>MW</i> heuristic.

        This instance variable contains the value specified by the
        user for the command line argument \c --mw_wordLen.  This
        value indicates the length of each word to be considered in
        the two sequences being compared.  For example if the value of
        wordLen is 16, then it compares two 16-nucleotide
        sub-sequences to determine if they have at least wordCount
        common words.  The default value for this parameter is 8.
    */
    int wordLen;
        
    /** Number of word matches that signify sufficiently
        similar fragments.

        This instance variable is used to track the number of
        sufficiently similar <i>wordLen</i>-word matches to be
        expected above which two fragments are declared as being
        sufficiently similar by this heuristic.  The default value for
        this instance variable is 4.  However, the user can change
        this value via the command-line parameter \c --mw_wordCount.
        Bear in mind that there is a strong relationship between this
        value and --mw_wordLen.
    */
    int wordCount;

    /** The map 

     */
    WordMap wordMap;
};
 
#endif
