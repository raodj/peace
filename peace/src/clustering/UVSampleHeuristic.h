#ifndef UV_SAMPLE_HEURISTIC_H
#define UV_SAMPLE_HEURISTIC_H

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

#include "arg_parser.h"
#include "Heuristic.h"
#include "HeuristicFactory.h"
#include "Utilities.h"
#include "HashMap.h"

#include <string>
#include <vector>

/** \typedef HashMap<int, std::vector<short> > UVHashTable;

    \brief A hash map to provide rapid access to the set of
    pre-computed hash values for a given EST.

    This typedef provides a convenient short cut to refer to a hash
    map that is used to cache pre-computed hash values for a given
    EST.  The zero-based index of the EST is used as the key into this
    hash map.  Each entry in the hash map contains a std::vector that
    essentially contains the array of hash values.

    \note Refer to the additional documentation with the instance
    variable uvCache for motivation for this hash map.
*/
typedef HashMap<int, std::vector<unsigned short> > UVHashTable;

/** Heuristic based upon the "u/v sample heuristic" used in WCD,
    a type of common word heuristic.  Considers all words of length v
    in the first sequence and every 16th word of length v in the second
    sequence.  Returns true if it finds at least u common words.
*/
class UVSampleHeuristic : public Heuristic {
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
        data in the base class.
    */
    virtual ~UVSampleHeuristic();

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
    UVSampleHeuristic(HeuristicChain* chain,
					  const std::string& name = "uv");
    
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
    
    /** Instance variable to track if a given word (of length \c v)
        appears in reference EST.

        This instance variable is created in the initialize() method
        to point to an array of 4<sup>v</sup>

    */
    char* s1WordMap;

    /** Instance variable to track if a given word (of length \c v)
        appears in reference EST.

        This instance variable is created in the initialize() method
        to point to an array of 4<sup>v</sup> 
    */
    char* s1RCWordMap;

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

    /** Instance variable to maintain the \em v parameter for the
		<i>u/v</i> heuristic.

        This instance variable contains the value specified by the
        user for the command line argument \c --uv_v.  This value
        indicates the length of each word to be considered in the two
        sequences being compared.  For example if the value of v is
        16, then it compares two 16-bp sub-sequences to determine if
        they have at least UVSampleHeuristic::u common words.  The
        default value is 8.
    */
    int v;

    /** Number of nucleotide bases to skip when searching for <i>v</i>
        words.
        
        This instance variable indicates the number of nucleotide
        bases to skip when comparing two sequences.  Skipping base
        pairs improves performance but impacts the quality of results.
        Consequently, there is a tradeoff between
        UVSampleHeuristic::v, UVSampleHeuristic::u, and this instance
        variable.  Larger this value faster is the heuristic but its
        sensitivity/quality will degrade (and vice versa).  The
        default value for this variable is 16.  However, the value can
        be changed by the user via the \c --uv_wordShift command-line
        parameter.
    */
    int wordShift;

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
    
private:
    /** Number of <i>v</i>-word matches that signify sufficiently
        similar fragments.

        This instance variable is used to track the number of
        sufficiently similar <i>v</i>-word matches to be expected
        above which two fragments are declared as being sufficiently
        similar by this heuristic.  The default value for this
        instance variable is 4.  However, the user can change this
        value via the command-line parameter \c --uv_u.  Bear in mind
        that there is a symbiotic relationship between this value and
        UVSampleHeuristic::v and UVSampleHeuristic::wordShift.
    */
    int u;

    /** A hash map to cache hash values (<i>v</i> base pairs in
        length) to seedup <i>u/v</i> heuristic.

        <p>The <i>u/v</i> heuristic \b used to generate hash values
        (in the \c runHeuristic method) by iterating over the base
        pairs in a given EST sequence.  However, this approach turned
        out to be rather slow. Furthermore, it was observed that the
        same set of hash values were recomputed for various pair-wise
        EST comparisons.</p>
        
        <p>Therefore, to improve the overall performance of the
        <i>u/v</i> heuristic, it was proposed that the hash values be
        cached to improve performance. Of course, this does increase
        the net amount of memory consumed. Consequently, to aid in
        caching only the required subset of EST hash values, this hash
        map was introduced to store the necessary sequences and
        rapidly access them when needed.</p>

        <p>The entries in this hash map are computed by the \c
        computeHash method, which is invoked from the runHeuristic
        method.</p>
    */
    UVHashTable uvCache;
};
 
#endif
