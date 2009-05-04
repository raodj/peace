#ifndef UV_SAMPLE_HEURISTIC_H
#define UV_SAMPLE_HEURISTIC_H

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

#include "arg_parser.h"
#include "Heuristic.h"
#include "HeuristicFactory.h"
#include "Utilities.h"
#include <string>
#include <vector>
#include <HashMap.h>

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
    /** Display valid command line arguments for this heuristic.
        
        This method must be used to display all valid command line
        options that are supported by this heuristic.  Note that
        derived classes may override this method to display additional
        command line options that are applicable to it.  This method
        is typically used in the main() method when displaying usage
        information.

        \note Derived heuristic classes <b>must</b> override this
        method to display help for their custom command line
        arguments.  When this method is overridden don't forget to
        call the corresponding base class implementation to display
        common options.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.
        
        This method is used to process command line arguments specific
        to this heuristic.  This method is typically used from the
        main method just after the heuristic has been instantiated.
        This method consumes all valid command line arguments.  If the
        command line arguments were valid and successfully processed,
        then this method returns \c true.

        \note Derived heuristic classes <b>must</b> override this
        method to process any command line arguments that are custom
        to their operation.  When this method is overridden don't
        forget to call the corresponding base class implementation to
        display common options.
        
        \param[inout] argc The number of command line arguments to be
        processed.
        
        \param[inout] argc The array of command line arguments.
        
        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.
    */
    virtual bool parseArguments(int& argc, char **argv);
    
    /** Method to begin heuristic analysis (if any).
        
        This method is invoked just before commencement of EST
        analysis.  This method typically loads additional information
        that may be necessary for a given heuristic from data files.
        In addition, it may perform any pre-processing as the case may
        be.

        \note Derived classes must override this method.
        
        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns with a
        non-zero error code.
    */
    virtual int initialize();
    
    /** Set the reference EST id for analysis.
        
        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  Setting the
        reference EST provides heuristics an opportunity to optimize
        certain operations, if possible.

        \note This method must be called only after the initialize()
        method is called.

        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns an error code.
    */
    virtual int setReferenceEST(const int estIdx);

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

        \param[in] name The human readable name for this heuristic.
        This name is used when generating errors, warnings, and other
        output messages for this heuristic.

	\param[in] refESTidx The reference EST's index in a given
	multi-FASTA file.  Index values start with 0 (zero).  The
	refESTidx is supplied as a global argument that is processed
	in the main() method.  This value is simply copied to the
	refESTidx member in this class.
    */
    UVSampleHeuristic(const std::string& name, const int refESTidx,
		      const std::string& outputFileName);
    
    /** Determine whether the analyzer should analyze, according to
	this heuristic.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.
        
        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.
        
        \return This method returns true if the heuristic says the
	EST pair should be analyzed, and false if it should not.
    */
    virtual bool runHeuristic(const int otherEST);

    /** Method to compute the <i>u/v</i> hash values for a given EST.

        This method is a utility method that was introduced to
        streamline the process of computing and caching <i>u/v</i>
        hash values for a given EST.  This method uses the variables
        \c u, \c v, and \c wordShift (all of them user configurable)
        along with a ESTCodec::NormalEncoder object to compute hash
        values into a std::vector. The vector is added to the \c
        uvCache hash map for future reference.

        \param[in] estIdx The zero-based index of the EST whose hash
        values is to be computed and cached.
    */
    void computeHash(const int estIdx);
    
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

    /** Instance variable to maintain the \i v parameter for the
	<i>u/v</i> heuristic.

    */
    static int v;

    static int wordShift;

    static int BitMask;

    /** Instance variable to store the number of bits to be shifted to
        create hash values.

        <p>This instance variable is set to the value of 2 * (\i v -
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
    /** The set of arguments specific to the UV heuristic.

        This instance variable contains a static list of arguments
        that are specific only to this analyzer class.  This argument
        list is statically defined and shared by all instances of this
        class.

        \note Use of static arguments and parameters renders this UV
        sample heuristic class not to be MT-safe.
    */
    static arg_parser::arg_record argsList[];    

    static int u;

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

    /** The hint key that is used to add hint for normal or
	reverse-complement D2 computation.

	This hint key is used to set a hint in the \c hints hash
	map. This string is defined as a constant to save compute time
	in the core \c runHeuristics method.
    */
    const std::string hintKey;
};
 
#endif
