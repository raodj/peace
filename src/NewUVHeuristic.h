#ifndef NEW_UV_HEURISTIC_H
#define NEW_UV_HEURISTIC_H

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
class NewUVHeuristic : public Heuristic {
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
        
        \param[in,out] argc The number of command line arguments to be
        processed.
        
        \param[in,out] argv The array of command line arguments.
        
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
    virtual ~NewUVHeuristic();

protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived Heuristic classes must be instantiated via the
        HeuristicFactory API methods.

        \param[in] name The human readable name for this heuristic.
        This name is used when generating errors, warnings, and other
        output messages for this heuristic.

		\param[in] outputFileName The output file to which any
		analysis information is to be written. Currently this
		parameter is unused.
    */
    NewUVHeuristic(const std::string& name,
				   const std::string& UNREFERENCED_PARAMETER(outputFileName));
    
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

	/** Instance variable to track the lenght of the reference EST.

		This instance variale is a convenience variable that is used
		to track the length of the reference EST.  Since EST length is
		a commonly used operation, possibly this instsance variable
		can be moved into the EST class.
	*/
    int refESTLen;

	/** Instance variable to track the lenght of the "other" EST with
		which the reference EST is being analyzed.

		This instance variale is a convenience variable that is used
		to track the length of a given EST.  Since EST length is a
		commonly used operation, possibly this instsance variable can
		be moved into the EST class.
	*/	
    int otherESTLen;
	
    /** The threshold for number of common words, used by the
        <i>u/v</i> heuristic.
		
        This parameter contains the minimum number of common words
        that two cDNA fragments must contain in order to "pass" the
        <i>u/i</i> heuristic.
    */	
    int u;

	/** Number of base pairs to skip when comparing words.

		The number of base pairs that are to be skipped (in each pass)
		when building hashes (in the computeHash() method) and checking
		for common words.
	*/
    int wordShift;

	/** The number of passes that determines the set of words being
		compared.

		The passes are useful for longer fragments. Rather than
		checking each word, the passes permit every other word (or
		every third word if passes == 3) and using the resulting
		number of common words to determine if further analysis is
		required.
	*/
    int passes;
    
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
    */
    static int v;

    static int BitMask;

    static int argumentU;

    static int argumentWordShift;

    static int argumentPasses;

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
    /** The set of arguments specific to the UV heuristic.

        This instance variable contains a static list of arguments
        that are specific only to this analyzer class.  This argument
        list is statically defined and shared by all instances of this
        class.

        \note Use of static arguments and parameters renders this UV
        sample heuristic class not to be MT-safe.
    */
    static arg_parser::arg_record argsList[];

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
