#ifndef PRIMES_HEURISTIC_H
#define PRIMES_HEURISTIC_H

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
// Authors:   Dhananjai M. Rao          raodm@miamiOH.edu
//
//---------------------------------------------------------------------

#include <string>
#include <vector>

#include "Heuristic.h"
#include "Utilities.h"
#include "PrimesHelper.h"
#include "HashMap.h"

/** \typedef HashMap<int, PrimesHelpers::ESTMetric> NearbyMap;

    \brief A hash map to provide rapid access to the set of
    nearby ESTs along with distance metrics.

    This typedef provides a convenient short cut to refer to a hash
    map that is used to quickly look-up nearby ESTs. The zero-based
    index of the EST is used as the key into this hash map.  Each
    entry in the hash map contains an object with pertinent
    information computed based on primes.
*/
typedef HashMap<int, PrimesHelper::ESTMetric> NearbyMap;

/** Heuristic that uses prime numbers to build k-features to detect
    similar reads.
*/
class PrimesHeuristic : public Heuristic, PrimesHelper {
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
        certain operations -- specifically, this heuristic builds the
        features for this reference strain for quick comparisons.

        \note This method must be called only after the initialize()
        method is called.
        
        \param[in] refEST The reference EST against which a large
        number of other EST entries are going to be compared via
        subsequent calls to the shoudAnalyze(const EST*) method.
        
        \return If the processing of the reference EST was sucessful,
        then this method returns 0.  Otherwise this method returns an
        error code.
    */
    virtual int setReferenceEST(const EST* refEST);

    /** The destructor.

        The destructor frees memory allocated for holding any dynamic
        data in the base class.
    */
    virtual ~PrimesHeuristic();

    /** Method to display statistics regarding operation of this
        heuristic.

        This method overrides the default implementation in the base
        class to report additional statistics recorded by this class.

        \note Derived heuristic classes may override this method to
        display additional statistics. However, the additional
        information must be displayed after the base class method has
        completed its task.
        
        \param[out] os The output stream to which the statistics
        regarding the heuristic is to be dumped.
    */
    virtual void printStats(std::ostream& os) const override;

protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead only the
        HeuristicFactory is permitted to instantiate this class as per
        the framework conventions.
        
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
    PrimesHeuristic(HeuristicChain *chain,
                    const std::string& name = "primes");
    
    /** Determine whether the analyzer should analyze, according to
        this heuristic.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.
        
        \param[in] otherEST An immutable pointer to the EST with which
        the reference EST is to be compared.
        
        \return This method returns true if the heuristic says the
        EST pair should be analyzed, and false otherwise.
    */
    virtual bool runHeuristic(const EST* otherEST);
    
private:
    /** The number of dimensions (aka features) in the primes-based
        heuristic to be generated.

        The instance variable tracks the number of features to be
        extracted for comparing two sequences.  This value is set via
        \c --pri-heur-features command-line argument.  The default value is
        arbitrarily set to 4.
    */
    int numFeatures;

    /** The prime number to be used for 'A' and 'T' nucleotides

        This isntance variable serves to hold the prime number to be
        used for generating the score for 'A' and 'T' base pair.  This
        value is set via \c --pri-heur-at command-line argument.  The
        default value is set to 71.
    */
    int atPrime;

    /** The prime number to be used for 'C' and 'G' nucleotides

        This instance variable serves to hold the prime number to be
        used for generating the score for 'C' and 'G' base pair.  This
        value is set via \c --pri-heur-cg command-line argument.  The
        default value is set to 113.
    */    
    int cgPrime;

    /** The distance threshold to determine if 2 sequences are
        sufficiently similar to warrant further analysis.

        This instance variable serves to hold a threshold value to be
        used to determine if two sequences are sufficiently similar.
        This value is set via \c --pri-heur-thresh command-line
        argument.  The default value is arbitrarily set to 1000,
        making most pairs similar.
    */
    float distThresh;

    /** Limit heuristic to use only the top top-N reads for further
        processing.  If this parameter is set then, this heuristic
        restricts itself to work only the top-N matches.  This is
        useful to restrict the number of entries to be checked,
        thereby improving overall performance of the system.  However,
        setting this to too small a value will degrade clustering
        quality.  Use the \c --pri-heur-topN command-line argument to
        set this parameter.
    */
    int topN;

    /** Limit heuristic to use only the top top-N-percentage reads for
        further processing.  This parameter is analogous to topN,
        except that 'N' is specified as a percentage of the reads
        below the given distThresh.  Use the \c --pri-heur-topN-per
        command-line argument to set this parameter.
    */    
    float topNper;

    /** An optional word length based on which position-weighted
        primes-based features are computed and used by this heuristic.
        If this value is -1, then weighed features are not used.  Use
        the \c --pri-heur-wordLen command-line argument to set this
        parameter.
    */
    int wordLen;
    
    /** The prime-numbers based features computed for the reference
        EST setup for this heuristic.  This value is computed in the
        setReferenceEST method in this class.
    */
    FloatVec refFeatures;

    /** The cached list of nearby ESTs along with metrics to speed-up
        heuristic processing.  This list is computed in the
        setReferenceEST and then used in the runHeuristic method.
    */
    NearbyMap nearest;

    /** The number of times the setReferenceEST method was called. */
    int numSetRef;

    /** The total time taken to complete calls to the setReferenceEST */
    double totSetRefTime;
};

#endif
