#ifndef HEURISTIC_CHAIN_H
#define HEURISTIC_CHAIN_H

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

#include <vector>
#include "Heuristic.h"
#include "HashMap.h"

/** Class that manages a list of heuristics.

    <p>This class represents a list of heuristics that are used to try
    and minimize the number of times EST analysis is performed on a
    given pair of ESTs.  Note that EST analysis (using d2 or clu) is a
    time consuming task and it must be performed sparingly to reduce
    runtimes (particularly when large sets of ESTs need to be
    clustered).  The heuristics try to elimintate performing EST
    analysis on pairs that are most likely unrelated (and will not be
    in the same cluster).</p>

    <p>Each heuristic in the chain implements a specific heuristic and
    returns a \c true or \c false value. If the result is \c false,
    then detailed EST analysis is not performed. However, if the
    result is \c true, then the next heuristic in the chain is used to
    further ensure that the comparison is really needed. It is the
    responsiblity of the caller to ensure that the chain of hueristics
    is correctly set so that the cheaper (possibly less aggressive)
    heuristic is run first followed by more computationally-expensive
    (and more aggressively filtering) heuristics are run.</p>
*/
class HeuristicChain {
public:
    /** Create the heuristics in the chain.

        This method must be used to establish the chain of heuristics.
        This method is typically invoked from the \c main method.
        This method parses the names of the heuristic specified in the
        parameter \c heuristicStr and instantiates suitable heuristics
        via the HeuristicFactory.

        \param[in] heuristicStr A string containing the list of
        heuristics to be created in this chain. This string is
        typically specified by the user as a command line parameter.
        If this parameter is \c NULL, then this method performs no
        specific action.

        \param[in] refESTidx The index of the EST that is to be used
        to serve as the starting EST for clustering.  In a
        Minimum-Spanning-Tree (MST) clustering scheme, this EST will
        be the root of the MST.

        \param[in] outputFile The file to which the results of
        clustering are to be written.
    */
    static HeuristicChain*
    setupChain(const char* heuristicStr, const int refESTidx,
               const std::string& outputFile);

    /** Get a pointer to the instance of the heuristic chain.
        
        Since this class is a singleton, the constructor is private
        and the only way to obtain an instance of the class is through
        this method.  The heuristic chain is available only after the
        \c setupChain method (that is invoked from main right after
        command line arguments are validated) has successfully
        completed its operation.  Until such time this method simply
        returns \c NULL.

        \return The process-wide unique pointer to the heuristic chain.
    */
    static inline HeuristicChain* getHeuristicChain() {
        return ptrInstance;
    }
    
    /** Display valid command line arguments for heuristics in the
	chain.

        This method simply calls the showArguments method on each
	heuristic in the chain.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Permits heuristics in the chain to process command line
       arguments.

       This method iterates over the heuristics that have been added
       to this chain and invokes parseArguments() method on each one
       of them.  This permits each heuristic in the chain to receive
       and process command line parameters targeted for the
       heuristics.

       \param[in,out] argc The number of command line arguments
       currently present in argv (the parameter list).

       \param[in,out] argv The list of command line arguments to be
       consumed by various heuristics, if they find parameters
       intended for their use.

       \return This method returns \c true if all the heuristics in
       the chain successfully processed command line arguments.  If an
       incorrect command line argument is received by any one of the
       heuristics then this method returns \c false to flag an error.
    */
    virtual bool parseArguments(int& argc, char **argv);

    /** Initializes all the heuristics in the chain.

        This method iterates over all the heuristics that have been
        added ot this chain and calls initialize() on each one of
        them.  If any one of the heuristics are unable to initialize
        correctly, then this method immediately returns an non-zero
        error code.

        \return This method returns zero on success. On errors this
        method returns a non-zero value.
    */
    virtual int initialize();

    /** Set the reference EST against which heuristics are to be run.

        This method is invokes the corresponding method on all the
        heuristics that have been added to this chain.  Setting the
        reference EST permits heuristics to build initial tables and
        other data structures for analyzing another EST (called via
        the shouldAnalyze) against the reference EST to determine if
        the heavy weight EST analyzer (d2 or CLU) must be run.

        \param[in] estIdx The index of the reference EST.
    */
    virtual int setReferenceEST(const int estIdx);

    /** Add the given heuristic to the heuristic chain.

        This method permits the heuristic chain to takes ownership of
        a given heuristic object by added it to its internal chain.

        \note The heuristic chain takes ownership of the object
        therefore that the heuristic pointer passed to this method
        must not be deleted by the caller.
        
        \param[in] heuristic The instance of class Heuristic that
        should be added to the heuristic chain. 
        
        \return This method returns \c true if the heuristic was
        successfully added. On errors this method returns \c false.
    */
    virtual bool addHeuristic(Heuristic* heuristic);

    /** Determine whether the analyzer should perform core
        (computationally intensive) analysis, according to this
        heuristic chain.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.
        
        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.
        
        \return This method returns \c true if all of the heuristics
        say the EST pair should be analyzed, and \c false if any
        heuristic does not.  (Conceivably a subclass could extend this
        class and use some sort of "consensus" among heuristics.)
    */
    inline bool shouldAnalyze(const int otherEST) {
        for (size_t i = 0; (i < chain.size()); i++) {
            if (!chain[i]->shouldAnalyze(otherEST)) {
                // Immediately stop when a heuristic says no
                return false;
            }
        }
        // If all heuristics say yes, return true
        return true;
    }

    /** Method to display statistics regarding operation of all the
        heuristics in this chain

        This method can be used to obtain a dump of the statistics
        gathered regarding the operation of all the heuristics in this
        chain.  The typical statistic generated by heuristics
        includes:

        <ul>

        <li>The number of times the heuristic was called.  More
        specifically this value indicates the number of times the \c
        shouldAnalyze() method was invoked.</li>

        <li>The number of successful matches reported by this
        heuristic.  This number indirectly indicates the number of
        times other heuristics or the actual heavy weight algorithm
        was invoked.</li>
        
        </ul>

        \param[out] os The output stream to which the statistics
        regarding the heuristics is to be dumped.

	\param[in] rank The rank of the process for which the
	statistics is being displayed.
    */
    void printStats(std::ostream& os, const int rank) const;

    /** Set a hint to be used by other algorithms.

        This method can be used by various heuristics to set hints for
        use by other algorithms after their operation.  The hints are
        typically set by heuristics at the end of their \c
        runHeuristics method.  The hints are then used by other
        heuristics or by the heavy weight EST analysis algorithms
        (such as D2).

        \param[in] hintKey The string identifier associated with the
        hint. Using string identifiers makes the use of hints a bit
        more programmer friendly and makes the code more readable.

        \param[in] hintValue The actual value to be associated with
        the hintKey.
    */
    inline void setHint(const std::string& hintKey, const int hintValue) {
        hints[hintKey] = hintValue;
    }

    /** Obtain the value associated with a given hint.

        This method must be used to obtain the hint value associated
        with a given hint key.  If the hint value is available, then
        it is stored in the \c hintValue parameter.

        \param[in] hintKey The string identifier associated with the
        hint whose value is to be retrieved.

        \param[out] hintValue The parameter into which the hint value
        (if available) must be stored.

        \return This method returns \c true if the requested hint was
        found and \c hintValue was updated.  Otherwise this method
        returns \c false.
    */
    inline bool getHint(const std::string& hintKey, int& hintValue) const {
        StringIntMap::const_iterator entry = hints.find(hintKey);
        if (entry != hints.end()) {
            hintValue = entry->second;
            return true;
        }
        // Hint not found.
        return false;
    }

    /** Method to obtain pointer to a given heuristic object.

	This method can be used to obtain a pointer to a specific
	heuristic class present in this chain.  If the heruistic does
	not exist then this method returns NULL.

	\note The caller must \b not delete the returned pointer.

	\param[in] name The name associated with a given heuristic.
	
	\return If the heursitic was found then this method returns a
	valid (non-NULL) pointer to the heuristic object. If the
	heuristic was not found, then this method returns NULL.
    */
    Heuristic* getHeuristic(const std::string& name) const;
    
    /** The destructor.

        The destructor frees up all the heuristics added to this
        heuristic chain.
    */
    virtual ~HeuristicChain();
    
protected:
    /** A hash map to store hints from heuristics.

        This hash map is used to rapidly store and retrieve hints that
        are generated by various heuristics as they operate.  These
        hints may be used by other heuristics in the chain or by heavy
        weight EST analyzers to further optimize their operations.
        Entries are added via the \c setHint method.  Entries may be
        accessed via the \c getHint method.

        \note Currently hints are restricted to be integers. Possibly
        at a later date this could be changed to store non-integer
        values as well.
    */
    StringIntMap hints;
    
private:
    /** The constructor.
        This is made private because the heuristic chain is a singleton,
        and should only be instantiated from the getHeuristicChain()
        static method.
    */
    HeuristicChain();
    
    /** The vector containing a list of heuristics in the chain.

        This vector contains the list of hueristics assocaited with
        this chain.  Heuristics are added to the list via the
        addHeuristic() method.  The heuristics are used by the
        shouldAnalyze() method.
    */
    std::vector<Heuristic*> chain;
    
    /** The pointer to the singleton instance of this class.

        Again, this is made private so that only methods of this class
        can access it. The getHeuristicChain() method in this class
        must be used to obtain an instance of this class.
    */
    static HeuristicChain* ptrInstance;
};

#endif
