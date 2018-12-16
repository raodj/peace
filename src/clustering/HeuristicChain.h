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

#include "Component.h"
#include <vector>
#include "Heuristic.h"
#include "HashMap.h"
#include "ParameterSetManager.h"

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
class HeuristicChain : public Component {
    friend class ClusteringSubSystem;
public:
    /** Enumeration of various hints provided by heuristics to other
        parts of PEACE.

        <p>This enumeration is used by heuristics to convey hints from
        the heuristics to other parts of PEACE (such as the heavy
        weight analyzers like D2 and TwoPassD2).  The hints are used
        by the analyzers to fine tune their operations thereby
        improving performance without compromising quality.</p>

        <p>Earlier (in PEACE 0.9) hints were conveyed using a
        hash_map.  However, with billions of calls the hash_map was
        adding considerable CPU cyles to the run.  Consequently, the
	x        hints were moved to a simple array (size based on \c LAST_HINT
        enumeration) to keep them as light weight as possible.</p>
    */
    enum HintKey {UNKNOWN_HINT, D2_DO_RC, MST_RC,

                  // Place new hints before this line.
                  LAST_HINT
    };

    /** The constructor.
		
	The heuristic chain that contain a set of heuristics and
	applies various operations on each heuristic in the chain.
    */
    HeuristicChain();
    	
    /** Create the heuristics in the chain.

        This method must be used to establish the chain of heuristics.
        This method is typically invoked from the \c main method.
        This method parses the names of the heuristic specified in the
        parameter \c heuristicStr and instantiates suitable heuristics
        via the HeuristicFactory.

        \param[in] heuristicStr A string containing the list of
        heuristics to be created in this chain. This string is
        typically specified by the user as a command line parameter.
        If this parameter is an empty string (\c ""), then this method
        performs no specific action.

        \return This method returns \c true if the chain was setup
        successfully.  On errors, it returns \c false.
    */
    bool setupChain(const std::string& heuristicStr);

    /** Initializes all the heuristics in the chain.

        This method is invoked when the SubSystem (that logically owns
        this component) is being initialized.  This method iterates
        over all the heuristics that have been added ot this chain and
        calls initialize() on each one of them.  If any one of the
        heuristics are unable to initialize correctly, then this
        method immediately returns an non-zero error code.

        \note Prior to invoking this method the set of cDNA fragments
        being processed must be set via a call to setESTList() method
        in this class.
        
        \return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    virtual bool initialize();

    /** Set the reference EST against which heuristics are to be run.

        This method is invokes the corresponding method on all the
        heuristics that have been added to this chain.  Setting the
        reference EST permits heuristics to build initial tables and
        other data structures for analyzing another EST (called via
        the shouldAnalyze) against the reference EST to determine if
        the heavy weight EST analyzer (d2 or CLU) must be run.

        \param[in] estIdx A pointer to the est that with which a large
        number of ESTs are going to be compared via the various
        heuristics that have been added to this chain.
    */
    virtual int setReferenceEST(const EST* est);

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

    /** Remove a given heuristic from this chain.

        This method permits a heuristic to be removed from this chain.
        The removed heuristic is deleted once it has been removed from
        the chain.

        \param[in] name The exact name of the heuristic to be removed
        from this chain.

        \return This method returns \c true if the heuristic was found
        and was successfully removed from the list.  Otherwise it
        returns \c false.
    */
    virtual bool removeHeuristic(const std::string& name);
    
    /** Determine whether the analyzer should perform core
        (computationally intensive) analysis, according to this
        heuristic chain.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.
        
        \param[in] otherEST The pointer to the EST with which the
        reference EST is to be compared.
        
        \return This method returns \c true if all of the heuristics
        say the EST pair should be analyzed, and \c false if any
        heuristic does not.  (Conceivably a subclass could extend this
        class and use some sort of "consensus" among heuristics.)
    */
    inline bool shouldAnalyze(const EST* otherEST) {
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
    inline void setHint(const HintKey hintKey, const int hintValue) {
        hints[int(hintKey)] = hintValue;
    }

    /** Obtain the value associated with a given hint.

        This method must be used to obtain the hint value associated
        with a given hint key.  The hint value is stored in the \c
        hintValue parameter.

        \param[in] hintKey The string identifier associated with the
        hint whose value is to be retrieved.
        
        \param[out] hintValue The parameter into which the hint value
        must be stored.  This method stores the last known hint value.
    */
    inline void getHint(const HintKey hintKey, int& hintValue) const {
        hintValue = hints[int(hintKey)];
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

    /** Obtain a reference to the parameter set manager.

        Obtain the parameter set manager that is used both by the
        Heuristic classes as well as the ESTAnalyzer classes.  Not all
        heuristics and analyzers use the parameter set manager.  They
        selectively use it depending on their needs.
        
        \return A pointer to the parameter set manager.  This pointer
        is always valid and is never NULL.
    */
    inline ParameterSetManager* getParamSetMgr() { return &psetMgr; }

    /** Set the ESTList to be used by this analyzer.

        This method must be used to set the list of cDNA fragments
        being processed by a given run of PEACE.  This list is
        typically a shared (between all components associated with a
        single instance of PEACE class) list of cDNA fragments that is
        used and processed by this class.

        \note This class maintains a pointer to this list internally
        until the finalize() method is called.  The methods in this
        class (and its child classes) do not add/remove entries in
        this list. They only update various attributes associated with
        each cDNA entry in the list.
        
        \param[in] estList The ESTList to be used by this analyzer.
        This pointer is used until the finalize() method is invoked.
    */
    void setESTList(ESTList* estList);

    /** Obtian a pointer to the list of cDNA fragments associated with
        all the heuristics.

        This method must be used to obtain a reference to the set of
        cDNA fragments to be processed by the various heuristics.

        \return The ESTList set for the various heuristics to
        operate. This pointer can be NULL if a valid ESTList has not
        yet been set.
    */
    inline ESTList* getESTList() const { return estList; }


    /** The run method that essentially calls the run method in all
	the heuristics in this chain.

	This method is a convinence method that has been introduced
	here to facilitate the process of dispatching run calls to all
	the filters in this chain.  This method is invoked from the
	ClusteringSubsystem::run() method which is called after
	initializations have successfully completed.

        \note This method assumes that the runtime context has been
        setup for the ClusteringSubSystem (as per the normal runtime
        operations).  The context is used by heuristics to obtain
        cross references to sub-components as needed.
	
	\return This method returns 0 (zero) on success. On errors it
	returns a non-zero error code.
    */
    virtual int run();
    
    /** The destructor.
        
        The destructor frees up all the heuristics added to this
        heuristic chain.
    */
    virtual ~HeuristicChain();
    
protected:
    /** An array to to store hints from heuristics.

        <p>This array is used to rapidly store and retrieve hints that
        are generated by various heuristics as they operate.  These
        hints may be used by other heuristics in the chain or by heavy
        weight EST analyzers to further optimize their operations.
        Entries are added via the \c setHint method.  Entries may be
        accessed via the \c getHint method.</p>
        
        <p>Earlier (in PEACE 0.9) hints were conveyed using a
        hash_map.  However, with billions of calls the hash_map was
        adding considerable CPU cyles to the run.  Consequently, the
        hints were moved to a simple array (size based on \c LAST_HINT
        enumeration) to keep them as light weight as possible.</p>
        
        \note Currently hints are restricted to be integers. Possibly
        at a later date this could be changed to store non-integer
        values as well.
    */
    int hints[LAST_HINT + 1];
    
    /** A shortcut reference to the shared list of cDNA fragments
        being analyzed.

        This list holds a pointer to the shared list of cDNA fragments
        currently being analyzed.  This pointer is initialized to NULL
        in the constructor.  A valid pointer is filled in when the
        initialize() method is invoked.  The pointer to the ESTList is
        set via a call to setESTList() method in this class.
        Typically, this pointer is set by the
        ClusteringSubSystem::initializeSubComponents() method.
    */
    ESTList* estList;
    
private:
    /** The vector containing a list of heuristics in the chain.

        This vector contains the list of hueristics assocaited with
        this chain.  Heuristics are added to the list via the
        addHeuristic() method.  The heuristics are used by the
        shouldAnalyze() method.
    */
    std::vector<Heuristic*> chain;

    /** The parameter set manager object.

        This instance variable contain the parameter set manager that
        is used both by the Heuristic classes as well as the
        ESTAnalyzer classes.  Not all heuristics and analyzers use the
        parameter set manager.  They selectively use it depending on
        their needs.
    */
    ParameterSetManager psetMgr;
};

#endif
