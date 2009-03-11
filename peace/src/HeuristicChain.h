#ifndef HEURISTIC_CHAIN_H
#define HEURISTIC_CHAIN_H

//---------------------------------------------------------------------------
//
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

#include <vector>
#include "Heuristic.h"

class HeuristicChain {
public:
  
    /** Get a pointer to the instance of the heuristic chain.
        
        Since this class is a singleton, the constructor is private
        and the only way to obtain an instance of the class is through
        this method.
    */
    static HeuristicChain* getHeuristicChain();

    /** Display valid command line arguments for this heuristic
	chain.  This method simply calls the showArguments method
	on each heuristic in the chain.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    virtual bool parseArguments(int& argc, char **argv);

    virtual int initialize();

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
    virtual bool shouldAnalyze(const int otherEST);
    
    /** The destructor.

        The destructor frees up all the heuristics added to this
        heuristic chain.
    */
    virtual ~HeuristicChain();
    
    
protected:
    // Currently this class has no protected members.
    
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
    static std::vector<Heuristic*> chain;
    
    /** The pointer to the singleton instance of this class.

        Again, this is made private so that only methods of this class
        can access it. The getHeuristicChain() method in this class
        must be used to obtain an instance of this class.
    */
    static HeuristicChain* ptrInstance;
};

#endif
