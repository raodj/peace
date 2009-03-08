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
// Authors: James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include <vector>

// Forward declaration to make compiler happy
class Heuristic;

class HeuristicChain {

 public:
  
  /** Get a pointer to the instance of the heuristic chain.
      Since this class is a singleton, the constructor is private
      and the only way to obtain an instance of the class is through
      this method.
  */
  static HeuristicChain* getHeuristicChain();

  /** Add the given heuristic to the heuristic chain.

      \param[in] h The instance of class Heuristic that should be
      added to the heuristic chain.
  
      \return This method returns 0 if everything went well.
  */
  virtual int addHeuristicToChain(Heuristic* h);

  /** Determine whether the analyzer should analyze, according to
      this heuristic chain.

      This method can be used to compare a given EST with the
      reference EST (set via the call to the setReferenceEST())
      method.

      \param[in] otherEST The index (zero based) of the EST with
      which the reference EST is to be compared.

      \return This method returns true if all of the heuristics say the
      EST pair should be analyzed, and false if any heuristic does not.
      (Conceivably a subclass could extend this class and use some sort
      of "consensus" among heuristics.)
  */
  virtual bool shouldAnalyze(const int otherEST);

  /** The destructor.
   */
  virtual ~HeuristicChain();


 protected:

 private:
  /** The constructor.
      This is made private because the heuristic chain is a singleton,
      and should only be instantiated from the getHeuristicChain()
      static method.
  */
  HeuristicChain();

  /** The vector containing a list of heuristics in the chain.
  */
  static std::vector<Heuristic*> chain;

  /** The pointer to the singleton instance of this class.  Again, this
      is made private so that only methods of this class can access it.
  */
  static HeuristicChain* ptrInstance;

};

#endif
