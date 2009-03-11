#ifndef HEURISTIC_FACTORY_H
#define HEURISTIC_FACTORY_H

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

#include <iostream>

// For declaration to make compiler happy.
class Heuristic;

class HeuristicFactory {
public:
    /** Method to instantiate a suitable heuristic.
        
        This method must be used to instantiate a suitable heuristic
	instance.  This method uses the name (parameter) to
        suitably instantiate a heuristic.  If the name is not
        valid, then this method returns NULL.
        
        \param[in] name The name of the heuristic to be instantiated.
        
        \param[in] refESTidx The reference EST index value to be used
        when instantiating a heuristic.  Index values start from 0
        (zero).  If this value is negative then this method returns
        NULL.

        \param[in] outputFileName The target file to which the
        analysis report is to be written (if any).  Note that this
        parameter may be ignored if this heuristic is used to generate
        clusters.  If the value of outputFileName is "" (empty string)
        then the outputs are streamed to standard out.
    */
    static Heuristic* create(const char* name, const int refESTidx,
                             const std::string& outputFileName);
    
    /** Method to display the list of heuristics available.

        This method is typically used to display the list of
        heuristics currently available.  This method is typically used
        in the main() method when displaying usage information.
        
        \param[out] os The output stream to which the list of heuristic
        names must be written.
    */
    static void displayList(std::ostream& os);
    
protected:
    // Currently this class has no protected members.
    
private:
    /** The default constructor.

        The default constructor has is private in order to ensure that
        this class is never instantiated.  Instead, the static methods
        in this class must be directly used to instantiate a heuristic.
    */
    HeuristicFactory() {}

    /** The destructor.

        The destructor is private to ensure that objects of this class
        are never deleted.
    */
    ~HeuristicFactory() {}
};

#endif
