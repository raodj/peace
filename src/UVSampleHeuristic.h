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
    UVSampleHeuristic(const int refESTidx, const std::string& outputFileName);
    
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

 private:

    static arg_parser::arg_record argsList[];
    
    /** A simple array to map characters A, T, C, and G to 0, 1, 2,
        and 3 respectively.

        This is a simple array of 255 entries that are used to convert
        the base pair encoding characters A, T, C, and G to 0, 1, 2,
        and 3 respectively to compute the hash as defined by CLU.
        This array is initialized in the constructor and is never
        changed during the life time of this class.
    */
    static char CharToInt[];

    static int u;

    static int v;

    static int wordShift;

};
 
#endif