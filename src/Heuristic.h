#ifndef HEURISTIC_H
#define HEURISTIC_H

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

/** The base class of all heuristics.

    This class must be the base class of all heuristics in the
    system. This class provides some default functionality that can be
    readily used by the heuristics. This class enables the
    HeuristicChain to manage a list of heuristics and dispatch method
    calls to various heuristics.
*/
class Heuristic {
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
    virtual void showArguments(std::ostream& os) = 0;

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
    virtual bool parseArguments(int& argc, char **argv) = 0;
    
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
    virtual int initialize() = 0;
    
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
    virtual int setReferenceEST(const int estIdx) = 0;
    
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
    bool shouldAnalyze(const int otherEST);

    /** Obtain human-readable name for this heuristic.

        This method must be used to obtain the human readable name set
        for this heuristic. This method essentially returns the value
        set when this heuristic class was instantiated.

        \return A human readable name associated with this heuristic.
    */
    const std::string& getName() const { return heuristicName; }
    
    /** The destructor.
        
        The destructor frees memory allocated for holding any dynamic
        data in the base class.
    */
    virtual ~Heuristic();

    /** Method to display statistics regarding operation of this heuristic.

        This method can be used to obtain a dump of the statistics
        gathered regarding the operation of this heuristic.  The
        typical statistic generated by heuristics includes:

        <ul>

        <li>The number of times the heuristic was called.  More
        specifically this value indicates the number of times the \c
        shouldAnalyze() method was invoked.</li>

        <li>The number of successful matches reported by this
        heuristic.  This number indirectly indicates the number of
        times other heuristics or the actual heavy weight algorithm
        was invoked.</li>
        
        </ul>

        \note Derived heuristic classes may override this method to
        display additional statistics. However, the additional
        information must be displayed after the base class method has
        completed its task.
        
        \param[out] os The output stream to which the statistics
        regarding the heuristic is to be dumped.
    */
    virtual void printStats(std::ostream& os) const;

    /** Method to obtain the count of times this heuristic was run.
    */
    inline int getRunCount() {
      return runCount;
    }

    /** Method to obtain the count of times this heuristic was
	successful (i.e. result indicated that heavy weight analysis
	should be undertaken).
    */
    inline int getSuccessCount() {
      return successCount;
    }
    
protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived Heuristic classes must be instantiated via the
        HeuristicFactory API methods.

        \param[in] name The human readable name for this heuristic.
        This name is used when generating errors, warnings, and other
        output messages for this heuristic.
    */
    Heuristic(const std::string& heuristicName);

    /** Determine whether the analyzer should analyze, according to
	this heuristic.
        
        This method is invoked from the shouldAnalyze() method to
        perform the actual heuristic analysis.  The analysis is
        performed between the refESTidx (set via earlier call to
        setReferenceEST) and otherEST.

        \note Derived heuristic classes must override this method and
        provide a proper implementation.
        
        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.
        
        \return This method must return \c true if the heuristic says
	the EST pair should be analyzed, and \c false otherwise.
    */    
    virtual bool runHeuristic(const int otherEST) = 0;

    /** The index of the reference EST in a given file. 

        This member object is used to hold the index of a reference
        EST in a given file.  The index values begin from 0 (zero).
        This member is initialized in the constructor and is changed
        by the setReferenceEST() id.
    */
    int refESTidx;
    
    /** The name of this heuristic.
        
        This instance variable contains the human recognizable name
        for this heuristic.  This value is set when the heuristic is
        instantiated (in the constructor) and is never changed during
        the life time of this heuristic.  This information is used when
        generating errors, warnings, and other output messages.
    */
    const std::string heuristicName;
    
private:
    /** Variable to track the number of times this heuristic was run.

        This instance variable is used to track the number of times
        this heuristic was run. This variable is initialized to zero
        in the constructor.  It is incremented each time the
        shouldAnalyze() method is invoked to run the the heuristic.
    */
    int runCount;

    /** Variable to track the number of times this heuristic passed.

        This instance vairable tracks the number of times the
        heuristic passed.  This number indirectly indicates the number
        of times other heuristics or the actual heavy weight algorithm
        was invoked.  This value is incremented in the \c
        shouldAnalyze() method each time the runHeuristic method
        returns a \c true value.
    */
    int successCount;
};

#endif
