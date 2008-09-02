#ifndef FMWSCA_H
#define FMWSCA_H

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
//
//---------------------------------------------------------------------------

#include "FWAnalyzer.h"

/** FMWSCA: Framed, Multi-Word String Compare Analyzer.

    <p>This analyzer provides a conventional EST analyzer that compares
    EST base pair data using a given number of base pairs.  The total
    number of base pairs to be compared is called a Frame.  A frame is
    broken into a sequence of words.  The frame size and word size (in
    terms of number of base pairs) is specified as command line
    arguments.  This analyzer compares the tail (3' end) of a given
    EST sequence with the beginning (5' end) of all other ESTs in a
    given file.  The file must be in FASTA format.</p>

    <p>This class has been implemented by extending the FWAnalyzer
    base class.  The base class hierarchy provides most of the
    standard functionality involved in reading FASTA files and
    generating formatted output.</p>
    
    \note An instance of this class is typically created via the
    ESTAnalyzerFactory class.
*/
class FMWSCA : public FWAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    ~FMWSCA();

    /** Display valid command line arguments for this analyzer.

        This method must be used to display all valid command line
        options that are supported by this analyzer (and its base
        classes).
        
        \note This method calls the base class's showArguments first.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.

        This method is used to process command line arguments specific
        to this EST analyzer.  This method is typically used from the
        main method just after the EST analyzer has been instantiated.
        This method consumes all valid command line arguments.  If the
        command line arguments were valid and successfully processed,
        then this method returns \c true.

        \note This method consumes its custom command line arguments
        first and then call's the base class's parseArguments() method.
        
        \param[inout] argc The number of command line arguments to be
        processed.

        \param[inout] argc The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.
    */
    virtual bool parseArguments(int& argc, char **argv);
    
protected:
    /** Flag to perform case sensitive comparisons.

	This flag indicates if the FMWSCA must perform case sensitive
	string comparisons.  By default this analyzer performs case
	insensitive comparisons.  The default can be changed by the
	user by sepcifying a suitable command line argument.  This
	member is initialized to false.  However, its value may be
	changed by the parseArguments() method.
    */
    static bool caseSensitive;

    /** Method to compare two frames and compute similarity.

        For each word in the refFrame, this method searches for number
        of occurrences of the word in otherFrame.  It then divides the
        count by frameSize^2.  If two frames are identical it results
        in 100.  If thw two frames are completely different the result
        would be a 0 (zero).

        \note This method overrides the default implementation in the
        FWAnalyzer base class.
        
        \param[in] refFrame The reference frame for comparison
        purposes.  Note that the reference frame is always a constant
        in a given set of caparisons.  Consequently, certain analyzers
        can pre-compute and reuse metrics to make analysis fast.

        \param[in] otherFrame The other frame for comparison.  This
        frame is always guaranteed to be from a different EST than the
        refFrame.

        \param[in] wordSize The size of a word within the given frame.
        This value is always greater than 0 (zero) and less than frame
        size.
        
        \return This method is return a similarity metric in the range
        0 to 100 between the given frame and the refFrame.
    */
    using FWAnalyzer::analyze;
    virtual float analyze(const std::string& refFrame,
			  const std::string& otherFrame,
			  const int wordSize);
    
private:
    /** The default constructor.

        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ESTAnalyzerFactory is a
        friend of this class, an object can be instantiated via teh
        ESTAnalyzerFactory::create() method.

	
	\param[in] refESTidx The reference EST index value to be used
	when performing EST analysis.  This parameter should be >= 0.
	This value is simply passed onto the base class.

	\param[in] outputFile The name of the output file to which the
	EST analysis data is to be written.  This parameter is ignored
	if this analyzer is used for clustering.  If this parameter is
	the empty string then output is written to standard output.
	This value is simply passed onto the base class.
    */
    FMWSCA(const int refESTidx, const std::string& outputFile);

    /** The set of c arguments for this EST analyzer.

        This instance variable contains a static list of arguments
        that are specific to the FMWSCA EST analyzers.
    */
    static arg_parser::arg_record argsList[];    
};


#endif
