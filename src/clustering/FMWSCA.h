#ifndef FMWSCA_H
#define FMWSCA_H

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
    using ESTAnalyzer::getMetric;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    ~FMWSCA();


    /** Add valid command line arguments for this analyzer.

        This method must be used to add all valid command line options
        that are supported by this analyzer.  Note that derived
        classes may override this method to add additional command
        line options that are applicable to it.  This method is
        invoked when the clustering sub-system is initialized.
        
        \note Derived EST analyzer classes may override this method to
        display help for their custom command line arguments.  When
        this method is overridden don't forget to call the
        corresponding base class implementation to add common options.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);
    
protected:
    /** Flag to perform case sensitive comparisons.

        This flag indicates if the FMWSCA must perform case sensitive
        string comparisons.  By default this analyzer performs case
        insensitive comparisons.  The default can be changed by the
        user by sepcifying a suitable command line argument.  This
        member is initialized to false.  However, its value may be
        changed via the \c --case command-line parameter.
    */
    bool caseSensitive;

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
    virtual float getMetric(const std::string& refFrame,
                            const std::string& otherFrame,
                            const int wordSize);
    
private:
    /** The default constructor.

        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ESTAnalyzerFactory is a
        friend of this class, an object can be instantiated via teh
        ESTAnalyzerFactory::create() method.
    */
    FMWSCA();
};


#endif
