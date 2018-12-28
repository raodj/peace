#ifndef RANDOMIZE_READS_H
#define RANDOMIZE_READS_H

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

#include "Tool.h"

/** \file RandomizeReads.h

    \brief This file contains the class definition for a tool to
    randomly shuffle reads in a given file for roubust testing.

    This file contains the class definition for the RandomizeReads
    tool.  This tool randomly shuffles the entries in a given input
    file.  Shuffling the data ranodmly ensures that testing with data
    sets is not invalidated due to any inherent ordering already
    present in the data set.
*/

/** A tool to randomize reads in a data set for robust testing.

    The overall operation of this tool is rather straightforward,
    and as follows:

    <ol>

    <li>Processes command-line arguments used by this tool.</li>
    
    <li> It loads all the ESTs from a given FASTA file.</li>

    <li>Randomly shuffles the reads.</li>

    <li>Writes the shuffled read to a given output file.</li>

    </ol>

    The generated ESTs (with alignment annotation) can be graphically
    viewed using the \c ShowAlignment tool.
*/
class RandomizeReads : public Tool {
public:
    /** The main method to perform tasks for this tool.
        
        <p>This method is invoked from the ::main() method associated
        with PEACE-tools once it has detected that the tool to be used
        is this RandomizeReads tool.  Any pending (unprocessed)
        command line parameters are passed to this method for its
        use.</p>

        This method uses an arg_parser object to process the pending
        command line arguments. If all the necessary arugments were
        supplied then it load reads, randomly shuffles them, and
        writes them to a given output file.
        
        \param[in] argc The number of remaining command line arguments
        for use by this class/method.
        
        \param[in,out] argv The actual command line arguments for use
        by this class/method.
    */
    static int main(int argc, char *argv[]);
    
protected:
    /** An internal refactored helper method to process command-line
        arguments.

        \param[in] argc The number of remaining command line arguments
        for use by this class/method.
        
        \param[in,out] argv The actual command line arguments for use
        by this class/method.

        \param[out] inputPath The input file name specified by the user.

        \param[out] outPath The output file name specified by the user.        
    */
    static int parseArgs(int argc, char *argv[], std::string& inputPath,
                         std::string& outputPath);
    
private:
    /** The default constructor for this class.

        This class is not meant to be directly instantiated.  Instead
        all the work is done by the static main() method.
        Consequently, the constructor has been made private to force
        users to use the main() method instead.  The constructor
        currently does not have a specific task to perform.
    */
    RandomizeReads() {}

    /** The destructor.

        The destructor currently does not have a specific task to
        perform and is merely present to adhere to coding conventions.
    */
    ~RandomizeReads() {}
};

#endif
