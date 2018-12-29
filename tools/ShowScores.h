#ifndef SHOW_SCORES_H
#define SHOW_SCORES_H

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
#include "PrimesHelper.h"

/** \file ShowScores.h

    \brief This file contains the class definition for the ShowScores
    tool that print scores with respect to a given read.

    This file contains the class definition for the ShowScores tool.
    This tool prints the scores with respect to a given read.  Note
    that different analyzers and heuristics supported by PEACE can
    be used to generate the scores.
*/

/** A tool to print scores with respect to a given read.

    The overall operation of this tool is rather straightforward,
    and as follows:

    <ol>

    <li>Processes command-line arguments used by this tool. Most of
    the command-line arguments are passed onto initialize the PEACE
    object in the base class.</li>
    
    <li>It loads all the ESTs from a given FASTA file.</li>

    <li>It uses clustering subs class in PEACE, it generates and prints
    the n-dimensional features for each read in the FASTA file.
    </li>

    </ol>

    The generated ESTs (with alignment annotation) can be graphically
    viewed using the \c ShowAlignment tool.
*/
class ShowScores : public Tool, PrimesHelper {
public:
    /** The main method to perform tasks for this tool.
        
        <p>This method is invoked from the ::main() method associated
        with PEACE-tools once it has detected that the tool to be used
        is this ShowScores tool.  Any pending (unprocessed) command
        line parameters are passed to this method for its use.</p>

        This method processes the "--ref-read" parameter first. It
        then passes on rest of the arguments to initialize the PEACE
        object in the base class, which does all the setup.  Then it
        uses the analyzer set to compare a given reference read
        against all the other reads. It prints the resulting metrics.
        
        \param[in] argc The number of remaining command line arguments
        for use by this class/method.
        
        \param[in,out] argv The actual command line arguments for use
        by this class/method.
    */
    static int main(int argc, char *argv[]);
    
protected:
    /** The default constructor for this class.

        This class is meant to be instantiated only from the public
        static main() method.  Consequently, the constructor has been
        made private to force users to use the main() method instead.
        The constructor currently does not have a specific task to
        perform and is merely present to adhere to coding conventions.

        \param[in] argc The number of remaining command line arguments
        to be passed to PEACE for initialization.
        
        \param[in,out] argv The actual command line arguments for use
        by PEACE for initialization.
    */
    ShowScores(int argc, char *argv[]);

    /** Print the scores with respect to given reference read. */
    void showScores(const int refID);
    
    /** The destructor.

        The destructor frees any dynamic memory allocated to store the
        data members in this class.  The destructor currently does not
        have a specific task to perform and is merely present to
        adhere to coding conventions.
    */
    ~ShowScores() {}

private:
};

#endif
