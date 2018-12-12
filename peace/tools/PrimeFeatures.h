#ifndef PRIME_FEATURES_H
#define PRIME_FEATURES_H

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

/** \file PrimeFeatures.h

    \brief This file contains the class definition for the
    PrimeFeatures tool that print n-dimensional features based on
    prime numbers.

    This file contains the class definition for the PrimeFeatures
    tool.  This tool prints the n-dimensional features based on prime
    numbers for each EST in a given data file.
*/

/** A tool to print n-dimensonal features based on prime numbers.

    The overall operation of this tool is rather straightforward,
    and as follows:

    <ol>

    <li>Processes command-line arguments used by this tool.</li>
    
    <li> It loads all the ESTs from a given FASTA file.</li>

    <li>Using the PrimesHelper class in PEACE, it generates and prints
    the n-dimensional features for each read in the FASTA file.
    </li>

    </ol>

    The generated ESTs (with alignment annotation) can be graphically
    viewed using the \c ShowAlignment tool.
*/
class PrimeFeatures : public Tool, PrimesHelper {
public:
    /** The main method to perform tasks for this tool.
        
        <p>This method is invoked from the ::main() method associated
        with PEACE-tools once it has detected that the tool to be used
        is this PrimeFeatures tool.  Any pending (unprocessed) command
        line parameters are passed to this method for its use.</p>

        This method uses an arg_parser object to process the pending
        command line arguments. If all the necessary arugments were
        supplied then it instantiates an PrimeFeatures object and uses
        various methods in the class (and its base class) to perform
        the operations of this tool.
        
        \param[in] argc The number of remaining command line arguments
        for use by this class/method.
        
        \param[in,out] argv The actual command line arguments for use
        by this class/method.
    */
    static int main(int argc, char *argv[]);

    /** Method to generate the n-dimensional features for all the
        reads.
        
        This method is invoked to generate the n-dimensional features
        for all the reads.

        \note This method assumes that the EST data has already been
        successfully loaded from a given FASTA file.

        \param[in] numFeatures The number of dimensions or features to
        be generated.

        \param[out] out The output stream to which the annotated ESTs
    */
    void generateFeatures(int numFeatures, std::ostream &os = std::cout);
    
protected:
    /** The default constructor for this class.

        This class is meant to be instantiated only from the public
        static main() method.  Consequently, the constructor has been
        made private to force users to use the main() method instead.
        The constructor currently does not have a specific task to
        perform and is merely present to adhere to coding conventions.

        \param[in] atPrime The prime number to be used for 'A' and 'T'
        nucleotides.

        \param[in] cgPrime The prime number to be used for 'C' and 'G'
        nucleotides.        
    */
    PrimeFeatures(const int atPrime, const int cgPrime);

    /** The destructor.

        The destructor frees any dynamic memory allocated to store the
        data members in this class.  The destructor currently does not
        have a specific task to perform and is merely present to
        adhere to coding conventions.
    */
    ~PrimeFeatures() {}

private:
};

#endif
