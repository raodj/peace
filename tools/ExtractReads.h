#ifndef EXTRACT_READS_H
#define EXTRACT_READS_H

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
#include "ArgParser.h"

/** \file ExtractReads.h

    \brief This file contains the class definition for a tool to
    extract a given list of reads from a given fasta file.
*/

/** A tool to extract reads in a data set into another data file.

    The overall operation of this tool is rather straightforward,
    and as follows:

    <ol>

    <li>Processes command-line arguments used by this tool.</li>
    
    <li> It loads all the ESTs from a given FASTA file.</li>

    <li>Extract the specified reads.</li>

    <li>Write the shuffled read to a given output file.</li>

    </ol>
*/
class ExtractReads : public Tool {
public:
    /** The main method to perform tasks for this tool.
        
        <p>This method is invoked from the ::main() method associated
        with PEACE-tools once it has detected that the tool to be used
        is this ExtractReads tool.  Any pending (unprocessed) command
        line parameters are passed to this method for its use.</p>

        This method uses an arg_parser object to process the pending
        command line arguments. If all the necessary arugments were
        supplied then it load reads, extracts specified reads, and
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

        \param[out] idxList The index positions of the ESTs to be
        extracted.

        \param[out] rndNum An optional random subset of reads to be
        extracted.

        \param[out] substrList An optional list of sub-strings to be
        used to extract entries.

        \return This method returns zero on success.
    */
    int parseArgs(int argc, char *argv[], std::string& inputPath,
                  std::string& outputPath,
                  ArgParser::StringList& idxList, int& rndNum,
                  ArgParser::StringList& substrList);

    /**
       Extracts a given list of random reads from the given input list.

       \param[in] rndNums The number of random reads to be extracted.

       \param[in] inputList The list of input EST reads.

       \param[out] output The output to where the data is to be written.
     */
    void extractRandomReads(const int rndNums, ESTList& inputList,
                            std::ostream& output) const;

    /**
       Extracts a given list of random reads from the given input list.

       \param[in] idxList The list of indexes in the inputList whose
       entries are to be extracted and written to output.

       \param[in] inputList The list of input EST reads.

       \param[out] output The output to where the data is to be written.
    */
    void extractIndex(const ArgParser::StringList& idxList,
                      ESTList& inputList, std::ostream& output) const;

    /**
       Extracts reads that contains any of the given set of substrings
       in their information.

       \param[in] substrList The list of sub-strings that are searched
       in each EST's info.  If any of the sub-strings are found in the
       info, then those ESTs are written to the output.

       \param[in] inputList The list of input EST reads.

       \param[out] output The output to where the data is to be written.
    */
    void extractSubstr(const ArgParser::StringList& substrList,
                       ESTList& inputList, std::ostream& output) const;

    /** Helper method to write out EST entry while honoring
        sub-fragment settings (if any).  This method uses the range in
        subSeqStart and subSeqEnd instance variables (set via
        command-ine arguments --sub-seq-beg and --sub-seq-end) to
        extract and write sub-fragements.

        \param[in] est The read to be written.

        \param[out] out The output stream to where the read is to be
        written.
    */
    void writeStrain(EST& est, std::ostream& out) const;

private:
    /** The default constructor for this class.

        This class is not meant to be directly instantiated.  Instead
        all the work is done by the static main() method.
        Consequently, the constructor has been made private to force
        users to use the main() method instead.  The constructor
        currently does not have a specific task to perform.
    */
    ExtractReads() {}

    /** The destructor.

        The destructor currently does not have a specific task to
        perform and is merely present to adhere to coding conventions.
    */
    ~ExtractReads() {}


    /** The starting index position of subsequence (if any). By
        default this value is -1 and is ignored.
    */
    int subSeqStart = -1;

    /** The ending index position of subsequence (if any). By default
        this value is -1 and is ignored.
    */
    int subSeqEnd = -1;
};

#endif
