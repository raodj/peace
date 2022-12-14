#ifndef MEGA_BLAST_RUNNER_H
#define MEGA_BLAST_RUNNER_H

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

/** \file MegaBlastRunner.h

	\brief A simple class to run Megablast and process its output.

	This file contains the declaration for the MegaBlastRunner class.
	It is used to run Megablast to identify best transcript-contig
	pairs for further analysis and scoring.
*/

#include "MegaBlastData.h"
#include "BlastCommon.h"

#include <cstdio>
#include <string>
#include <set>

/** \typedef std::set<MegaBlastData, LessMegaBlastData> MBDSet

    \brief A convenience typedef for a set of MegaBlastData entries.

    This set set contains the best matching entry for each generated
    contig reported by megablast.  This set is returned by call to the
    MegaBlastRunner::runMegaBlast() method.
*/
typedef std::set<MegaBlastData, LessMegaBlastData> MBDSet;


/** Helper class to run Megablast and process its output.

	This is a helper class that is used to run Megablast and process
	its outputs.  This class is used only when cDNA fragments have
	been generated from a known transcript or the transcript is
	known. Megablast is used to determine the best matching transcript
	to a given contig (generated by an assembler). This matching
	process is needed because various genomic-assemblers generate
	contigs in different orders and with different IDs.  Consequently,
	this phase of using Megablast to identify best transcript-contig
	pairs is essential.
*/
class MegaBlastRunner : BlastCommon {
public:
    /** The default constructor.

        The constructor does not perform any special operation other
        than to initialize members to default initial value.
    */
    MegaBlastRunner();

    /** The destructor.

        The destructor closes the connection to the Megablast process
        in case it was not already closed by the runMegaBlast() method.
    */
    ~MegaBlastRunner();

    /** Primary method run Megablast and get is processed outputs.

        This method runs Megablast as a child process (after obtaining
        command-line to run via BlastCommon::getMegablastCmdLine())
        and directly reads the outputs generated by it. It processes
        the outputs to create MegaBlastData objects. Next, for each
        contig (identified by unique IDs), the best matching
        transcript entry is maintained and returned in the supplied
        vector.

        \param[in] transBlastDB The full path to the source transcript
        BLAST database that contains the source transcripts to which
        the contigs (generated by an assembler) is to be matched.

        \param[in] genContigFile The full path to the FASTA file that
        contains the contigs generated by the assembler. Each contig
        is processed by Megablast to identify a suitable matching
        source transcript.

        \param[out] results This set is populated with entries created
        from results obtained from running Megablast. The set will
        contain the processed results from Megablast.  This set will
        contain only a single, best-matching entry for each
        contig. Existing entries (if any) in the set are not cleared
        by this method.

        \exception DecagonException This method throws a decagon
        exception if errors occur when running Megablast and
        processing its outputs.
    */
    void
    runMegaBlast(const std::string& transBlastDB,
                 const std::string& genContigFile,
                 MBDSet& results)
        throw(DecagonException);
                 
protected:
    /** Helper method to parse Megablast output.
        
        This method is repeatedly called from the runMegaBlast()
        method to process each line of output generated by
        Megablast. This method parses the line of output (assuming it
        is generated using a \c "-m 9" command-line argument to
        Megablast and creates a MegaBlastData object.

        \param[in] line The line of Megablast output to be parsed by
        this method.

        \return This method returns a MegaBlastData object that has
        been populated with necessary information from the line of
        Megablast output passed-in as the parameter.

        \exception DecagonException This method throws a decagon
        exception if the line of data could not be successfully
        parsed.
    */
    MegaBlastData parseLine(const std::string& line) const
        throw(DecagonException);

    /** Helper method to maintain the best entry for a given contig.

        This method is a helper method that is repeatedly invoked (for
        each line of Megablast output parsed) to retain the best
        matching source transcript entry for a given
        assembler-generated contig.  This method performs the
        following operations:

        <ol>

        <li>If the given set of entries (namely bestSoFar parameter)
        does not already have an entry for the given contig
        (entry.getGenContigID()) then this method adds the entry to
        the set and returns.</li>

        <li>If an entry already exists for the generated contig, then
        this method replaces the existing entry if the e-value of the
        new entry is better.  If e-values are the same then the entry
        with the better bit-score is retained.</li>

        </ol>

        \param[in,out] bestSoFar The set of per-generated-contig
        MegaBlastData entries that are the best matches found thusfar.

        \param[in] entry The new entry that has just been parsed to be
        used for updating the values in bestSoFar parameter.
    */
    void setBestEntry(MBDSet& bestSoFar, const MegaBlastData& entry) const;
    
private:
    /** The pipe to the Megablast process.

        This pointer maintains the OS-level stream to the outputs
        generated by Megablast.  The stream is created via a suitable
        call to popen() function.  This value is initialized to NULL
        in the constructor.  It is set to a valid value by the
        runMegaBlast() mehtod.  If it is not closed by the
        runMegaBlast() method (due to exceptions) then the destructor
        closes this pipe.
    */
    FILE* pipe;
};

#endif
