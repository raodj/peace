#ifndef SCORE_CALCULATOR_H
#define SCORE_CALCULATOR_H

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

#include "DecagonException.h"
#include "Contig.h"
#include "MegaBlastRunner.h"
#include "OnDemandESTList.h"
#include <fstream>

/** The top-level class that controls and coordinates all the
    clustering and assembly activities.

    This is the top-level class constituting PEACE.  The primary
    method for triggering the various processing performed by this
    system is the run() method. The run() method (along with other
    helper methods in this class) perform the following core tasks:

    <ol>

	<li>The sub-systems constituting PEACE are instantiated when this
	class is created.</li>

	<li>It gathers the various top-level command line arguments that
	are accepted by various sub-systems by invoking the
	SubSystem::addCommandLineArguments() method.</li>
	
    <li>The command-line arguments supplied to PEACE (and passed-in
    via the global main() method) are processed for the initial,
    top-level options.  These options provide necessary information
    about various components of sub-systems to be used, such as:
    assembler to use, cluster generator to use, the set of filters to
    be applied, and the cDNA data files to be processed. Additional
    command-line arguments to the components of sub-systems are not
    processed and are deferred for processing by the components when
    the SubSystem::run() method is invoked.</li>

    <li>Assuming the top-level arguments parsed in the previous step
    are valid, the top-level components constituting the filtering,
    assembly, and clustering sub-systems are instantiated when the
    SubSystem::initialize() method is invoked on each sub-system. This
    method also adds command-line arguments that are specific to each
    component.  The sub-systems are initialized in the following
    order: <ol>

    <li>First, the output sub-system is initialized to create any
    debugging streams or logs (in one of many supported output
    formats) as indicated by the user.</li>

    <li>Next the file management sub-system is initialized to load any
    data files that the user may have specified.</li>
    
    <li>Next, the filtering sub-system is initialized.</li>
    
    <li>The cluster maker sub-system is initialized.</li>
    
    <li>Finally, the assembler sub-system is initialized.</li>
    
    </ol>
    
    </li>

    <li>Next, command-line arguments (based on parameters collected
    from each one of the sub-systems in the previous step during
    initialization) are processed and validated. The command-line
    arguments are processed in a 2-phased approach because the actual
    parameters vary depending on the type of sub-system(s) being
    utilized for a specific run of PEACE.</li>

    <li>If any errors occur in any of the aforementioned steps, then
    all valid command line arguments are displayed and no further
    processing is performed.</li>
    
    <li>Once all the sub-systems have been successfully initialized,
    they are all permitted to run in the same order as enumerated in
    the previous step.</li>
    
    <li>Finally, all the sub-systems are finalized in the following
    order:

	<ol>
        <li>The input sub-system is finalized. </li>

        <li>Next, the filtering sub-system is finalized.</li>

        <li>The cluster maker sub-system is finalized.</li>

        <li>The assembler sub-system is finalized.</li>

        <li>Lastly, the output sub-system is finalized closing all
        outputs.</li>
        		
        </ol>

	</li>
    
    </ol>
*/
class ScoreCalculator {
public:
    /** The constructor.

        The constructor does not have very many tasks to perform (as
        all the member objects appropirately initialize themselves).
        It establishes the various command-line arguments supported by
        this tool.
    */
  ScoreCalculator();
	
    /** The primary method that performs the various tasks for this
        tool.

        This method is typically invoked from the main()
        function. This acts as the top-level coordinator for the
        various activities performed by the score calculator.

        \param[in,out] argc The number of command-line arguments that
        were passed-in to the process.  This instance variable
        determines the number of entries in the argv parameter.  This
        value is changed to reflect the number of arguments
        successfully parsed and consumed by this class.

        \param[in,out] argv The array of command-line arguments to be
        processed. Processed arguments are removed from this list.

        \return This method returns \c true if the score calculations
        were successfully completed.  Otherwise (if routine errors are
        detected) then this method returns \c false.

        \exception DecagonException This method throws an exception
        when serious errors are encountered.
    */
    bool run(int& argc, char *argv[]) throw(DecagonException);
    
    /** The destructor.
        
        The destructor does not have very many operations to perform
        and is primarily present to adhere to coding conventions.
    */
    virtual ~ScoreCalculator();

protected:
    /** Helper method to actually parse command line arguments.

        This is a helper method that is invoked from the run() method
        actually parse the command-line arguments.  This method
        creates a temporary ArgParser object to actually parse the
        command-line arguments.  The arguments are ultimately
        populated in the various instance variables in this class.

        \return This method returns \c true if the command-line
        arguments were parsed successfully and all the necessary
        arguments were supplied. On errors or missing arguments, this
        method returns \c false.
    */
    bool parseCmdLineArgs(int& argc, char *argv[]);

	/** Helper method to load transcripts from a given FASTA file into
		a given list.

		This a helper method that is used to load reads from a given
		FASTA file. For example, it is invoked from from the run()
		method to load the source transcripts from the given source
		transcript FASTA file.

        \param[in] fileName The FASTA file from where the reads are to
        be loaded.

        \param[out] targetList The target EST list into which the
        reads are to be loaded by this method.

        \throws DecagonException This method throws an exception if
        the file is not a FASTA file or if errors occur when loading
        the reads.
    */
	void loadReads(const std::string& fileName, ESTList& targetList)
        throw(DecagonException);

	/** Helper method to load just the core contig information from a
		given ACE file into a given contig list.

		This a helper method that is used to load contigs from an
		assembler-generated contig file. For example, it is invoked
		from from the run() method to load the generated contigs from
		the given contig file.  This method checks to see if the
		contig file is in FASTA file format. If so, it loads it as a
		FASTA file. Otherwise it assumes the file is an ACE file and
		loads it as an ACE file.

        \param[in] contigFileName The contig file (ACE or FASTA file
        format) from where the contigs are to be loaded.

        \param[out] contigList The target contig list into which the
        contigs are to be loaded by this method.

        \note In order to reduce memory footprint, this method does
        not retain alignment information and individual reads. Only
        the aggregate contig information (such as: ID, consensus
        sequence, nucleotide distribution) are retained.

        \throws DecagonException This method throws an exception if
        the contigs could not be successfully read from the specified file.
    */
	void loadContigs(const std::string& contigFileName, ContigList& contigList)
        throw(DecagonException);

    /** Helper method to write the contigs to a FASTA file.

        <p>This is a helper method that is invoked from the main run()
        method. It is used to write the contigs (loaded from an
        assembler-generated output file) to a FASTA file. Writing the
        contigs to a FASTA file is needed in order to run Megablast,
        as Megablast requires generated-contigs to be in a FASTA file
        format (to compare with source transcripts).</p>

        <p>This method writes each contig to the file. The contigIDs
        are used as the FASTA header. The consensus sequence for each
        contig is written to the FASTA file. Quality values for the
        consensus sequence is not written.</p>
        
        \param[in] contigList The list of contigs from which the
        contigs are to be written in FASTA format.

        \param[in] fastaFileName The output FASTA file name to which
        the contig information is to be written.
    */
    void writeContigsToFASTAFile(const ContigList& contigList,
                                 const std::string& fastaFileName)
        throw (DecagonException);
    
    /** Find a contig whose ID begins with a given prefix.

        Megablast generates the first word (delimited by a space) in
        the header of a given contig/read as an identifier. This
        method is a helper method that searches the list of contigs in
        a given list to find a matching contig.

        \param[in] prefix The prefix string to search for in the
        ID of each contig.

        \param[in] list The list of contigs to search in.

        \return This method returns a reference to the contig whose ID
        starts with the given prefix.

        \exception DecagonException This method throws an exception if
        a contig with the given prefix was not found in the list.
    */
    const Contig& findContig(const std::string& prefix,
                             const ContigList& list) const
        throw(DecagonException);

    /** Find a read whose ID begins with a given prefix.

        Megablast generates the first word (delimited by a space) in
        the header of a given contig/read as an identifier. This
        method is a helper method that searches the list of reads in a
        given list to find a matching entry.

        \param[in] prefix The prefix string to search for in the
        info/headers of the reads.

        \param[in] list The list of reads to search in.

        \return This method returns a reference to the read whose ID
        starts with the given prefix.

        \exception DecagonException This method throws an exception if
        a read with the given prefix was not found in the list.
    */
    const EST& findRead(const std::string& prefix,
                        const ESTList& list) const
        throw(DecagonException);
    
    /** Generate scores for each pair in the given set.

        This method a refactored method that was introduced to
        streamline the code in the run() method.  This method iterates
        over the list of entries in given set (generated by processing
        the output from Mega Blast) and generates score for each pair
        of entries.

        \param[in] goodPairs The set of Megablast Data entries (that
        list highlight matching transcript for each contig) to be used
        to generate scores.

        \param[out] os The output stream to which the socres are to be
        written.        
    */
    void generateScores(const MBDSet& goodPairs, std::ostream& os);
    
private:
    /** Command-line argument value indicating the contig (ACE or
        FASTA) file generated by a genomic-assembler.

        This instance variable is populated (in the parseCmdLineArgs()
        method) to refer to the value supplied with the
        "--gen-contig-file" command-line parameter. This value
        indicates the path to an ACE or FASTA file containing contigs
        generated by a genomic assembler.  This program automatically
        detects the file format (ACE or FASTA) and loads the contigs.
    */
    std::string genContigFile;
    
    /** An optional output file stream to which the results from
        analysis are to be written.

        This stream is set ot refer to an output file when if the user
        specifies an alternative output file via the
        <code>--output-file</code> command-line parameter. The stream
        is setup in the parseCmdLineArgs() method.  This stream is
        used in the run() method.
    */
    std::ofstream outputFile;
    
    /** Command-line argument value indicating the source FASTA file
        containing transcripts for validation.

        This instance variable is populated (in the parseCmdLineArgs()
        method) to refer to the value supplied with the
        "--src-trans-file" command-line parameter. This value
        indicates the path to the FASTA file containing source
        transcripts against which the generated contigs are to be
        compared for scoring.  Typically, the source transcript file
        is the file from where cDNA fragments were generated.  This
        file is used to create a BLAST DB (as needed) via the
        MakeBlastDB class.
    */
    std::string srcTransFile;

	/** The list of fragments in the source transcript file.

		This list is populated with the fragments from the source
		transcript file.  The reads are loaded from the srcTransFile
		via the loadReads() method in this class.
	*/
	OnDemandESTList srcTransList;

	/** The list of contigs in the assembler-generated output file.

		This list is populated with the contigs read an
		assembler-generated output file. The contigs are typically
		read from an ACE file.  The contigs are loaded from the given
		genContigFile via the loadContigs() method in this class.
	*/
	ContigList genContigList;
};

#endif
