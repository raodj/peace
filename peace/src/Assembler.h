#ifndef ASSEMBLER_H
#define ASSEMBLER_H

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

#include "arg_parser.h"

// Some forward declarations to keep compiler fast and happy
class BatonAnalyzer;

/** \file Assembler.h
    
    \brief The common base class and primary Application Program
    Interface (API) for all gene assemblers.

    This file contains the class declaration for the Assembler class.
    This class is the primary API for all assemblers in PEACE.
    Additional documentation regarding the API provided by this class
    is available with the class and method documentation.
*/

/** The base class of all assemblers in PEACE.

    <p>This class must be the base class of all assemblers in the
    system. This class provides some default functionality that can be
    readily used by each assembler.  However, this class is very
    generic, thereby providing the necessary flexibility for
    implementing different types of assemblers.</p>

	<p>This class cannot be directly instantiated.  Instead a suitable
    derived assembler can be created via the AssemberFactory.  Refer
    to the documentation on the AssemblerFactory for additional
    details.</p>

    @see AssemblerFactory
    @see BatonAssembler
*/
class Assembler {
public:
    /** Display valid command line arguments for this assembler.

        This method must be used to display all valid command line
        options that are supported by this assembler.  Note that
        derived classes may override this method to display additional
        command line options that are applicable to it.  This method
        is typically used in the system core when displaying usage
        information.

        \note Derived assembler classes <b>must</b> override this
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
        to this assembler.  This method is typically used from the
        core system just after an derived assembler has been
        instantiated.  This method must consumes all valid command
        line arguments applicable to the implementation.  If the
        command line arguments were valid and successfully processed,
        then this method returns \c true.

        \note Derived assembler classes <b>must</b> override this
        method to process any command line arguments that are custom
        to their operation.  When this method is overridden don't
        forget to call the corresponding base class implementation to
        display common options.
        
        \param[in,out] argc The number of command line arguments to be
        processed.

        \param[in,out] argv The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.  This method returns true if all arguments
        are consumed successfully.
    */
    virtual bool parseArguments(int& argc, char **argv);
    
    /** Primary method to perform cDNA assembly.

        This method is the primary interface method that is invoked
        (on each and every MPI process) to perform the complete
        assembly and output generation process.  This method is
        invoked once all the command line parameters are processed by
        the various subsystems constituting PEACE.  This method is a
        pure-virtual method.  Therefore all derived assembler classes
        must override this method to perform all the necessary
        operations.

		\note This method must be invoked only after the initialize()
		method has been invoked.

        \return This method returns zero if the assembly process was
        successfully completed.  Otherwise this method returns a
        non-zero error code.
    */
    virtual int assemble() = 0;

	/** A method to handle initialization tasks for the Assembler.
		
		This method is called after an assembler object has been
		instantiated but \e before the cDNAs (to be assembled) have
		been loaded.  This method performs the following tasks:

        <ol>

        <li> First it loads all the ESTs and other cDNA fragments to
		be processed from the specified input file.</li>

        <li>Next it creates the output file (if specified, otherwise
        output is dumped to standard out) to write assembler
        data.</li>
        
        </ol>
        
		\return This method returns zero on success. On errors, this
		method returns a non-zero value.
    */
    virtual int initialize();
	
    /** The destructor.

        The destructor frees memory allocated for holding any data in
        this base class.
    */
    virtual ~Assembler();
    
protected:
    /** The default constructor.

        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived assembler classes must be instantiated via the
        AssemblerFactory API methods.

        \param[in] name The human readable name for this assembler.
        This name is used when generating errors, warnings, and other
        output messages for this object.

        \param[in] outputFileName The file name to which output must
        be written.  If a valid output file is not specified (that is
        outputFile is an empty string), then results are written to
        standard output.  The outputFileName is simply copied to the
        outputFileName member object and is later on used in the
        initialize method.
    */
    Assembler(const std::string& name, const std::string& outputFileName);

    /** The name of this assembler.

        This instance variable contains the human recognizable name
        for this assembler.  This value is set when this object is
        instantiated (in the constructor) and is never changed during
        the life time of this object.  This information can be used
        when generating errors, warnings, and other output messages.
    */
    const std::string name;
    
    /** The file name to which results are to be written.
        
        This member object is used to hold the file name to which the
        primary output from the assembly process is to be written.
        This member is initialized to NULL.  However, the value is
        changed by the parseArguments method depending on the actual
        value specified by the user via the \c --output command line
        parameter.  The file stream to which the data is to be written
        is created in the
    */
    const std::string outputFileName;

    /** The Standard Flowgram Format (SFF) input file from where cDNA
        fragments are to be loaded.

        This static member is used to hold the SFF file name from
        where all the cDNA fragments are to be loaded.  The actual
        file name is set by the parseArguments method if the the user
        has specified a file name via the \c --sffFile comamnd line
        argument.

        \note The user must specify an input file either via the \c
        --sffFile or the \c --fastaFile (but not both) command line
        options.
    */
    static char* sffFileName;

    /** The FASTA file from where cDNA data is to be loaded.

        This member object is used to hold the FASTA file name from
        where all the cDNA data is to be loaded.  The actual file name
        is set by the parseArguments method if the user has specified
        a file name via the \c --fastaFile comamnd line argument.

        \note The user must specify an input file either via the \c
        --sffFile or the \c --fastaFile (but not both) command line
        options.
    */
    static char* fastaFileName;
    
    /** Flag to indicate if lower-case characters must be masked out of
        reads.
        
        Typically lower-case characters ('a', 't', 'c', 'g') are used
        to indicate bases that must be masked out of reads. This
        notation is used by DUST (part of NCBI BLAST) utility that
        identifies and tags low complexity regions with lower-case
        letters. If this flag is \c false (default) then these
        lower-case characters are converted to 'N' causing them to
        ignored by PEACE. If this flag is \c true, then these bases
        are converted to upper-case equivalents. This flag is passed
        to EST::create which actually does the conversions.  The
        command line parameter to change the default value from \c
        false to \c true is \c --no-mask-bases.
    */
    static bool noMaskBases;

	/** Flag to indicate if all 'N's must be randomly converted to
        a valid base.
	
        Typically 'N' characters in reads are retained and are masked
        over by the analyzers. However, in certain analyzers, the
        assumption is that these 'N' characters are randomly converted
        to a suitable base character.  This flag controls the behavior
        whether 'N' characters should be randomly converted.  If this
        value is \c false (the default) the characters are \b not
        converted.  However, if the flag is set to \c true (either via
        a command line parameter or through the static method), then
        the 'N' characters are randomly converted to regular bases.
        The command line parameter to enable random conversion of \c
        'N's is \c --rand-N.
    */
    static bool randomizeNbases;
    
private:
    /** The set of common arguments for all assemblers.

        This instance variable contains a static list of arguments
        that are common all the assembler classes.  The common
        argument list is statically defined and shared by all
        assembler instances.

        \note The use of the static argument list may make the
        assembler class hierarchy not to be MT-safe.
    */
    static arg_parser::arg_record commonArgsList[];
    
    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \note This method is intentionally has a dummy implementation.
        
        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.
        
        \return Reference to this.
    */
    Assembler& operator=(const Assembler& src);
};

#endif
