#ifndef INPUT_FILE_FACTORY_H
#define INPUT_FILE_FACTORY_H

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

#include <iostream>
#include "Component.h"
#include "ArgParser.h"
#include "Utilities.h"

// For declaration to keep compiler happy & fast
class InputFile;

/** Primary component of InputSubSystem for instantiating a suitable
	InputFile object to read cDNA data from a file.

	This class is a component of the input sub-system.  It provides a
	centralized and convenient mechanism to instantiate a suitable
	class (derived from InputFile base class) to read cDNA data from a
	given cDNA file.
 */
class InputFileFactory : public Component {
    friend class InputSubSystem;
public:
    /** The destructor.

        The destructor is merely present to adhere to coding
        conventions.  Currently, this does not have any special tasks
        to perform.
    */
    ~InputFileFactory();
    
    /** Method to instantiate a suitable InputFile

        This method must be used to instantiate a suitable InputFile
        object to read EST data from a given file.  This method uses
        the fileName (parameter) to automatically detect the file type
        and suitably instantiate an InputFile object to read EST data.
        If the file could not be read by any of the available readers
        then this method returns NULL.

        \param[in] fileName The name of the EST data file along with
        optional path (the path can be relative or absolute) from
        where the data is to be read.

        \return If the file type could be detected and the initial
        processing of the data was successful then this method returns
        a valid InputFile pointer. Otherwise it generates an error
        message and returns NULL.
    */
    InputFile* create(const std::string& fileName);

    /** Helper method to create an input file of a specific type.

        This method can be used to open a given file with a specific
        file type format.

        \param[in] fileName The name of the EST data file along with
        optional path (the path can be relative or absolute) from
        where the data is to be read.

        \param[in] format The file type format for the specified
        file. The currently valid formats are: \c fasta, \c sff.  If
        the input file is not of the specified format then this method
        generates an error.

        \return If the file format is compatible with the data file
        and and the initial processing of the data was successful then
        this method returns a valid InputFile pointer. Otherwise it
        generates an error message and returns NULL.
    */
    InputFile* create(const std::string& fileName, const std::string& format);
    
    /** Add the set of command line parameters for this component.

        This method is typically invoked by the InputSubSystem when it
        is being requested for command line arguments.  This method
        adds the various command line arguments that can be used to
        load cDNA entries from one or more data files.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    void addCommandLineArguments(ArgParser& argParser);

    /** Initialize the operations of this component.
        
        This method is invoked when the InputSubSystem is being
        initialized.  This method (along with helper methods) performs
        the necessary operations to load cDNA data from the input
        files (if any) specified by the user.
        
        \return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    bool initialize();


    /** Helper method to load list of files specified in a file list
        with a given optional file format.

        This method is a helper method that is used to load files from
        a given file list. This method is currently called from the
        initialize() method. This method loads one file at a time from
        the given file into the list of ESTs maintained in the runtime
        context (RuntimeContext::estList). If any one file fails to
        load successfully this method immediately returns with \c
        false without processing the remainder of the list.
        
        \note The cDNA fragments that are not logically owned by this
        process (if running on a multi-process MPI job) are
        unpouplated by default to minimize memory footprint.
	
        \param[in] fileNames The list of strings containing the names
        of the files to be loaded by this method.

        \param[in] startIndex The starting index position of the EST
	from where the full data is to be retained. For other entries,
	the header information and nucleotide sequences are loaded
	on-demand.  This is done to reduce peak memory footprint
	(thereby enabling the processing of large data files). This
	value must be less than endIndex.  This value is with
	reference to the full list of ESTs (not just relative to one
	file).
        
        \param[in] endIndex The ending index position of the EST upto
	(and <b>not</b> including) which the full data is to be
	retained. For other entries, the header information and
	nucleotide sequences are loaded on-demand.  This value must be
	greater than startIndex. This value is with-reference-to the
	full list of cDNA fragments.  By default the startIndex and
	endIndex are set such that all cDNA fragments are loaded
	on-demand.
        
        \return If all the files in the list are successfully loaded,
        then this method returns \c true.  If an error occurs when
        processing one of the files, this method immediately returns
        with \c false.
    */
    bool loadFiles(ArgParser::StringList fileNames,
                   const std::string& format = "",
                   const long startIndex = MAX_READS,
                   const long endIndex   = MAX_READS);
	
    /** Wind-up the operations of this component.
        
        This method is invoked when the InputSubSystem is being
        finalized.  Currently, this method does not have any special
        tasks to perform.
    */
    void finalize();
    
protected:
    /** The default constructor.

        The default constructor has is private in order to ensure that
        this class is instantiated only by its friends.  Currently
        only the InputSubSystem class can instantiate this component.
        Instead, the static methods in this class must be directly
        used to instantiate a suitable InputFile reader object.
    */
    InputFileFactory();
    
private:
    /** The list of FASTA files specified by the user at the command
        line.

        This instance variable is passed to the command line argument
        parser and is used to hold the list of FASTA file names (with
        optional absolute or relative path) as specified by the user.
        The list is set in the addCommandLineArguments() method and is
        filled-in when PEACE performs the initial round of command
        line argument processing.  The list is processed by the
        initialize() method.
    */
    ArgParser::StringList fastaFileNames;

    /** The list of SFF files specified by the user at the command
        line.

        This instance variable is passed to the command line argument
        parser and is used to hold the list of SFF file names (with
        optional absolute or relative path) as specified by the user.
        The list is set in the addCommandLineArguments() method and is
        filled-in when PEACE performs the initial round of command
        line argument processing.  The list is processed in the
        initialize() method.
    */
    ArgParser::StringList sffFileNames;

    /** The list of data files specified by the user at the command
        line.

        This instance variable is passed to the command line argument
        parser and is used to hold the list of generic data file names
        (with optional absolute or relative path) as specified by the
        user.  The list is set in the addCommandLineArguments() method
        and is filled-in when PEACE performs the initial round of
        command line argument processing.  The list is processed in
        the initialize() method.

        \note The primary difference between this list and others is
        that the format of the data files in list list are
        automatically detected.
    */
    ArgParser::StringList estFileNames;
    
    /** Flag to indicate if lower-case characters must be masked out of
        reads.
        
        Typically lower-case characters ('a', 't', 'c', 'g') are used
        to indicate bases that must be masked out of reads. This
        notation is used by DUST (part of NCBI BLAST) utility that
        identifies and tags low complexity regions with lower-case
        letters. If this flag is \c false (default) then these
        lower-case characters are converted to 'N' causing them to
        ignored by PEACE. If this flag is \c true, then these bases
        are converted to upper-case equivalents. This flag is used to
        set options in InputFile objects created by this class.
    */
    bool noMaskBases;

    /** Flag to indicate if all \c 'N's must be randomly converted to
        a valid base.
        
        Typically \c 'N' characters in reads are retained and are
        masked over by the analyzers. However, in certain analyzers,
        the assumption is that these \c 'N' characters are randomly
        converted to a suitable base character.  This flag controls
        the behavior whether \c 'N' characters should be randomly
        converted.  If this value is \c false (the default) the
        characters are \b not converted.  However, if the flag is set
        to \c true (either via a command line parameter or through the
        static method), then the \c 'N' characters are randomly
        converted to regular bases.
    */
    bool randomizeNbases;

    /** Flag to indicate if ESTs should be loaded and retained in
        memory instead of loading on-demand.
        
        By default, EST information is not retained in memory but
        loaded on-demand.  This is done to minimize peak memory
        footprint when processing large files, but at the cost of
        increase in I/O.  However, with RAM becoming abundant, it may
        make sense to hold gigabytes of reads in RAM.  This flag is
        used to override the default on-demand operation and just
        retain all reads in RAM to improve overall performance (at the
        expense of memory usage).  This flag is overridden via the
        command-line argument --no-ondemand
    */
    bool noOnDemand;    
};


#endif
