#ifndef FILE_LOADER_H
#define FILE_LOADER_H

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
// as a result of using, modifying or distributing this software or
// its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include <string>

/** \file FileLoader.h
    
    \brief A utility class to read sequences from various file formats
    and populate the EST list.

    This file contains the class declaration for the FileLoader. The
    FileLoader is used to read cDNA sequence data from various file
    formats, such as: SFF or FASTA.  This class forms the glue between
    various file formats and the list of cDNA fragments maintained in
    the EST class.  Additional documentation regarding the API
    provided by this class is available with the class documentation.
*/

/** A class to read cDNA sequences from various file formats into
    memory.

    <p>This class is a helper class that is used to read cDNA sequence
    data from various file formats, such as: SFF and FASTA.  For each
    cDNA sequence read, this class creates a suitable EST object for
    it and populates the static list mainained in the EST class.  Note
    that this class does not directly read the data from various
    formats but delegates the tasks to different file readers.
    However, it forms the "glue" between the various file readers and
    the EST class.  The API provided by this class is used by both the
    clustering hierarchy of classes and the assembler hierarchy of
    classes.</p>

    <p>Note that the API methods in this class are all static methods.
    Consequently, the methods can be directly invoked without needing
    to create an instance of this class.  Moreover, one could argue
    that it may make sense to fold these methods into various reader
    classes.  However, it seems cleaner to separate the process of
    creating in-memory representations from the process of loading
    data from various file formats.</p>
*/
class FileLoader {
public:
    /** Method to load EST information from a FASTA file.

        This method can be used to load information regarding ESTs
        from a FASTA file.  The file name from where the data is to be
        loaded must be passed in as the parameter.

        \param[in] toolName The name of the tool that is attempting to
        load the file.  This value is used only to generate a
        meaningful error message.
        
        \param[in] fileName The file name of the FASTA file from where
        the EST information is to be uploaded.  If the input file name
        is literally the string \c "<none>" then, this method returns
        with \c true immeidately, without any further processing.

        \param[in] noMaskBases A boolean flag to indicate if
        lower-case characters must be masked out of reads.  Typically
        lower-case characters ('a', 't', 'c', 'g') are used to
        indicate bases that must be masked out of reads. If this flag
        is \c true (default) then these lower-case characters are
        converted to 'N' causing them to ignored by PEACE. If this
        flag is \c false, then these bases are converted to upper-case
        equivalents causing them to be processed normally. This flag
        is passed to EST::create which actually does the conversions.

        \param[in] randomizeNbases Boolean flag to indicate if all \c
        'N's must be randomly converted to a valid base.  Typically \c
        'N' characters in reads are retained and are masked
        over. However, in certain analyzers, these \c 'N' characters
        are randomly converted to a suitable base character.  This
        flag controls the behavior whether \c 'N' characters should be
        randomly converted.  If this value is \c false (the default)
        the characters are \b not converted.  However, if the flag is
        set to \c true then the \c 'N' characters are randomly
        converted to regular bases.
        
        \return This method returns true if all the ESTs were
        successfully loaded from the given file.
    */
    static bool loadFASTAFile(const std::string& toolName, const char* fileName,
                              const bool maskBases = true,
                              const bool randomizeNbases = false);

    /** Method to load fragments from a Standard Flowgram Format (SFF)
        file.

        This method can be used to load information regarding
        ESTs/sequences from a SFF file.  This method uses an instance
        of the SFFReader class to load the fragments.  The file name
        from where the data is to be loaded must be passed in as the
        parameter.

        \param[in] toolName The name of the tool that is attempting to
        load the file.  This value is used only to generate a
        meaningful error message.
        
        \param[in] fileName The file name of the SFF file from where
        the ESTs/454 sequence information is to be loaded.  If the
        input file name is literally the string \c "<none>" then, this
        method returns with \c true immeidately, without any further
        processing.

        \param[in] noMaskBases A boolean flag to indicate if
        lower-case characters must be masked out of reads.  Typically
        lower-case characters ('a', 't', 'c', 'g') are used to
        indicate bases that must be masked out of reads. If this flag
        is \c true (default) then these lower-case characters are
        converted to 'N' causing them to ignored by PEACE. If this
        flag is \c false, then these bases are converted to upper-case
        equivalents causing them to be processed normally. This flag
        is passed to EST::create which actually does the conversions.

        \param[in] randomizeNbases Boolean flag to indicate if all \c
        'N's must be randomly converted to a valid base.  Typically \c
        'N' characters in reads are retained and are masked
        over. However, in certain analyzers, these \c 'N' characters
        are randomly converted to a suitable base character.  This
        flag controls the behavior whether \c 'N' characters should be
        randomly converted.  If this value is \c false (the default)
        the characters are \b not converted.  However, if the flag is
        set to \c true then the \c 'N' characters are randomly
        converted to regular bases.

        \return This method returns true if all the ESTs were
        successfully loaded from the given file.
    */
    static bool loadSFFFile(const std::string& toolName, const char *fileName,
                            const bool maskBases = true,
                            const bool randomizeNbases = false);
    
protected:
    /** The file name used by API to indicate that there is no input
        file.

        This static member object contains a sentinel string \c
        "<none>" that is used by library API to indicate the absence
        of an actual input file to be loaded.
    */
    static const std::string IgnoreFileName;
    
private:
    /** The only constructor.

        The default constructor has is private in order to ensure that
        this class is never instantiated.  Instead, the static methods
        in this class must be directly used to load data from various
        file formats.
    */
    FileLoader() {}
    
    /** The destructor.
        
        The destructor is private (symmetric with the constructor) to
        ensure that objects of this class are never deleted.
    */
    ~FileLoader() {}
};

#endif
