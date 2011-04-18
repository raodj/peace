#ifndef FASTA_FILE_H
#define FASTA_FILE_H

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

#include <fstream>
#include "InputFile.h"

/** The class for reading EST entries from a given FASTA text file.

    This class serves as a concrete implementation of an InputFile to
    read EST entries from a given FASTA text file.  This class
    implements the various API methods stipulated by the InputFile
    base class.

    \note This class is not really meant to be directly instantiated.
    Instead, the InputFileFactory class is designed to serve as the
    primary entry point for instantiating a suitable, derived
    InputFile object.
*/
class FASTAFile : public InputFile {
public:
    /** API method to detect if a given file is a FASTA file.

        This method can be used to determine if the a given file
        appears to be a valid FASTA file. This method performs only a
        basic check (and does not read the full file) to ensure that
        the given file appears to the a FASTA file.

        \note The check for a FASTA file is very simple as there is no
        special identifiers within the data to indicate if a file is a
        FASTA file. Consequently, for auto-detecting file types, it is
        best to check for a FASTA file format the last in the set of
        available file format options.

        \param[in] fileName The path (can be relative or absolute) to
        the file to be verified if it is a valid FASTA file.
        
        \return If a single nucleotide sequence is successfully read
        from the file then this method returns \c true.  method
        returns \c true.  If the given file is not readable or if the
        it is not in a FASTA file format then this method returns \c
        false.
    */
    static bool isFASTA(const std::string& fileName);
    
    /** Constructor to create an InputFile for reading data from a
        FASTA text file.

        \param[in] fileName The full path to the FASTA file from where
        the data is to be read.
    */
    FASTAFile(const std::string& fileName);

    /** Determine if the input file is still good or if it has
        encountered any errors.
        
        This method can be used to determine if this input file is
        still in "good" state or if an error has occured preventing it
        from operating correctly.
        
        \return This method returns true if the file is good. If an
        error had occured this method returns false.
    */
    bool good() const { return !fastaFile.fail(); }
    
    /** Determine if another entry is pending to be read from this
        InputFile.

        This method can be used to determine if another <i>new</i>
        entry is pending to be read from this input file.  This method
        returns true, if the input file is not at EOF.

        \return This method returns \c true if another entry is
        pending to be read.  If all entries have been read then this
        method returns \c false.
    */
    bool hasNextEntry();

    /** Destructor.

        The destructor for this class.  The destructor merely closes
        the file handle that was opened in the constructor.
    */
    ~FASTAFile();
    
protected:
    /** Obtain the current position from where the data for the next
        entry is going to be read.

        \return The current file pointer offset location from where
        the next entry is going to be read.
        
        \see setCurrPos()
    */
    size_t getCurrPos();

    /** Set the current position from where the next entry is to be
        read.

        \param[in] position The position that was reported by the
        getCurrPos() implementation (and was saved by the InputFile)
        in an earlier call.  This is essentially the file offset
        within the specified input file.

        \return This method returns \c true if the position for next
        entry to be read has been reset by this method. On errors or
        if the position could not be reset then this method returns
        false.

        \see getCurrPos()
    */
    bool setCurrPos(const size_t position);
    
    /** Helper method to read information about a given entry.

        This method is appropriately implemented by this derived class
        to read the information for the current entry from the FASTA
        file.

        \param[out] info The header information about the EST to be
        read from the file.

        \param[out] sequence The nucleotide sequence information that
        has been read from the file.

        \param[out] quality The nucleotide quality information for
        each of the nucleotide bases in the sequence.  Currently, this
        method simply clears out this vector.

        \return This method returns true if the data was successfully
        read.  On errors it returns false.
    */
    bool readEntry(std::string& info, std::string& sequence,
                   QualityVector& quality);
    
private:
    /** The OS-level file pointer used to read data from the FASTA
        file.
        
        This instance variable that used to hold the file pointer from
        where the FASTA data is being read.  This instance variable is
        initialized in the constructor once the FASTA file has been
        successfully opened for reading. The file is closed in the
        destructor.
    */
    std::ifstream fastaFile;
};

#endif
