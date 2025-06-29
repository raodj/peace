#ifndef INPUT_FILE_H
#define INPUT_FILE_H

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

#include <map>
#include <string>
#include <vector>

#include "EST.h"

/** The base class for various types of input files supported by
    PEACE.

    This class serves as the base class for the various input files
    supported by PEACE.  In addition to providing a generic interface
    to other components, the base class also provides the following
    features:

    <ul>

    <li>Common infrastructure to normalize nucleotide bases read from
    an input file.  See the normalizeBases method for additional
    details.</li>
    
    <li>Common infrastructure to randomize 'N' nucleotide entries in
    the sequence and randomly convert it to one of \c ATCG.</li>

    <li>Repopulate the header information and sequence information for
    a given EST.</li>
    
    </ul>

    \note This class is merely a base class and is not meant to be
    directly instantiated.  Instead, one of derived classes must be
    instantiated and used.  The InputFileFactory class is designed to
    serve as the primary entry point for instantiating a suitable,
    derived InputFile object.
*/
class InputFile {
public:    
    /** Obtain the path to the file from where the EST data is being
        read.

        This method merely returns the file name set when this base
        class was instantiated.
        
        \return The path to the file from where the data is being
        read.
     */
    const std::string& getFileName() const { return fileName; }

    /** Determine if the input file is still good or if it has
        encountered any errors.

        This method can be used to determine if this input file is
        still in "good" state or if an error has occured preventing it
        from operating correctly.

        \return This method returns true if the file is good. If an
        error had occured this method returns false.
    */
    virtual bool good() const = 0;
    
    /** Determine if another entry is pending to be read from this
        InputFile.

        This method can be used to determine if another <i>new</i>
        entry is pending to be read from this input file.

        \note This method is a pure virtual method in this interface
        class and is suitably implemented by one of the derived
        classes.
    */
    virtual bool hasNextEntry() = 0;
    
    /** Read the next entry from the input file.

        This method must be used to read the next entry from this
        InputFile and create a suitable EST entry for it. This method
        returns a valid EST entry only when the hasNextEntry() method
        returns \c true.

        \param[in] estID The unique integer ID value to be set for the
        new EST entry being created. It is best to ensure that the
        estIDs are created consecutively.

        \return On successfully reading the data this method
        instantiates a new EST entry (allocates memory using the \c
        new operator) and populates with the necessary data.

        \see repopulate()
    */
    EST* getNextEntry(const int estID);

    /** The EST entry to be repopulate by this method.

        This method can be used to repopulate the following
        information in this EST entry:

        <ul>

        <li>The header information about the EST (this information is
        returned by the EST::getInfo() method.</li>

        <li>The nucleotide sequence information about the EST.</li>

        <li>The nucleotide quality sequence (if available).</li>

        </ul>

        \param[in,out] est The EST entry to be repopulated by this
        method. If the repopulation is successful then this method
        returns true. Otherwise it returns false.
    */
    bool repopulate(EST& est);

    /** Set options to be used when processing EST entries.

        This method can be used to set up options for processing
        nucleotide bases in the reads. Typically, the options are set
        right after an input file is created. The options are applied
        the next time the getNextEntry() or repopulate() methods are
        called.
        
        \param[in] maskBases If this flag is \c true, then all
        lowercase bases are converted to 'N' rather than uppercase
        characters, causing them to be ignored by downstream
        processing.
        
        \param[in] randomizeNbases If this value is \c false the \c
        'N' characters are \b preserved.  However, if the flag is set
        to \c true then the 'N' characters are randomly converted to a
        regular base character from the list \c ATCG.

        \param[in] normalizeNTs If this value is \c false then the
        reads are left untouched and the previous 2 flags have no
        effect.
    */
    void setOptions(bool maskBases, bool randomizeNbases, bool normalizeNTs);

    /** Destructor.

        The destructor for the input file. Currently the destructor
        does not have any specific tasks to perform.
    */
    virtual ~InputFile();
    
protected:
    /** The constructor.

        The constructor merely initializes the instance variables to
        their default initial value.  Since this class is merely a
        base class (that is not meant to be directly instantiated) the
        constructor is protected.  Typically derived classes are
        instantiated via the InputFileFactory class.

        \param[in] fileName The path to the file from where the data
        for ESTs is being read.  This value is stored in the base
        class for future reference.
    */
    InputFile(const std::string& fileName);

    /** Obtain the current position from where the data for the next
        entry is going to be read.

        This method must be appopriately implemented by one of the
        derived classes to return the position from where the next
        entry is going to be read.  The value returned by this method
        is passed back in to the \c setCurrPos() method.  It is upto
        the derived class to return a suitable value that can be
        interpreted by these two methods.
        
        \note This method must be suitably implemented by the derived
        class.
        
        \return The current location from where the next entry is
        going to be read.

        \see setCurrPos()
    */
    virtual size_t getCurrPos() = 0;

    /** Set the current position from where the next entry is to be
        read.

        This method must be appropriately implemented by the derived
        classes to set the logical position from where the next entry
        is to going to be read.  It is upto the derived class to
        suitably implement this method.

        \param[in] position The position that was reported by the
        getCurrPos() implementation (and was saved by the InputFile)
        in an earlier call.

        \return This method returns \c true if the position for next
        entry to be read has been reset by this method. On errors or
        if the position could not be reset then this method returns
        false.

        \see getCurrPos()
    */
    virtual bool setCurrPos(const size_t position) = 0;
    
    /** Helper method to read information about a given entry.

        This method must be appropriately implemented by the derived
        class to read the information for the current entry.

        \param[out] info The header information about the EST to be
        read from the file.

        \param[out] sequence The nucleotide sequence information that
        has been read from the file.

        \param[out] quality The nucleotide quality information for
        each of the nucleotide bases in the sequence.

        \return This method returns true if the data was successfully
        read.  On errors it returns false.
    */
    virtual bool readEntry(std::string& info, std::string& sequence,
                           QualityVector& quality) = 0;
    
private:
    /** The full path to the file from where the data is being read
        for this input file.

        This instance variable is initialized by the constructor to
        the path of the file from where the EST data is being read.
        This value can be accessed via the getFileName() method.
     */
    const std::string fileName;

    /** Flag to indicate if lower case nucleotide base entries must be
        considered as mask bases.

        If this flag is \c true, then all lowercase \c "atcg" bases
		are converted to \c 'N'. Otherwise they are converted to
		uppercase letters. This value is set via the setOptions()
		method.
    */
    bool maskBases;

    /** Flag to indicate if \c 'N' characters in sequence reads must
        be converted randomly converted to one of \c "ATCG" values.

        If this value is \c false the \c 'N' characters are \b
        preserved.  However, if the flag is set to \c true then the
        'N' characters are randomly converted to a regular base
        character from the list \c ATCG.  This value is set via the
        setOptions() method.
    */
    bool randomizeNbases;

    /**
       If this flag is true then the reads are normalized to be only
       nucleotide sequences consisting of ATCG.  If this flag is
       false, then the this normalization is not done.
     */
    bool normalizeNTs;
    
    /** A map to maintain the offsets of each EST entry in this
        InputFile.

        This map is used to quickly look-up the file offset of a given
        EST entry. The ID of the EST is the key into this map.
        Entries are added each time the getNextEntry() method is
        called.  This array is used to repopulate information for a
        given EST entry.
    */
    std::map<int, size_t> entryOffset;

    /** A dummy operator=
            
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.

        \return Reference to this.
    */
    InputFile& operator=(const InputFile& src);
};

#endif
