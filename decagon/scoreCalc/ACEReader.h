#ifndef ACE_READER_H
#define ACE_READER_H

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

#include "HashMap.h"
#include "DecagonException.h"
#include "NtDistr.h"
#include <fstream>

/** \typedef HashMap< std::string, std::pair<int, int>, StringHasher > BSLookupMap

	\brief A typedef for a hash map whose key is std::string and
    contains pairs of integers as value.
    
    <p>The following typedef provides a short cut for using a hash map
    whose key is a std::string and contains pairs of integers. It is
    used to store Base Segment \c BS entries in an ACE file. The \c BS
    entry is in the following form:</p>

    <p><tt>BS 1 515 K26-572c</tt></p>

    <p>The last word, which is the name of the read is used as the key
    in the hash map. The start-offset and end-offset values are stored
    as the pair of integers. These values are used later on when the
    actual reads are loaded from the ACE file to create the necessary
    BS entries in contig.  The entries are created in a delayed manner
    as we need to know the index where a specific entry is being added
    to the list of reads prior to creating the entries.</p>
*/
typedef HashMap< std::string, std::pair<int, int>, StringHasher > BSLookupMap;

// Some forward declarations to keep compiler happy
class Contig;
class ESTList;

/** The class for reading Contigs from a given ACE text file.

    This class provides the necessary implementation to read contigs
    from an ACE file into memory for further analysis.  The ACE file
    is typically generated as the output of a genomic-assembler (such
    as: EAST, BATON, MiRA, etc.). The data read from the ACE file is
    stored in an assembler-neutral Contig object.  The contigs created
    from the file are stored into ContigList vector. Details on the
    ACE file format are available at
    http://www.animalgenome.org/bioinfo/resources/manuals/ace.filefmt.html
*/
class ACEReader {
public:    
    /** Constructor to create an ACEFile for reading data from an ACE
        text file.

        The constructor attempts to open the given file name for
        reading as an ACE file. It attempts to read the first \c AS
        header from the file so that the next read will commence with
        reading a contig from the \c CO header.

        \param[in] fileName The full path to the ACE file from where
        the contig data is to be read.

        \exception DecatonException This method throws an exception if
        the ACE file could not be opened or if the \c AS header was
        not successfully read.
    */
    ACEReader(const std::string& fileName) throw (DecagonException);

    /** Determine if the input file is still good or if it has
        encountered any errors.
        
        This method can be used to determine if this input file is
        still in "good" state or if an error has occured preventing it
        from operating correctly.
        
        \return This method returns \c true if the file is good. If an
        error had occured this method returns \c false.
    */
    inline bool good() const { return !aceFile.fail(); }

    /** Determine if another entry is pending to be read from this
        reader.

        This method can be used to determine if another <i>new</i>
        entry is pending to be read from this reader.  This method
        returns true, if the input file is not at EOF.

        \return This method returns \c true if another entry is
        pending to be read.  If all entries have been read then this
        method returns \c false.
    */
    inline bool hasNextEntry() const { return !aceFile.eof(); }

    /** Destructor.

        The destructor for this class.  The destructor merely closes
        the file handle that was opened in the constructor.
    */
    ~ACEReader();

    /** Read the next Contig from the ACE file.

        This is a top-level method that must be reapeatedly invoked
        (until this method returns \c false) to load all the contigs
        from this ACE reader.  This method reads all the information
        associated with the next contig and stores it in the given
        contig object.  Specifically this method operates in the
        following manner:

        <ol>

        <li>It checks and reads the next \c CO contig header.  If a \c
        CO header is not found then this method returns \c false.</li>

        <li>Next the consensus sequence for the contig is read (via
        call to readMultiLineData() helper method).</li>

        <li>Next the checkReadContigQuality() method is invoked to
        load any quality (\c BQ) reads for the contig.</li>

        <li>Next the checkReadAFBSEntries() method is used to read all
        the \c AF and \c BS entries for the current contig. </li>

        <li>Next this method calls repeatedly calls the loadRead()
        method (until it returns \c false) to load all the reads
        associated with the current contig.</li>

        <li>Finally the any readCtRtWaEntries() method is repeatedly
        invoked (until it returns false) to read and skip over all \c
        CT, \c RT, and \c WA entries (if any) associated with the
        contig.</li>

        <li>Lastly, this method reads any additional reads that maybe
        associated with a contig.  This operations is performed only
        because Mira seems to generate ACE file with contigs at this
        point (possibly be due to a bug in Mira?).</li>
        
        </ol>

        \param[out] contig The contig into which the next set of
        contig information in the ACE file is to be read and loaded.

        \param[out] reads The cDNA reads associated with the given
        contig.

        \return This method returns if the contig was successfully
        loaded from the ACE file. If a contig inforation was not found
        in the ACE file this method returns \c false (indicating end
        of the ACE file).

        \exception DecagonException This method throws a suitable
        exception whenever errors occur during contig reading.
    */
    bool readContig(Contig& contig, ESTList& reads) throw (DecagonException);

    /**
       A simple unit-testing method.

       This method is a simple unit testing method that can be used to
       test the functionality of this class.  This method creates a
       local ACEReader object and attempts to load a given ACE file.
       This method is pretty simple and is defined as:

       \code

       try {
           ACEReader ar(aceFileName);
           Contig contig;
           ESTList reads;
           while (ar.readContig(contig, reads));
       } catch (const DecagonException& de) {
           std::cerr << de << std::endl;
       }       
       
       \endcode

       \param[in] aceFileName The path to the ACE file to be used for
       testing.
     */
    static void testReader(const std::string& aceFileName);
    
protected:
    /** Helper method to read the contig header line from the ACE file.

        This is a refactored method (introduced to streamline the code
        in readContig() method) that is invoked from the readContig()
        method to read the contig header (\c CO) for the next contig.

        \return This method returns the contig ID from the CO header.

        \exception DecagonException This method throws an exception if
        the contig header was not successfully read from the file.
    */
    std::string readCOHeader() throw(DecagonException);

    /** Helper method to skip over blank lines and return the next
        non-blank line.

        ACE file format can have blank lines in it. However, the blank
        lines don't provide any useful information. Consequently, this
        method is used to conveniently skip over blank lines and
        return the next non-blank line form the ACE file.

        \param[out] line The string into which the next non-blank line
        (if any) is to be read.

        \return This method returns \c true if a valid non-blank line
        was successfully read. Otherwise this method returns \c false.
    */
    bool getNextLine(std::string& line);

    /** Helper method to read multi-line data from the ACE file.

        This is a convenience method that can be used to read multiple
        lines of data from the ACE file. Multiple lines of data may be
        used to store nucleotide sequences and quality data. This
        method repeatedly reads lines until a blank-line is
        encountered. Note that it first skips over any blank lines
        before it starts reading multi-line data.

        \param[out] data The string to contain multiple lines
        concatenated together. Line breaks and leading/trailing white
        spaces are not included in the string.

        \param[in] eol An optional string to be added to the end of
        each line of data read from the file. This is useful when
        newlines are meangiful in the data (like in quality scores).
        
        \return This method returns \c true if a valid non-blank line
        was successfully read. Otherwise this method returns \c false.        
    */
    bool readMultiLineData(std::string& data, const std::string& eol = "");

    /** Check and read quality values for the contig.

        This is a helper method that is used to read the quality
        values for a given contig from the ACE file, if available. The
        availability check for quality values is performed by reading
        the next non-blank line and checking if it is the \c BQ
        header. If the next non-blank line is not \c BQ then this
        method resets the file pointer (to the position prior to this
        method call) and exits immediately.  On the other hand, if the
        \c BQ header was read, then this method reads the quality
        values and sets the vlaues in the specified contig.

		\param[in,out] contig The contig into which the quality values
		are to be read and set. This method checks to ensure that the
		size of the consensus sequence for the contig matches the
		number of quality values read.

		\exception DecagonException This method throws an exception if
		the length of the conensus sequence for the contig does not
		match the number of quality values read from the ACE file.
    */
    void checkReadContigQuality(Contig& contig) throw (DecagonException);

    /** Helper method to read the set of \c AF and \c BS entries
        associated with a given contig.

        This is a ehlper method that is invoked from the readContig()
        method to read the set of \c AF and \c BS entries associated
        with a given contig.  These entries are of the form:

        <p>
        <tt><b>AF</b> <i>&lt;read name&gt;</i> <i>&lt;<b>C</b> or
<b>U</b>&gt;</i> <i>&lt;padded start consensus position&gt;</i></tt><br/>

        <tt><b>BS</b> <i>&lt;padded start consensus position&gt;</i> <i>&lt;padded end consensus position&gt;</i> <i>&lt;read name&gt;</i></tt>

        </p>

        \param[out] afEntries The hash map to which the \c AF entries
        from the ACE file are to be added.  These entries are used to
        look-up adidtional information about individual cDNA reads
        associated with a given contig.

        \param[out] bsEntries The hash map to wich the \c BS entries
        read from the ACE file are to be added.  These entries are
        used to look-up adidtional information about individual cDNA
        reads associated with a given contig.
    */
    void checkReadAFBSEntries(StringIntMap& afEntries, BSLookupMap& bsEntries)
        throw (DecagonException);

    /** Helper method to load a read (RD) entry associated with a
        given contig from the ACE file.

        <p>This method is a helper method that is repeatedly invoked
        from the readContig() method to load all the reads associated
        with a contig.  Note that a single contig is composed from one
        or more cDNA fragments/reads.  This method loads the next \c
        RD line from the ACE file. If the next line is not an \c RD
        entry, then this method undoes file read and returns \c false
        (indicating an \c RD entry was not found).</p>

        <p>However, if the next line happens to be an \c RD entry,
        then this method reads the components of the \c RD entry that
        is of the form:

        <p>
        <tt><b>RD</b> <i>&lt;read name&gt;</i> <i>&lt;# of padded
        bases&gt;</i> <i>&lt;# of whole read info items&gt;</i>
        <i>&lt;# of read tags&gt;</i></tt><br/>
        
        <i>Multi-line Nucleotide sequence for the read</i><br/>

        <tt><b>QA</b> <i>&lt;qual clipping start&gt;</i> <i>&lt;qual
        clipping end&gt;</i> <i>&lt;align clipping start&gt;</i>
        <i>&lt;align clipping end&gt;</i><br/>
        
        <tt><b>DS CHROMAT_FILE:</b> <i>&lt;name of chromat
        file&gt;</i> <b>PHD_FILE:</b> <i>&lt;name of phd file&gt;</i>
        <b>TIME:</b> <i>&lt;date/time of the phd file&gt;</i>
        <b>CHEM:</b> <i>&lt;prim, term, unknown, etc&gt;</i>
        <b>DYE:</b> <i>&lt;usually ET, big, etc&gt;</i>
        <b>TEMPLATE:</b> <i>&lt;template name&gt;</i>
        <b>DIRECTION:</b> <i>&lt;fwd or rev&gt;</i>
        
        </p>

        This method extracts the aforementioned information and adds
        it to the reads parameter passed ot this method. Next is
        extracts the \c AF entry for the given read (from the
        afEntries passed to this method) and updates the distribution
        of nucleotides in the \c ntDistr parameter passed to this
        method.

        It then loads the optional quality values from the \c QA entry
        associated with the read.  Currently the \c QA entry is
        discarded by this method. Next, it checks and reads the \c DS
        entry (if any).  Similar to the \c QA entry, the \c DS entry
        is not really used by this class.

        Finally, this method creates a suitable AlignmentInfo object
        and adds it to the supplied contig.
        
        \param[out] contig The contig to which the next read is to be
        added.

        \param[out] reads The list of reads associated with the
        current contig. The next read from the ACE file is added ot
        this list.

        \param[out] ntDistr The overall nucleotide distribution (that
        indicates the number of bases occurring at each nucleotide
        position in the final consensus contig) that is being
        currently read and is to be updated by this method.

        \param[in] afEntries The \c AF entries that provide the
        initial starting position (in number of nucleotides) for the
        reads associated with the current contig being loaded.

        \exception DecagonException This method exposes any exceptions
        that may occur when reading data from the ACE file.

        \return This method returns \c true to indicate that a read
        for a given contig was successfully read and another call to
        this method is necessary.  If a read was not found then this
        method returns \c false.
    */
	bool loadRead(Contig& contig, ESTList& reads, NtDistrList& ntDistr,
                  StringIntMap& afEntries)
        throw (DecagonException);

    /**
       Helper method to read and ignore \c CT, \c RT, and \c WA
       entries (for a contig) in an ACE file.

       This method is a simple helper method that reads a \c CT, \c
       RT, or a \c WA header entry. This method was introduced because
       these headers have nested scopes delimited by \c '{' and \c '}'
       characters. This method reads one or more lines from the ACE
       file until a full header is successfully read.

       \exception DecagonException This method throws an eception when
       an error occurs when attempting to read a given header entry.

       \return This method returns \c false if a \c CT, \c RT, or a \c
       WA entry was not read.  Otherwise this method returns \c true
       indicating a header was successfully read and another call to
       this method is needed to read any remaining headers.
    */
    bool readCtRtWaEntries() throw (DecagonException);
    
private:
    /** The OS-level file pointer used to read data from the FASTA
        file.
        
        This instance variable that used to hold the file pointer from
        where the FASTA data is being read.  This instance variable is
        initialized in the constructor once the FASTA file has been
        successfully opened for reading. The file is closed in the
        destructor.
    */
    std::ifstream aceFile;
};

#endif
