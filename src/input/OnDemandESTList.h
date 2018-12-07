#ifndef ON_DEMAND_EST_LIST_H
#define ON_DEMAND_EST_LIST_H

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

#include "ESTList.h"

// Some forward declarations to keep the compiler happy and fast
class InputFile;

/** A list of cDNA fragments that are loaded on demand.

    This class extends the shared ESTList class and provides the
    following additional functionality:

    <ul>

    <li>The class can automatically load entries from a given
    InputFile (a base class for various input files such as: FASTA and
    SFF). </li>
    
    <li>When data is read from an InputFile, the the headers and
    character sequence information for the cDNA fragments are cleared
    and then reloaded on-demand. This feature strives to reduce the
    peak memory-footprint for each process and thereby enabling the
    processing of large sets of data files in a distributed
    manner.</li>

    </ul>
*/
class OnDemandESTList : public ESTList {
public:
    /** The default constructor.
        
        The default constructor creates an empty list without any EST
        entries in it. Entries can be added to the list via the add()
        methods in the base class or loaded from an InputFile via the
        add method provided by this class.
    */
    OnDemandESTList();

    /** The destructor.

        The destructor is typically called only once when the process
        terminates and the shared list of ESTs is to be deleted.  The
        destructor calls the reset() method to clear out all the
        contents of this list and closes any open InputFile objects
        (maintained by this class to load cDNA information on-demand).
    */
    ~OnDemandESTList();

    /** Obtain a mmutable pointer to an EST entry in this list.

        This method overrides the default implementation in the base
        class to check and load the necessary cDNA information (if it
        is not already loaded).

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; OnDemandESTList::size().

        \note This method also loads the necessary header information
        and nucleotide sequence from an input file if
        needed. Consequently, calling this method for non-cached
        entries can be time consuming.  To get whatever information is
        in the cache, use the overloaded get(const int, const bool)
        method.

        \return A pointer to the EST entry at the specified index. The
        return value for an invalid index is NULL.
    */
    EST* get(const int index, const bool populate = false) {
        // If EST is populated return it directly. Otherwise
        // repopulate entry and then return it.
        return ((estVector[index]->isPopulated() || (!populate)) ?
                estVector[index] :
                repopulate(index));
    }

    /** Interface method to add entries ESTList from a given input
        file.

	This method is typically used by the InputFileFactory to load
	data files specified by the user.  This method assumes that
	the input file has been suitably configured to handle masking
	of bases and handling of \c 'N' nucleotide entries.
		
        \param[in] inputFile Reference to a valid InputFile object
        from where the cDNA fragments are to be loaded into this list.

        \param[in] startIndex The starting index position of the EST
	from where the full data is to be retained. For other entries,
	the header information and nucleotide sequences are loaded
	on-demand.  This is done to reduce peak memory footprint
	(thereby enabling the processing of large data files). This
	value must be less than endIndex.  This value with reference
	to this list of ESTs (not just relative to this file).
        
        \param[in] endIndex The ending index position of the EST upto
	(and <b>not</b> including) which the full data is to be
	retained. For other entries, the header information and
	nucleotide sequences are loaded on-demand.  This value must be
	greater than startIndex. This value is with-reference-to this
	list of cDNA fragments.  By default the startIndex and
	endIndex are set such that all cDNA fragments are loaded
	on-demand.

        \note This class maintains an open handle to the file (until
        reset() method is called).  
        
        \return On successful creation, this method returns \c
        true. The reference to the global list can be obtained via a
        call to the get() method.
    */
    bool add(InputFile* inputFile, const long startIndex = 0x7ffffffL,
             const long endIndex = 0x7ffffffL);
    
    /** Clears out all EST entries, open input files, and resets
        internal values to default initial value.

        This method can be used to reset the list to an empty list
        without any EST entries in it. This method also closes any
        open input file. In addition, it resets all the internal
        counters and variables to their default initial value(s).
    */
    void reset();

    /** Repopulate a given EST entry in this list.

        This method can be used to repopulate information (that was
        earlier unpopulated to save memory footprint) for a given EST.
        The EST data is reloaded from the appropriate InputFile.

        \param[in] index \param index The index of the entry for which
        a pointer is to be returned by this method.  This index value
        must be in the range 0 &le; ESTList::size().

        \return A pointer to the EST entry at the specified index. The
        return value for an invalid index is NULL.
    */
    EST* repopulate(int index) const;

    /** Repopulate a given EST entry in this list.

        This method can be used to repopulate information (that was
        earlier unpopulated to save memory footprint) for a given EST.
        The EST data is reloaded from the appropriate InputFile.

        \param[in] est The entry for which a pointer is to be returned
        by this method.  If this entry is already populated, then this
        method performs no operations.
    */
    virtual void repopulate(const EST* est) const;
    
    /** Obtain the information associated with a given EST (load it if
        not available).
        
        This method returns the name and other information associated
        with the EST.  This information is typically the first header
        line read from a FASTA file. If the information is not
        available then just the information is loaded from

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; ESTList::size().
        
        \return Any information available for this EST. Return the
        information
    */
    virtual std::string getESTInfo(const int index) const;

    /** Method to clear general information and sequence data.

        This method can be used to unpopulate the FASTA header and
        actual sequence (base pairs) information a given entry.  This
        frees up memory allocated to hold this data thereby minimizing
        the memory footprint for this EST.  This enables holding a
        large number of skeleton EST's in memory.

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; ESTList::size().
    */
    virtual void unpopulate(int index) const {
	estVector[index]->unpopulate();
    }
    
protected:
    /** The input source from where some (or all) of the EST data has
        been read.

        This instance variable is used to track the file from where
        the data for the ESTs has been read. This pointer is
        maintained so that unpopulated EST data can be repopulated on
        demand.  The pointer is initialized to NULL.  A new InputFile
        object is created by the add(const std::string&, const long,
        const long) method and added to this vector once the data from
        this file has been successfully loaded.
    */
    std::vector<InputFile*> inputFileList;

private:    
    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.

        \return Reference to this.
    */
    OnDemandESTList& operator=(const OnDemandESTList& src);
};

#endif
