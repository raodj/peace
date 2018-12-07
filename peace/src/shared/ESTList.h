#ifndef EST_LIST_H
#define EST_LIST_H

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

#include <vector>
#include <limits>
#include <functional>

#include "Utilities.h"
#include "EST.h"

// Forward declarations
class InputFile;

/** A list of cDNA fragments.

    This class is essentially used to encapsulate a list of ESTs that
    are being processed by the various components constituting
    PEACE. In addition, this class provides the following features:

    <ul>

    <li>It provides a shared data structure for sharing cDNA fragments
    between multiple sub-systems constituting PEACE.</li>

    <li>It provides a simple list that can contain fragments from
    multiple files.</li>
	
    <li>The class serves as a base class for OnDemandESTList (in the
    InputSubSystem) that can automatically load entries from a given
    InputFile (InputFile is the base class for various input files
    such as: FASTA and SFF). </li>

    </ul>

    \note Currently, we assume that each PEACE job will be used to
    process a single list of ESTs.  However, the fragments in this
    list can be loaded from multiple files or be directly added by the
    user.
*/
class ESTList {
public:
    /** The default constructor.
        
        The default constructor creates an empty list without any cDNA
        fragments in it.  Additional entries can be added using the
        add() method in this class.
    */
    ESTList();
    
    /** The destructor.
        
        The destructor is typically called only once when the process
        terminates and the shared list of ESTs is to be deleted.  The
        destructor calls the reset() method to clear out all the
        contents of this list.
    */
    virtual ~ESTList();
  
    /** Obtain a immutable pointer to an EST entry in this list.

        This method must be used to obtain a valid pointer to a given
        entry in this list.  This method simply returns the EST entry
        without any further checks to verify if the entry is fully
        populated.  EST entries may have been unpopulated (so
        nucleotide sequence information, FASTA header, and quality
        vector will be empty) to reduce memory space.  For many
        operations, the standard entries in an EST object will be
        sufficient.  However, if complete, fully populated entries are
        needed then the overloaded ESTList::get(const int, const bool)
        method must be used.

        \note The EST object returned by this method may be
        unpopulated and may not have all the data filled-in.

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; ESTList::size().

        \return A pointer to the EST entry at the specified index. The
        return value for an invalid index is NULL.
    */
    inline const EST* get(const int index) const {
        // Simply return the entry (derived class checks to repopulate)
        return estVector[index];
    }

    /** Obtain a mutable pointer to an EST entry in this list.

        This method must be used to obtain a valid pointer to a given
        entry in this list.

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; ESTList::size().

	\param repopulate If this flag is \c true and the requested
	EST is not fully populated (see EST::isPopulated()) then the
	EST entry is populated with all the information (such as
	nucleotide sequence, quality vector, and header GI
	information) before it is returned by this method.
		
        \note The method in this class does not repopulate entries
        because this class does not unpoplate EST entries to begin
        with. However, This method may be overridden in derived
        classes to load the necessary header information and
        nucleotide sequence from an input file if
        needed. Consequently, calling this method for non-cached
        entries can be time consuming.  To get whatever information is
        in the cache, use the overloaded get(const int) method(s) or
        operator[].

        \return A pointer to the EST entry at the specified index. The
        return value for an invalid index is NULL.
    */
    virtual EST* get(const int index,
                     const bool UNREFERENCED_PARAMETER(populate) = false) {
        return estVector[index];
    }
    
    /** Obtain an immutable (constant object) pointer to an EST entry
	in this list.

        This is a convenience overload of operator[] that can be used
        to obtain a valid pointer to a given entry in this list.  This
        method simply returns the EST entry without any further checks
        to verify if the entry is fully populated.  EST entries may
        have been unpopulated (so nucleotide sequence information,
        FASTA header, and quality vector will be empty) to reduce
        memory space.  For many operations, the standard entries in an
        EST object will be sufficient.  However, if complete, fully
        populated entries are needed then the overloaded
        ESTList::get(const int, const bool) method must be used.

        \note The EST object returned by this method may be
        unpopulated and may not have all the data filled-in.

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; ESTList::size().

        \note This method merely calls ESTList::get(const int). Refer
        to the documentation on ESTList::get(const int) for additional
        details.
        
        \return A pointer to the EST entry at the specified index. The
        return value for an invalid index is NULL.
    */    
    inline const EST* operator[](const int index) const
    { return get(index); }

    /** Obtain a mutable (non-constant object) pointer to an EST
	entry in this list.

        This is a convenience overload of operator[] that can be used
        to obtain a valid pointer to a given entry in this list.  This
        method simply returns the EST entry without any further checks
        to verify if the entry is fully populated.  EST entries may
        have been unpopulated (so nucleotide sequence information,
        FASTA header, and quality vector will be empty) to reduce
        memory space.  For many operations, the standard entries in an
        EST object will be sufficient.  However, if complete, fully
        populated entries are needed then the overloaded
        ESTList::get(const int, const bool) method must be used.

        \note The EST object returned by this method may be
        unpopulated and may not have all the data filled-in.

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; ESTList::size().

        \note This method merely calls ESTList::get(const int, const
        bool) [with second parameter being set to \c false to ensure
        that time is not spent on repopulating the object]. Refer to
        the documentation on ESTList::get(const int, const bool) for
        additional details.
        
        \return A pointer to the EST entry at the specified index. The
        return value for an invalid index is NULL.
    */    
    inline EST* operator[](const int index)
    { return get(index, false); }
	
    /** Determine the number of EST entries in this list.
        
        This method can be used to determine the current number of
        ESTs that have been added to this list.

        \return The number of ESTs in this list. If the list is empty
        then this method returns 0.
    */
    inline int size() const { return estVector.size(); }
    
    /** Create and add a valid EST (without quality).
        
        This method must be used to create a valid EST in the system.
        The information required to create the EST must be passed in
        as the parameter.  The EST names are expected to be unique in
        a given file.

        \note If the new EST is successfully instantiated, then this
        method adds the newly created EST to the end of the list of
        ESTs maintianed by this class.  Consequenlty, the parameter \c
        id must be equal to estList.size().
	
        \param[in] id The unqiue ID value to be set for this EST.
        
        \param[in] info The name and other information associated with
        the EST.  This information is typically the first header line
        read from a FASTA file.  This information can be an empty
        string (\c "").
        
        \param[in] sequence The actual sequence of base pairs
        associated with this EST.  The sequence information that must
        be used to create this EST.  The sequence information can be
        an empty string (\c "").  Any nucleotide masking and base
        normalization must be performed pior to invoking this
        method. See InputFile::normalizeBases() method for details.
        
        \return If the id is valid and a duplicate EST with the same
        ID is not present, then this method creates a new EST and
        returns a pointer to that EST back to the caller.
    */
    EST* add(const int id, const std::string& info,
             const std::string& sequence);

    /** Create and add a valid EST (with quality)
        
        This method must be used to create a valid EST in the system.
        The information required to create the EST must be passed in
        as the parameter.  The EST names are expected to be unique in
        a given file.

        \note If the new EST is successfully instantiated, then this
        method adds the newly created EST to the end of the list of
        ESTs maintianed by this class.  Consequenlty, the parameter \c
        id must be equal to estList.size().
	
        \param[in] id The unqiue ID value to be set for this EST.
        
        \param[in] info The name and other information associated with
        the EST.  This information is typically the first header line
        read from a FASTA file.  This information can be an empty
        string (\c "").
        
        \param[in] sequence The actual sequence of base pairs
        associated with this EST.  The sequence information that must
        be used to create this EST.  The sequence information can be
        an empty string (\c "").  Any nucleotide masking and base
        normalization must be performed pior to invoking this
        method. See InputFile::normalizeBases() method for details.

        \param[in] quality The quality vector to be set for this
        entry.  The number of entries in this vector must match the
        length of the nucleotide sequence parameter.
        
        \return If the id is valid and a duplicate EST with the same
        ID is not present, then this method creates a new EST and
        returns a pointer to that EST back to the caller.
    */
    EST* add(const int id, const std::string& info,
             const std::string& sequence, const QualityVector& quality);

    /** Interface method to add entries ESTList from a given input
        file.

        This method is typically used to load daa files specified by
        the user.  This method assumes that the input file has been
        suitably configured to handle masking of bases and handling of
        \c 'N' nucleotide entries.
		
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

        \note This method is actually implemented in a derived
        class. The base class simply returns \c false indicating
        failure.
        
        \return On successful creation, this method returns \c
        true. The reference to the global list can be obtained via a
        call to the overloaded get() method(s).
    */
    virtual bool add(InputFile* inputFile, const long startIndex = 0x7ffffffL,
		     const long endIndex = 0x7ffffffL);
	
    /** Obtain count of ESTs that have been flagged as being processed.
        
        This method can be used to determine the number of ESTs that
        have been flagged as being processed. Subtracting this number
        from the total number of ESTs indicates the number of ESTs to
        be processed.
        
        \note This method iterates over the list of ESTs to determine
        the current number of processed ESTs. So use this method
        sparingly.
	    
        \return The number of ESTs that have been flagged as having
        been processed.
    */
    int getProcessedESTCount() const;

    /** Helper method to determine the longest EST.

        This method can be used to determine the length of the longest
        EST loaded thus far.  This information is typically used to
        allocate buffers and other data structures for analysis.

        \note This method computes the length of the longest EST the
        first time it is invoked. Consequently, it should be called
        only after all the ESTs have been loaded.

        \return The length of the longest EST to be processed.
    */
    size_t getMaxESTLen() const;
    
    /** Dump currently loaded ESTs in FASTA format.

        This method can be used to dump the currently loaded EST's in
        FASTA file format to a given output stream.

        \param[out] os The output stream to which EST data is to be
        dumped.
    */
    void dumpESTList(std::ostream& os) const;

    /** Dump currently loaded and (un)processed ESTs in FASTA format.

        This method can be used to dump the currently loaded EST's in
        FASTA file format to a given output stream.

        \param[out] os The output stream to which EST data is to be
        dumped.

        \param[in] processed If this flag is \c true, then this method
        dumps only those ESTs that have been flagged as having been
        processed. If this flag is \c false, then this method dumps
        only un-processed ESTs.
    */
    void dumpESTList(std::ostream& os, const bool processed) const;

    /** Delete and clear out the last EST in the list.
	    
        This method can be used to delete the last EST in the
        list. This method rests the maximum EST length instance
        variable as needed.  This method is typically used to remove
        dummy ESTs that are added to the end of the list by some
        filters.
        
        \param[in] count The number of ESTs to be removed from the
        list.
    */
    void deleteLastESTs(const int count);

    /** Clears out all EST entries, open input files, and resets
        internal values to default initial value.

        This method can be used to reset the list to an empty list
        without any EST entries in it. This method also closes any
        open input file. In addition, it resets all the internal
        counters and variables to their default initial value(s).
    */
    virtual void reset();

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
    virtual void unpopulate(int UNUSED(index)) const {
	// Intentionally blank -- overridden in derived class(es).
    }
    
    /** Repopulate a given EST entry in this list.

        This method is overridden in derived class(es) to repopulate
        information (that was earlier unpopulated to save memory
        footprint) for a given EST.  Typically, the EST data is
        reloaded from the appropriate InputFile.

        \note The base class implementation does not perform any
        special task as it assumes that all the necessary information
        is available in memory.  It simply returns the entry at the
        given index.

        \param[in] index \param index The index of the entry for which
        a pointer is to be returned by this method.  This index value
        must be in the range 0 &le; ESTList::size().

        \return A pointer to the EST entry at the specified index. The
        return value for an invalid index is NULL.
    */
    virtual const EST* repopulate(int index) const;

    /** Repopulate a given EST entry in this list.

        This method is overridden in derived class(es) to repopulate
        information (that was earlier unpopulated to save memory
        footprint) for a given EST.  Typically, the EST data is
        reloaded from the appropriate InputFile.

        \note The base class implementation does not perform any
        special task as it assumes that all the necessary information
        is available in memory.

        \param[in] est The entry to be repopulated.  This entry must
        be a valid (non-NULL) entry obtained from this list via
        suitable API method(s).
    */
    virtual void repopulate(const EST* UNUSED(est)) const {}
    
    /** Obtain the information associated with a given EST (load it if
        not available).
        
        This method returns the name and other information associated
        with the EST.  This information is typically the first header
        line read from a FASTA file. If the information is not
        available then just the information is loaded from

        \param index The index of the entry for which a pointer is to
        be returned by this method.  This index value must be in the
        range 0 &le; ESTList::size().
        
        \note The base class implementation does not perform any
        special task as it assumes that all the necessary information
        is available in memory.  It simply returns the entry at the
        given index.
        
        \return Any information available for this EST. Return the
        information
    */
    virtual std::string getESTInfo(const int index) const {
        return get(index)->getInfo();
    }
    
protected:
    /** The list of EST's currently being used.

        This list contains the complete set of ESTs that are currently
        defined.  This list includes partially loaded ESTs as well.
        New entries are added to the list by the add() method.
    */
    std::vector<EST*> estVector;

    /** Instance variable to track the longest EST in this list.

        This instance variable is used to track the length of the
        longest EST in this list.  This value is initialized to zero.
        It is computed once when the getMaxESTLen() method is called.
        After that, the various overloaded add() methods in this class
        update this value when new entries are added.
    */
    mutable int maxESTLen;

private:    
    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.

        \return Reference to this.
    */
    ESTList& operator=(const ESTList& src);
};

#endif
