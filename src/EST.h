#ifndef EST_H
#define EST_H

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
#include <functional>
#include <memory>
#include "ESTCustomData.h"
#include "Utilities.h"

/** A single EST.

    This class is used to represent a single EST.  An EST object
    instance consists of the following information:

    <ul>

    <li>\b id: A unique identifier (usually a number) for this
    EST.</li>
    
    <li>\b info: The name and other information associated with the
    EST.  This information is typically the first header line read
    from a FASTA file.</li>

    <li>\b sequence: The actual sequence of base pairs associated with
    this EST.<li>

    <li>\b offset: The offset of in the FASTA file from where this EST
    was read.  This information can be used to conditionally and
    rapidly load ESTs from a file.</li>

    </ul>
    
*/
class EST {
public:
    /** EST Constructor.

        This constructor is used to instantiate an EST method.
	
        \param[in] id The unqiue ID value to be set for this EST.
        
        \param[in] info The name and other information associated with
        the EST.  This information is typically the first header line
        read from a FASTA file.  This information can be NULL.
        
        \param[in] sequence The actual sequence of base pairs
        associated with this EST.  The sequence information that must
        be used to create this EST.  The sequence information can be
        NULL.
        
        \param[in] offset The offset of in the FASTA file from where
        this EST was read.  This information can be used to conditionally
        and rapidly load EST's from a file.        
    */
    EST(const int id, const char *info,
        const char* sequence = NULL, const int offset = -1);
  
    /** Create a valid EST.
        
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
        read from a FASTA file.  This information can be NULL.
        
        \param[in] sequence The actual sequence of base pairs
        associated with this EST.  The sequence information that must
        be used to create this EST.  The sequence information can be
        NULL.
        
        \param[in] offset The offset of in the FASTA file from where
        this EST was read.  This information can be used to conditionally
        and rapidly load EST's from a file.

        \return If the id is valid and a duplicate EST with the same
        ID is not present, then this method creates a new EST and
        returns a pointer to that EST back to the caller.
    */
    static EST* create(const int id, const char *info,
                       const char* sequence = NULL,
                       const long offset = -1);

    /** Loads data from a FASTA file to create an EST.

        This method provides a convenient interface for loading
        information regarding an EST from a given FASTA file and using
        the information to create either a fully populated or
        partially populated EST.

        \param[in,out] fastaFile The FASTA file from where the EST data
        is to be currently loaded.  If this pointer is NULL then this
        method perform no action and returns immediately with NULL.

        \param[in,out] lineNum A line number counter to be updated to
        provide the user with a more meaningful error message.

        \note At the end of this method the fastaFile's file pointer
        will point at the beginning of the next EST (if any) in the
        file.
    */
    static EST* create(FILE* fastaFile, int& lineNum);

    /** Obtain the list of ESTs.

        This method may be used to obtain a reference to the list of
        ESTs currently defined.

        \return The list of ESTs currently defined.
    */
    static std::vector<EST*>& getESTList() { return estList; }

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
	static int getProcessedESTCount();
	
    /** Obtain the number of ESTs in this list.

        This method may be used to determine the number of ESTs that
        have been defined and added to this list.

        \return The number of ESTs currently defined.
    */
    static int getESTCount() { return (int) estList.size(); }

    /** Helper method to determine the longest EST.

        This method can be used to determine the length of the longest
        EST loaded thus far.  This information is typically used to
        allocate buffers and other data structures for analysis.

        \note This method computes the length of the longest EST the
        first time it is invoked. Consequently, it should be called
        only after all the ESTs have been loaded.

        \return The length of the longest EST to be processed.
    */
    static size_t getMaxESTLen();
    
    /** Obtain a given EST from the EST list.

        This method is a convenience method that can be used to obtain
        a given EST from the list of ESTs.

        \param[in] estIdx The zero-based index of the EST that is
        desired from the list of ESTs in this class.  If this index is
        invalid then the behavior of this method is undefined.

        \return A mutable pointer to the EST at the provided estIdx
        index position in the EST list.
    */
    static EST* getEST(const int estIdx) { return estList[estIdx]; }
    
    /** Dump currently loaded ESTs in FASTA format.

        This method can be used to dump the currently loaded EST's in
        FASTA file format to a given output stream.

        \param[out] os The output stream to which EST data is to be
        dumped.
    */
    static void dumpESTList(std::ostream& os);

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
    static void dumpESTList(std::ostream& os, const bool processed);

    /** Dump this EST information in FASTA format.

        This method can be used to dump the information associated
        with the EST in FASTA format to a given output stream.

        \param[in] os The output stream to which the EST's information
        must be written in FASTA format.
    */
    void dumpEST(std::ostream& os);
    
    /** Delete and clear all ESTs.

        This method can be used to delete and clear all the EST's from
        the internal list of EST's currently loaded.
    */
    static void deleteAllESTs();

	/** Delete and clear out the last EST in the list.

		This method can be used to delete the last EST in the
		list. This method rests the maximum EST length instance
		variable as needed.  This method is typically used to remove
		dummy ESTs that are added to the end of the list by some
		filters.

		\param[in] count The number of ESTs to be removed from the
		list.
	*/
	static void deleteLastESTs(const int count);
	
    /** Obtain the ID of this EST.
        
        \return The ID of the EST that was set when this EST was
        created.
    */
    inline int getID() const { return id; }
    
    /** Obtain the information associated with this EST.
        
        The name and other information associated with the EST.  This
        information is typically the first header line read from a
        FASTA file.  This information can be NULL if the EST is only
        partially loaded.

        \return Any information available for this EST.
    */
    inline const char* getInfo() const { return info; }

    /** Obtain the actual sequence of base pairs for this EST.

        Note that sequence ifnoramtion for an EST can be null if itis
        only partially loaded from a file.  Entries are parially
        loaded to reduce memory foot print when processing large data
        sets.
        
        \return The actual sequence of base paris for this EST if
        available. Otherwise this method returns NULL.
    */
    inline const char* getSequence() const { return sequence; }

    /** Obtain the similarity metric for this EST.

        The similarity metric is a quantitative representation of the
        similarity between two ESTs.  The similarity metric is
        generated during analysis when one EST is compared with
        another.  The similarity value is initialized to -1.

        \return The similarity metric for this EST.
    */
    inline float getSimilarity() const { return similarity; }

    /** Set the similarity metric for an EST.

        This method must be used to change the similarity metric for
        this EST.  The similarity metric is a quantitative
        representation of the similarity between two ESTs.  The
        similarity metric is generated during analysis when one EST is
        compared with another.

        \param[in] sim The similarity metric value to which this EST's
        similarity much be changed.
    */
    inline void setSimilarity(const float sim) { similarity = sim; }
    
    /** Method to clear general information and sequence data.

        This method can be used to unpopulate the FASTA header and
        actual sequence (base pairs) information from this EST.  This
        frees up memory allocated to hold this data thereby minimizing
        the memory footprint for this EST.  This enables holding a
        large number of skeleton EST's in memory.
    */
    void unpopulate();

    /** Repopulate necessary information from a given fasta file.

        This method can be used to request an EST to repopulate its
        FASTA header and actual sequence (base pair) information from
        a given FASTA file.  This method uses the offset (saved when
        this EST was originally loaded) to load the information from
        the file.

        \param[in,out] fastaFile The file from where the EST
        information is to be loaded.  If the file changes during EST
        analysis the behavior of this method is undefined.

        \return This method returns true if the repopulating the data
        was successfully completed. On errors this method returns
        false.
    */
    bool repopulate(FILE *fastaFile);

    /** Change the custom data associated with this EST.

        This method can be used to change (or set) the custom data
        associated with this EST.  Note that any earlier custom data
        associated with this EST is lost (and deleted if necessary by
        auto_ptr) before the new value is set.

        \param[in,out] src The new custom data to be set for this EST.
        After this call, this EST owns the data referred by src.
    */
    void setCustomData(std::auto_ptr<ESTCustomData>& src) { customData = src; }

    /** Change the custom data associated with this EST.

        This method can be used to change (or set) the custom data
        associated with this EST.  Note that any earlier custom data
        associated with this EST is lost (and deleted if necessary by
        auto_ptr) before the new value is set.

        \param[in,out] src The new custom data to be set for this EST.
        After this call, this EST owns the data referred by src.
    */    
    void setCustomData(ESTCustomData* src)
    { customData = std::auto_ptr<ESTCustomData>(src); }


    /** Obtain a mutable reference to custom data associated with this EST.

        This method can be used to obtain a mutable reference to the
        custom data associated with this EST.  This method essentially
        returns the custom value set by the last successful call to
        one of the polymorphic setCustomData() methods in this class.
        By default this method returns NULL.

        \note The custom data set in this class is returned as a
        auto_ptr.

        \return The custom data (if any) associated with this EST.
    */
    inline std::auto_ptr<ESTCustomData>& getCustomData() { return customData; }

    /** Obtain an immutable reference to custom data associated with
        this EST.

        This method can be used to obtain an immutable reference to
        the custom data associated with this EST.  This method
        essentially returns the custom value set by the last
        successful call to one of the polymorphic setCustomData()
        methods in this class.  By default this method returns NULL.

        \note The custom data set in this class is returned as a
        auto_ptr.

        \return The custom data (if any) associated with this EST.
    */    
    inline const std::auto_ptr<ESTCustomData>& getCustomData() const
    { return customData; }
    
    /** The destructor.
        
        The destructor for the EST essentially releases the memory
        used to hold the information and sequence data for a given EST.
    */
    ~EST();

    /** Functor for EST sorting.

        This Functor is used when sorting ESTs based on similarity
        metric at the end of analysis just prior to generating the
        final report.
    */
    struct LessEST : public std::binary_function<EST, EST, bool> {
        inline bool operator()(const EST* x, const EST* y) {
            return (x->similarity > y->similarity);
        }
    };

    /** Helper method to read a line from a given file.
		
		This is a helper method that can be used to read a long line
		from a given file.
		
		\param[in] fp The file from where the line is to be read.
		
		\return The string read from the file.
    */
    static std::string getLine(FILE *fp);

    /** Determine if this EST has already been processed.

        This method exposes a generic flag that is provided as a
        convenience for algorithms to mark if this EST has gone
        through their processing.

        \return This method returns \c true if this EST has been
        flagged as having been processed. Otherwise this method
        returns \c false.
     */
    inline bool hasBeenProcessed() const { return processed; }

    /** Set if this EST has already been processed.

        This method provides a generic flag as a convenience for
        algorithms to mark if this EST has gone through their
        processing. By default ESTs are marked has processed when they
        are instantiated.
        
        \param[in] processedFlag If this flag is \c true then this EST
        is flagged as having been processed. If this flag is \c false
        then the EST is flagged as not-processed (and requiring
        processing).
    */
    inline void setProcessed(const bool processedFlag)
    { processed = processedFlag; }
    
protected:
    /** The unique ID for this EST.

        This member holds the unique ID for this EST.  The ID is set
        when the EST is instantiated and is never changed during the
        life time of this EST.  The id is used to access and extract
        EST information.
    */
    const int id;

    /** The name and other information associated with the EST.  This
        information is typically the first header line read from a
        FASTA file.  The information may be dynamically loaded on
        demand to reduce memory footprint when processing large data
        sets.
    */
    char *info;

    /** The actual sequence of base pairs associated with this EST.
        This information is typically read from a FASTA file.  The
        information may be dynamically loaded on demand to reduce
        memory footprint when processing large data sets.
    */
    char *sequence;

    /** The offset in the FASTA file to load the data from.

        The offset of in the FASTA file from where this EST was read.
        This information can be used to conditionally and rapidly load
        EST's from a file.  This value is initialized when the EST is
        insantiated and is never changed during the life time of an
        object.
    */
    const long offset;

    /** A similarity value for this EST with respect to another EST.

        This instance variable is used to hold a similarity metric for
        this EST.  The similarity metric is generated during analysis
        when one EST is compared with another.  The similarity value
        is initialized to -1.  It is accessed via the getSimilarity()
        method and changed via the setSimilarity() method.
    */
    float similarity;

    /** Instance variable to track if EST has gone through some
        processing.

        This is a generic flag that is provided as a convenience for
        algorithms to mark if this EST has gone through their
        processing.  By default this instance variable is intialized
        to \c false. Once it has been processed, the setProcessed()
        method can be used to set/reset this flag. The
        hasBeenProcessed() method can be used to determine if this EST
        has already been processed.
    */
    bool processed;

    /** Size of the longest EST.

        This static instance variable is used to track the size (in
        number of nucleotides) of the longest EST ever
        instantiated. The size of the longest EST can be used by
        algorithms to optimally allocate memory for processing ESTs.
    */
    static size_t maxESTlen;
    
    /** Place holder for some other custom data.
        
        This pointer acts as a convenient place holder for other
        classes to associate some uninterpreted user data (or data
        structure).  This member is initialized to NULL in the
        constructor.  Note that this pointer is managed using an
        auto_ptr that automatically deletes the data when the auto_ptr
        loses ownership of the data object.
    */
    std::auto_ptr<ESTCustomData> customData;
    
private:    
    /** The default constructor.

        The default constructor has been made private to ensure that
        EST's are never directly created.  Instead, a valid EST must
        be created using other constructors.
    */
    EST();

    /** A utility method to duplicate a c-string.

        This msethod is a simple utililty method that can be used to
        duplicate a given C-string.  This method uses the stsandard
        C++ new operator to duplicate the given C-string.

        \return This method simply returns NULL if src is
        NULL. Otherwise this method returns a pointer to a duplicate
        version of the specified string.
    */
    static char* duplicate(const char *src);
    
    /** The list of EST's currently being used.

        This list contains the complete set of ESTs that are currently
        defined.  This list includes partially loaded ESTs as well.
        New entries are added to the list by the create method.
    */
    static std::vector<EST*> estList;

    /** A dummy operator=

        The operator=() is supressed for this class as it has constant members
        whose value is set when the object is created.  These values cannot be
        changed during the lifetime of this object.

        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.

        \return Reference to this.
    */
    EST& operator=(const EST& src);
};

#endif
