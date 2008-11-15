#ifndef EST_H
#define EST_H

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

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

        \param[inout] fastaFile The FASTA file from where the EST data
        is to be currently loaded.  If this pointer is NULL then this
        method perform no action and returns immediately with NULL.

        \param[inout] lineNum A line number counter to be updated to
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
    
    /** Delete and clear all ESTs.

        This method can be used to delete and clear all the EST's from
        the internal list of EST's currently loaded.
    */
    static void deleteAllESTs();
    
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

	\param[inout] fastaFile The file from where the EST
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

        \param[inout] src The new custom data to be set for this EST.
        After this call, this EST owns the data referred by src.
    */
    void setCustomData(std::auto_ptr<ESTCustomData>& src) { customData = src; }

    /** Change the custom data associated with this EST.

        This method can be used to change (or set) the custom data
        associated with this EST.  Note that any earlier custom data
        associated with this EST is lost (and deleted if necessary by
        auto_ptr) before the new value is set.

        \param[inout] src The new custom data to be set for this EST.
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
    /** EST Constructor.

        This constructor is used to instantiate an EST from the create
        method.

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
