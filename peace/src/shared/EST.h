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

/** \typedef std::vector<int> QualityVector

    \brief Shortcut typedef for std::vector<int>

    This type definition is a shortcut to define the data type that is
    used to manage the quality data associated with the various bases
    in a given cDNA read.  The type definition is primarily used to
    ease changing the native data type used to store the information
    when the need arises.
*/
typedef std::vector<unsigned char> QualityVector;

/** A single EST.

    This class is used to represent a single EST.  An EST object
    instance consists of the following information:

    <ul>

    <li>\b id: A unique identifier (usually a number) for this EST.
    Typically entries are assigned consecutive IDs and some parts of
    PEACE count on fragments having consecutive integer ID
    values.</li>
    
    <li>\b info: The name and other information associated with the
    EST.  This information is typically the first header line read
    from a FASTA file.</li>

    <li>\b sequence: The actual sequence of base pairs associated with
    this EST.<li>
	
    <li>\b quality: The quality vector associated with each nucleotide
    entry in the fragment.  The quality is typically stored as a
    log-normal integer value in the range 0 to 255.</li>

    </ul>

	\note Not all of the information is always stored in this
	object. Some of the information is loaded on-demand to reduce
	memory footprint.
*/
class EST {
    friend class InputFile;
public:
    /** EST Constructor (without quality information).

        This constructor is used to instantiate an EST object
	
        \param[in] id The unqiue ID value to be set for this EST.
        
        \param[in] info The name and other information associated with
        the EST.  This information is typically the first header line
        read from a FASTA file.  This information can be an empty
        string.
        
        \param[in] sequence The actual sequence of base pairs
        associated with this EST.  The sequence information that must
        be used to create this EST.  The sequence information can be
        an empty string (\c ""). Any nucleotide masking and base
        normalization must be performed pior to invoking this
        method. See EST::normalizeBases() method for details.

		\param[in] maskBases If this flag is \c true, then all the
		lower-case characters in the input sequence are first
		converted to \c 'N' base characters.  If this flag is \c
		false, then all lower-case nucleotides are converted to
		upper-case equivalents.

		\param[in] randomizeNbases If this flag is \c true, then all
		entries with 'N' are randomly changed to one of \c ATCG.		
    */
    EST(const int id, const std::string& info,
        const std::string& sequence = "");

    /** EST Constructor (with per-nucleotide quality information).

        This constructor is used to instantiate an EST object with
        additional quality information.
     
        \param[in] id The unqiue ID value to be set for this EST.
        
        \param[in] info The name and other information associated with
        the EST.  This information is typically the first header line
        read from a FASTA file.
        
        \param[in] sequence The actual sequence of base pairs
        associated with this EST.  The sequence information that must
        be used to create this EST.  The sequence information can be
        an empty string.  Any nucleotide masking and base
        normalization must be performed prior to calling the
        constructor via EST::normalizeBases() method.
        
        \param[in] quality The quality vector to be set for this
        entry.  The number of entries in this vector must match the
        length of the nucleotide sequence parameter.

		\param[in] maskBases If this flag is \c true, then all the
		lower-case characters in the input sequence are first
		converted to \c 'N' base characters.  If this flag is \c
		false, then all lower-case nucleotides are converted to
		upper-case equivalents.

		\param[in] randomizeNbases If this flag is \c true, then all
		entries with 'N' are randomly changed to one of \c ATCG.		
    */    
    EST(const int id, const std::string& info,
        const std::string& sequence, const QualityVector& quality);

    /** Copy constructor.

        This method can be used to make a copy of an EST. This method
        copies all the currently available information for the given
        EST.  However, if some information has been unpopulated in the
        source, that information will also be absent in the cloned
        object returned by this method.

        \param[in] src The source EST from where the information is to
        be copied into this EST.
    */    
    EST(const EST& src);
    
    /** Dump this EST information in FASTA format.

        This method can be used to dump the information associated
        with the EST in FASTA format to a given output stream.

        \param[in] os The output stream to which the EST's information
        must be written in FASTA format.
    */
    void dumpEST(std::ostream& os);

    /** Obtain the ID of this EST.
        
        \return The ID of the EST that was set when this EST was
        created.
    */
    inline int getID() const { return id; }
    
    /** Obtain the information associated with this EST.
        
        The name and other information associated with the EST.  This
        information is typically the first header line read from a
        FASTA file.  This information can be an empty string if the
        EST is only partially loaded.

        \return Any information available for this EST. The returned
        information must not be modified or deleted.
    */
    inline const char* getInfo() const { return info.c_str(); }

    /** Obtain the actual sequence of base pairs for this EST.

        Note that sequence information for an EST can be empty if it
        is only partially loaded from a file.  Entries are parially
        loaded to reduce memory foot print when processing large data
        sets.
        
        \return The actual sequence of base paris for this EST if
        available. Otherwise this method returns an empty string.
    */
    inline const char* getSequence() const { return sequence.c_str(); }

	/** Obtain the reverse-complementary sequence for the base pairs
		for this EST.

        This is a convenience method that can be used to build the
        reverse complementary representation of a given nucleotide
        sequence.  In other words, given the sequence \c AATCGG this
        method returns \c CCGATT.

        \note This method computes the RC-sequence each time it is
        called.  The time-complexity is O(n).
        
        \return The reverse complementary nucleotide sequence for this
        EST.  The length of the returned sequence will be exactly the
        same as the original nucleotide sequence.
	*/
	std::string getRCSequence() const;
	
	/** Convenience method to detect if this cDNA fragment has quality
		values associated with each nucleotide.

		This method can be used to detect if this cDNA fragment has
		quality values associated with each nucleotide. 

		\return This method returns \c true if the quality values are
		present for this cDNA fragment. If quality values are not
		available, then this method returns \c false.
	*/
	inline bool hasQuality() const { return !quality.empty(); }
	
	/** Obtain the quality vector associated with this EST entry.

		This method can be used to obtain the quality values
		associated with each nucleotide in the cDNA fragment. If
		quality values are available, each entry in the vector
		provides the quality value for the corresponding nucleotide
		entry. That is, <code>getQuality[i]</code> (0 &le; i &lt
		<code>getSequenceLength()</code>) provides the quality value
		for the nucleotide at <code>getSequence[i]</code>.  The
		quality values for each nucleotide is the phred-scaled base
		error probability value which is calculated as
		-10*log<sub>10</sub> Pr{base is wrong}.

		\return A vector containing the phred-scaled base error
		probability for each nucleotide in this cDNA fragment. If
		quality values are not available (or are not currently
		populated) then this method returns an empty vector.
	*/
	inline const QualityVector& getQuality() const { return quality; }
	
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

    /** Obtain the length of the sequence.

        This method is the preferred approach for obtaining the length
        of the EST sequence.  This method merely returns the
        pre-computed (computed when the entry was created) length,
        thereby saving CPU cycles.

        \return The pre-computed length for this EST.
    */
    inline int getSequenceLength() const { return sequenceLen; }

    /** Determine if this EST has all the pertinent information.

        EST objects can be unpopulated (or cleared) by calling the
        unpopulate() API method.  This method can be used to determine
        if an EST object has been unpopulated or if it contains all
        the necessary information.

        \return This method returns true if the EST object is fully
        populated. If it returns false, that means this entry has been
        unpopulated.
    */
    inline bool isPopulated() const { return populated; }

    /** Helper method to normalize a given nucleotide sequence based
        on options.

        This method is used to normalize fragments typically read from
		a InputFile file. This method normalizes the sequences such
		that the resulting sequence is over the set {'A', 'T', 'C',
		'G', 'N'} in the following manner:
		
		<ul>
		
		<li>If the maskBases flag is true, then all lowercase
		nucleotides are converted to 'N'. Otherwise they are converted
		to uppercase equivalents.</li>
		
		<li>All nucleotides that are not in \c "ATCG" are converted to
		'N'.</li>
		
		</ul>
		
		\param[in] srcSequence The source sequence of nucleotides to
		be normalized by this method.

		\param[in] maskBases If this flag is \c true, then all the
		lower-case characters in the input sequence are first
		converted to \c 'N' base characters.  If this flag is \c
		false, then all lower-case nucleotides are converted to
		upper-case equivalents.

		\param[in] randomizeNbases If this flag is \c true, then all
		entries with 'N' are randomly changed to one of \c ATCG.
		
        \param[in] randSeed The seed to be used for randomizing
		bases. Setting the same seed for a given sequence ensures that
		the characters are consistently pseudo-random.  Having a
		consistent pseudo-random characters enables consistent
		processing of cDNA fragments.

		\return This method returns a normalized version of the given
		source sequence.
    */
    static std::string normalizeBases(const std::string& srcSequence,
									  const bool maskBases,
									  const bool randomizeNbases,
									  const int randSeed = 0);
	
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
    std::string info;

    /** The actual sequence of base pairs associated with this EST.
        This information is typically read from a FASTA file.  The
        information may be dynamically loaded on demand to reduce
        memory footprint when processing large data sets.
    */
    std::string sequence;

    /** The length of the sequence in number of characters.

        This instance variable essentially maintains the string length
        of the sequence.  This value is initialized in the constructor
        method and is never changed during the life time of the
        object.
    */
    const int sequenceLen;
    
    /** The per-nucleotide quality information.

        This vector contains the per-nucleotide quality
        information. This information is typically read from a FASTAQ
        or SFF file.  The information may be dynamically loaded on
        demand to reduce memory footprint when processing large data
        sets.
    */
    QualityVector quality;

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
   
    /** Place holder for some other custom data.
        
        This pointer acts as a convenient place holder for other
        classes to associate some uninterpreted user data (or data
        structure).  This member is initialized to NULL in the
        constructor.  Note that this pointer is managed using an
        auto_ptr that automatically deletes the data when the auto_ptr
        loses ownership of the data object.
    */
    std::auto_ptr<ESTCustomData> customData;

    /** Flag to indicate if the data for this EST has been unpopulated.

        This flag is used to track if the information in this EST has
        been unpoulated (or cleared out) via a call to the
        unpopulate() method.  If so this flag is set to false. By
        default this flag is initialized to true.
    */
    bool populated;

    /** Repopulate the information in this EST entry.

        This method is currently invoked only from the InputFile class
        to repopulate an EST entry. This method resets the #populated
        flag back to true.

        \param[in] info The name and other information associated with
        the EST.  This information is typically the first header line
        read from a FASTA file.
        
        \param[in] sequence The actual sequence of base pairs
        associated with this EST.  The sequence information that must
        be used to create this EST.  The sequence information can be
        an empty string.  Any nucleotide masking and base
        normalization must be performed pior to invoking this
        method. See InputFile::normalizeBases() method for details.
        
        \param[in] quality The quality vector to be set for this
        entry.  The number of entries in this vector must match the
        length of the nucleotide sequence parameter.
    */
    void repopulate(const std::string& info, const std::string& sequence,
                    const QualityVector& quality);
                    
private:    
    /** The default constructor.

        The default constructor has been made private to ensure that
        EST's are never directly created.  Instead, a valid EST must
        be created using other constructors.
    */
    EST();

    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.

        \return Reference to this.
    */
    EST& operator=(const EST& src);
};

#endif
