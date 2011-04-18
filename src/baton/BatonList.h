#ifndef BATON_LIST_H
#define BATON_LIST_H

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

#include "Baton.h"
#include "WindowPair.h"
#include "EST.h"

#include <vector>
#include <queue>

// Forward declarations to keep compiler happy and fast


/** \file BatonList.h

    \brief Class to encapsulate the set of Batons associated with a
    given cDNA fragment.

    This class essentially contains the list of Baton objects
    associated with a given cDNA fragment.  The list is implicitly
    organized based on the <i>n</i>-mer heads/ends of the batons
*/

/** \typedef std::vector<Baton> NmerBatonList

    \brief A vector that contains the set of batons that have the same
    <i>n</i>-mer sequence as their baton heads.

    This typedef provides a convenient short cut to refer to a vector
    that is used to contain the set of batons with the same \e n -mer
    sequences at both ends.  This type is used to declare another
    aggregate vector that contains a list of NmerBatonList objects.
*/
typedef std::vector<Baton> NmerBatonList;

/** \typedef std::vector<NmerBatonList> WindowBatonList

	\brief A vector that contains all the <i>n</i>-mer batons in a
	given window.  The size of the window and the number of windows is
	determined when the BatonList class is created.  Given \e n, this
	list contains 4<sup>\e n</sup> entries (the base 4 arises from the
	use of the four nucleotides \c ATCG).  Each entry (from the
	4<sup>\e n</sup> entries) correspond to a unique encoding of \e n
	consecutive nucleotides.  Each entry in this list is a vector of
	Baton objects that have the same <i>n</i>-mer sequence as their
	baton heads/ends.

    \note Currently this typedef is not used.
*/
typedef std::vector<NmerBatonList> WindowBatonList;

/** \typedef std::vector<int> IntVector

    \brief A vector that contains a list of integers.

    This typedef provides a convenient shortcut to refer to a list of
    integers.  This vector is currently used to hold/return encoded
    values of <i>n</i>-mer sequences.
*/
typedef std::vector<int> IntVector;

/**
    <p>Class to encapsulate the set of Batons associated with a given
    cDNA fragment.  This class essentially contains the list of Baton
    objects associated with a given cDNA fragment.  The list is
    implicitly organized based on the <i>n</i>-mers constituting the
    ends of the batons.  The <i>n</i>-mers constituting a Baton are
    encoded into 2*\e n bit patterns using the ESTCodec class (base
    pairs \c A, \c T, \c C, and \c G (both upper and lower case) is
    encoded into 2-bits codes \c 00, \c 11, \c 10, and \c 01
    respectively).  For example, \c AAA is encoded to
    000000<sub>2</sub> or 0<sub>10</sub>.  Similarly the sequence \c
    TCA is encoded as 111000<sub>2</sub> or 56<sub>10</sub>.</p>

    <p>Note that the same encoding is used for normal and reverse
    complement sequences.  However, this class uses
    ESTCodec::NormalEncoder or ESTCodec::RevCompEncoder to
    appropriately encode each Baton's <i>n</i>-mer heads.  The choice
    of the encoder to be used is made when the BatonList is created
    via one of the constructors.  Refer to the API documentation on
    the constructors for further details.</p>
*/
class BatonList {
    // Insertion operator to make dumping Batons easier.
    friend std::ostream& operator<<(std::ostream&, const BatonList&);  
public:
    /** Constructor to create a BatonList for a given cDNA index.

        This constructor is used primarily to create the list of
        batons for a given cDNA fragment, given the index of the
        fragment.  This method assumes that the cDNA fragment has
        already been loaded into the EST class.
        
        \param[in] est A pointer to an immutable EST object containing
        the cDNA fragment for which the BatonList is to be built and
        maintained by this object.
        
        \param[in] nMerSize The length (in number of nucleotides) of
        the Baton ends.  This ultimately determines the number of
        unique baton ends to be considered and maintained by this
        object.
        
        \param[in] makeRC If this flag is \c true, then this method
        builds the batons for the reverse complementary representation
        of the specified cDNA fragment.  On the other hand, if this
        flag is \c false (the default setting) then this method builds
        batons for the regular cDNA sequence.

        \param[in] windowSize The size (in nucleotides) of the window
        into which the cDNA fragment must be subdivided in order
        identify regions of sufficient similarity.  Note that the
        window sizes set such that there is a 50% overlap between
        adjacent windows.
    */
    BatonList(const EST* est, const int nMerSize, const bool makeRC,
              const int windowSize = 100);

    /** Constructor to create a BatonList for a given cDNA fragment.

        This constructor is used primarily to create the list of
        batons for a dynamically generated cDNA fragment.  Such
        fragments are typically consensus sequences that are generated
        during assembly or other analysis.  These cDNA fragments
        cannot be found in the list of fragments (in the EST class)
        that were loaded from a given data file.
        
        \param[in] sequence The nucleotide sequence for which batons
        are to be generated and maintained by this object.  This
        parameter must point to a valid C string (terminated by '\\0')
        that is at least one nucleotide long.
        
        \param[in] nMerSize The length (in number of nucleotides) of
        the Baton ends.  This ultimately determines the number of
        unique baton ends to be considered and maintained by this
        object.
        
        \param[in] makeRC If this flag is \c true, then this method
        builds the batons for the reverse complementary representation
        of the specified cDNA fragment.  On the other hand, if this
        flag is \c false (the default setting) then this method builds
        batons for the regular cDNA sequence.

        \param[in] windowSize The size (in nucleotides) of the window
        into which the cDNA fragment must be subdivided in order
        identify regions of sufficient similarity.  Note that the
        window sizes set such that there is a 50% overlap between
        adjacent windows.
    */
    BatonList(const char *sequence, const int nMerSize,
              const bool makeRC = false, const int windowSize = 100);
    
    /** Determine if the baton list corresponds to a reverse
        complement nucleotide sequence.

        This method can be used to determine if this baton list
        contains batons for the the normal or the reverse
        complementary representation of the given cDNA fragment.  This
        value is initialized in the constructor and is never changed
        during the life time of an object.
        
        \return This method returns \c true if the baton correspond to
        a reverse complementary representation of a given nucleotide
        sequence.  For normal fragment, this method returns \c false.
     */
    inline bool isRCBatonList() const { return isRC; }

    /** Obtain the index (within the list of cDNA fragments) of the
        fragment within the EST class.

        \return The index of the cDNA fragment (to be used in
        EST::getEST() method) for which this object conaints batons.
        If an index is not available (because the sequence has been
        dynamically generated) then this method returns -1.
    */
    inline int getIndex() const { return est->getID(); }

	/** Determine the number of windows into which the batons has been
		subdivided.

		This method can be used to determine the number of windows
		into which the cDNA fragment associated with this list has
		been sub-divided.  The number of windows is set based on the
		average window size that was specified when this object was
		instantiated.

		\return The number of windows into which the batons have been
		logically sub-divided.

        \see getWindowSize
        \see getWindow
	*/
	inline int getWindowCount() const { return windowCount; }

	/** Obtain the size (in number of nt) into which the cDNA fragment
		has been logically sub-divided.

		The window size value is set when the baton list is created in
		the buildBatons() method.  The window size is set such that
		the number of nucleotides are evently distributed between the
		various windows.  The suggested/typical window size is 100 nt.
		
		\return The window size (in number of nucleotides) of each
		window.

        \see getWindowCount
        \see getWindow
	*/
	inline int getWindowSize() const { return windowSize; }

	/** Determine logical window number for a given index position.

        The set of batons contained in this baton list are logically
		subdivided into a set of overlapping windows for baton-based
		analysis.  This is a convenience method that can be used to
		determine the logical (zero-based) window number for a given
		nucleotide (nt) index position.  The window size value is set
		when the baton list is created.  The window size is set such
		that the number of nucleotides are evently distributed between
		the various windows.

        \param[in] ntIndex The zero-based index position of a
        nucleotide (or a baton starting index position) to be
        translated to its corresponding logical window number.

        \return The logical (zero-based) window number into which the
        given index position falls.
        
        \see getWindowCount
        \see getWindowSize
        \see Baton::getStartIndex
    */
	inline int getWindow(const int ntIndex) const {
        return ntIndex / windowSize;
    }
	
    /** Obtain the nucleotide sequence for which this class is
        maintaining batons.

        This method correctly returns the nucleotide sequence (for
        which batons are being maintained by this object) by handling
        the following two cases:

        <ol>

        <li>In cases when the index is known, this method obtains the
        fragment from the shared list maintained by the EST class.</li>

        <li>If this class contains batons for a dynamically generated
        fragment (that is, the getIndex() method returns -1), the
        sequence is obtained from the local nucleotide sequence
        maintained in this class.

        </ol>

        \return The nucleotide sequence for which this object is
        maintaing the batons.
    */
    inline const char* getSequence() const {
        return (est == NULL) ? sequence.c_str() : est->getSequence();
    }

    /** Obtain the length of the nucleotide sequence for which this
        class is maintaining batons.

        This method correctly returns the nucleotide sequence (for
        which batons are being maintained by this object) by handling
        the following two cases:

        <ol>

        <li>In cases when the index is known, this method obtains the
        fragment from the shared list maintained by the EST class.</li>

        <li>If this class contains batons for a dynamically generated
        fragment (that is, the getIndex() method returns -1), the
        sequence is obtained from the local nucleotide sequence
        maintained in this class.

        </ol>

        \return The nucleotide sequence for which this object is
        maintaing the batons.
    */
    inline int getSequenceLength() const {
        return (est == NULL) ? sequence.size() : est->getSequenceLength();
    }

    /** Obtain the list of batons that all have the same <i>n</i>-mer
        heads/ends.

        This is a convenience method that can be used to obtain an
        immutable reference to the list of batons maintained
        internally by this class.  This method does not perform any
        special checks or processing but simply returns the
        pre-computed vector of batons.

        \param[in] nMerCode The encoded value of the \e n -mer baton
        ends corresponding to which the batons are to be returned.
        The encoding of values is assumed to follow the encodings
        provided by the ESTCodec class. No special checks are done on
        this parameter.  Consequently, it is the callers
        responsibility to ensure that the parameter is within valid
        range.

        \return A vector containing the set of batons that all have
        the same <i>n</i>-mer ends as specified by the encoded
        nMerCode parameter value.

        \see ESTCodec::NormalEncoder
        \see ESTCodec::RevCompEncoder
        \see BatonAnalyzer::align
    */
    inline const NmerBatonList& operator[](const int nMerCode) const {
        return batons[nMerCode];
    }
    
    /** Search/analysis method to obtain list of window-pairs that
        contain more than a given threshold of similar batons.

        This method essentially calls the other, overloaded method to
        obtain the necessary information.  It then walks the sets of
        windows and adds those pairs of windows that exceed the given
        threshold to the supplied vector.  Note that the information
        about the identical baton frequencies in all possible windows
        (including those that are below the specified threshold) is
        not preserved by this method.

        \param[in] other The other baton list with which the search
        for identical batons in windows must be performed.  Note that
        \c this and \c other baton list may have slightly different
        window sizes depending on the length of their cDNA fragments.

        \param[out] pairList This vector is populated by this method
        and is the primary return value from this method. This vector
        contains the pairs of windows whose identical baton count
        exceeds the specified threshold. The information is stored in
        the form of WindowPair objects.

        \param[in] threshold The number of identical batons that must
        be present in a given window pair in order for the pair to be
        added to the pairList.
    */
    void getWindowPairs(const BatonList& other,
                        std::vector<WindowPair>& pairList,
                        const int threshold = 7) const;

    /** Search/analysis method to obtain list of window-pairs that
        contain more than a given threshold of similar batons.

        This method essentially calls the other, overloaded method to
        obtain the necessary information.  It then walks the sets of
        windows and adds those pairs of windows that exceed the given
        threshold to the supplied priority queue.  The priority queue
        uses operator<() on WindowPair objects to sort them such that
        window-pairs with highest number identical batons are at the
        top of the priority queue.  Note that the information about
        the identical baton frequencies in all possible windows
        (including those that are below the specified threshold) is
        not preserved by this method.

        \param[in] other The other baton list with which the search
        for identical batons in windows must be performed.  Note that
        \c this and \c other baton list may have slightly different
        window sizes depending on the length of their cDNA fragments.

        \param[out] pairList This priority queue is populated by this
        method and is the primary return value from this method. This
        queue contains the pairs of windows whose identical baton
        count exceeds the specified threshold. The information is
        stored in the form of WindowPair objects.  The entries are
        sorted with the window-pair with highest number of identical
        batons at the top of the queue.

        \param[in] threshold The number of identical batons that must
        be present in a given window pair in order for the pair to be
        added to the pairList.
    */
    void getWindowPairs(const BatonList& other,
                        std::priority_queue<WindowPair>& pairList,
                        const int threshold = 7) const;
    
    /** The search/analysis method to obtain list of window-pairs that
        contain more than a given threshold of similar batons.

        This method performs the core task of comparing similar (with
        the same n-mer heads) batons from \c this baton list and \c
        other baton list to determine logical window-pairs whose baton
        frequency exceeds the given threshold.  Such window pairs are
        the candidate windows for alignment.  The window-pairs are
        integer values where the first value represents the logical
        window in \c this baton list and the second number indicates
        the logical window in the other baton list.  This method
        operates as follows:

        <ol>

        <li>The window sizes of the two baton lists (\c this and \c
        other) are used to compute the number of windows and
        identicalBatonCount is initialized to create a 2-D array of
        zeros.</li>

        <li>Next, the baton lists for each n-mer are searched to find
        batons of the same length.  The search proceeds assuming that
        the lists are sorted based on the baton lengths (sorting
        avoids a n<sup>2</sup> approach as this method is a frequently
        used one). For each identical baton the tallyBatons() method
        is used to update the identicalBatonCount 2-D array and add
        window-pairs to the pairList if the baton count exceeds the
        threshold.

        </ol>

        \param[in] other The other baton list with which the search
        for identical batons in windows must be performed.  Note that
        \c this and \c other baton list may have slightly different
        window sizes depending on the length of their cDNA fragments.

        \param[out] identicalBatonCount This is the primary return
        value from this method.  This vector-of-vectors is used to
        conveniently implement a 2-D array of integers.  The first
        dimension is the logical window number in \c this batons while
        the second dimension is the logical window number in \c other
        batons.  The value \c identicalBatonCount[win1][win2] contains
        the number of identical batons in the window pair \c <win1,
        win2>. The overhead of using the vector-of-vectors is
        sufficiently minimal when compared to that of a vanilla 2-D
        array.  However, the vector-of-vectors approach has its own
        adavntages (no need to worry about dynamic memory allocation).
        
        \param[in] threshold The number of identical batons that must
        be present in a given window pair in order for the pair to be
        added to the pairList.
    */
    void getWindowPairs(const BatonList& other,
						std::vector< IntVector >& identicalBatonCount,
                        const int threshold = 7) const;	

	/** Utility method that can be used to get codedNmers.

		This is a convenience utility method that can be used to get
		encoded \e n -mers used by this baton list.  This method
		provides encoded sequence for either the normal or reverse
		complement version of the cDNA fragment associated with this
		baton list.

		\param[out] codedNMers The list of coded \e n -mers are added
		to this vector as they are computed.
		
		\note This method performs the encoding each time it is
		called. So use it sparingly.  Thee encoded <i>n</i>-mer list
		is not cached because its usage lifetime is rather short -- it
		is needed only during alignment, once it has been determined
		that two cDNA fragments must be aligned.  It is not needed for
		general searching
	*/
	void getCodedNMers(IntVector& codedNMers) const;

    /** Obtain the number of <i>n</i>-mers that constitute the baton
        heads.

        This method can be used to determine the number of
        <i>n</i>-mers that constitute the head/ends of each baton in
        this list.  This method essentially returns the value that was
        specified when this class was instantiated.

        \return The number of nucleotides that determine the size of
        the baton heads.
    */
    inline int getNmerSize() const { return nMerSize; }

    /** Convenience method to determine if a given baton lies within
        the bounds of a given window.

        This method can be used to determine if a baton from this
        baton list lies within the bounds of a valid window.

        \param[in] baton The baton to be checked to see if it lies
        within the bounds of the given window.

        \param[in] window The logical index number of the window
        within this baton list whose containing bounds need to be
        checked.

        \return This method returns \c true if the \c baton fits
        within the logical bounds of the specified \c window.
        Otherwise this method returns \c false.
    */
    inline bool isBatonInWindow(const Baton& baton, const int window) const {
        return ((baton.getStartIndex() >= (window * windowSize / 2)) &&
                (baton.getStartIndex() <= ((window + 1) * windowSize / 2)));
    }
    
protected:
    /** Helper method to compute the list of encoded <i>n</i>-mers
        associated with the given nucleotide sequence.

        This method is indirectly invoked from the non-templatized
        getCodedNMers method.  This method performs the core task of
        building a series of <i>n</i>-mers and adding them to the
        specified vector.  The list is used for creating batons and
        later on for initial alignment.  The encoding of <i>n</i>-mers
        is performed using the supplied encoder.
        
		\tparam Encoder The type of encoder to be used.  This value
		can be ESTCodec::NormalEncoder or
		ESTCodec::RevCompEncoder.
        
        \param[in] seq The sequence of nucleotides to be processed to
        construct the batons.

        \param[out] codedNMers The vector to which the encoded
        <i>n</i>-mers are to be added.
        
		\param[in] encoder The encoder object to be used to generate
		encoding for the <i>n</i>-mers generated by this table.
    */
    template <typename Encoder>
    void getCodedNMers(const char* seq, IntVector& codedNMers,
                       Encoder encoder) const {
        // Compute the nMerCode for the first n-1 'mers using the
        // encoder. The encoder may do normal or reverse-complement
        // encoding for us.
        int nMerCode = 0;
        for(int i = 0; ((*seq != 0) && (i <  nMerSize - 1)); i++, seq++) {
            nMerCode = encoder(nMerCode, *seq);
        }
        // Reserve sufficient space in codedNMers array to minimize
        // memory reallocation.
        codedNMers.reserve(getSequenceLength());
        // Now compute the nMerCode for each n-mer and add to the vector
        for(int idx = 0; (*seq != 0); idx++, seq++) {
            nMerCode = encoder(nMerCode, *seq);
            codedNMers.push_back(nMerCode);
        }
    }
    
    /** Helper method to build batons.

        This is a helper method that is invoked from the various
        constructors to perform the task of actually building the
        batons.  This method performs the following tasks:

        <ol>

        <li>It creates a suitable encoder depending on whether a
        normal or a reverse complement nucleotide sequence is to be
        used for creating Batons.</li>

        <li>It then uses the templatized buildBatons() method to
        actually populate the \c batons vector with all the batons for
        the given nucleotide sequence.</li>

        <li>Finally, it sorts the batons associated with <i>n</i>-mer
        based on the length of the batons to optimize the frequently
        used baton comparison operations.</li>

        </ol>
        
        \note The constant instance variables used by this method are
        assumed to be initialized by the time this method is invoked
        by the constructor(s).  Since this method uses the values in
        the instance variables, it does not require any parameters.

        \param[in] winSize The suggested (average) size (in number of
        nucleotides) of a window to be used when searching for batons.
    */
    void buildBatons(const int winSize);

    /** Helper method to track identical batons and update various
        data structures (passed in as parameters).

        This method is a helper method and operates only using the
        parameters passed-in to the method.  This method is invoked
        from the main getWindowPairs() method each time an identical
        baton (batons have same n-mer ends, and are of the same
        length) is encountered.  This method uses a series of if-tests
        to suitably update the identicalBatonCount 2-D vector based on
        the value of baton1Win and baton2Win.

        \param[in] baton1Win The logical window in the first cDNA
        fragment (with which \c this baton list is associated).  This
        value also identifies the first dimension in the
        identicalBatonCount 2-D vector.

		\param[in] numWin1 The total number of windows in the first
		baton.
		
        \param[in] baton2Win The logical window in the second cDNA
        fragment.  This value also identifies the second dimension in
        the identicalBatonCount 2-D vector.

		\param[in] numWin2 The total number of windows in the second
		baton.
		
        \param[out] identicalBatonCount The 2-D vector-of-vectors to
        be updated with the number of matching baton counts.
    */
    void tallyBatons(const int baton1Win, const int numWin1,
                     const int baton2Win, const int numWin2,
                     std::vector< IntVector > &identicalBatonCount) const;
        
private:
    /** The pointer to an immutable EST object containing the cDNA
        fragment for which this class is maintaining the list of
        batons.

        This instance variable provides a convenient reference to the
        cDNA fragment for which this class is maintaining the list of
        Batons.  Note that this value can be NULL, if a valid object
        is not available (and this is the case when dealing with
        generated, partial consensus sequences).  This value is
        initialized in the constructor and is never changed during the
        life time of an object.
    */
    const EST* const est;

    /** The number of nucleotides constituting the baton heads.

        This instance variable contains the number of of nucleotides
        constituting the baton heads (or ends).  For an <i>n</i>-mer
        baton head, this instance variable contains <i>n</i>.  This
        value is initialized in the constructor and is never changed
        during the life time of an object.
    */
    const int nMerSize;

    /** Flag to indicate if batons correspond to regular or reverse
        complement sequence.

        This flag is used to determine if this baton list contains
        batons for the the normal or the reverse complementary
        representation of the given cDNA fragment.  This value is
        initialized in the constructor and is never changed during the
        life time of an object.
    */
    const bool isRC;
    
    /** The nucleotide sequence for which this baton list is
        maintaining batons.

        This instance variable provides a convenient reference to the
        nucleotide sequence for which this class is maintaining the
        list of batons.

        \note The nucleotide sequence is only set if the estIndex is
        -1 and the object is maintaining batons for a dynamically
        generatedc consensus sequence.  This is done to minimize the
        memory footprint of BatonList objects.
    */
    const std::string sequence;

    /** The complete list of batons associated with a given cDNA
        fragment.

        This vector contains the complete set of batons associated
        with the given cDNA fragment.  Each entry in this list is a
        vector which is indexed using the encoded values for baton
        heads.  For example, given <i>3</i>-mer batons, the baton head
        sequence \c TCA is encoded as 111000<sub>2</sub> or
        56<sub>10</sub> and batons[56] contains a vector of batons
        (with different lengths and starting index positions) that all
        have \c TCA as their baton heads. This vector is populated by
        the templatized buildBatons() method in this class.
    */
    std::vector<NmerBatonList> batons;

    /** The number of base-pairs that constitute a window for the set
        of batons maintained by this class.

        This value is set when the baton list is created in the
        buildBatons() method.  The window size is set such that the
        number of nucleotides are evently distributed between the
        various windows.  The suggested/typical window size is 100 nt.
    */
    int windowSize;

    /** The number of windows into whcih the batons have been
        subdivided.

        The number of windows is computed in the buildBatons() method
        (that is called from the constructors). The window count is
        set based on the average window size that was specified when
        this object was instantiated.
    */
    int windowCount;
    
    /** Instance variable to store the number of bits to be shifted to
        create <i>n</i>-mer encodings.

        <p>This instance variable is set to the value of 2. This is
        because each base pair \c ATCG is encoded using a 2-bit
        encoding.</p>

        <p>This instance variable is actually passed on to the
        ESTCodec::NormalEncoder or ESTCodec::RevCompEncoder when
        computing hash values.  Since this is value is passed in a
        template parameter, it is defined to be static (to ensure that
        it has external linkage as per the ISO/ANSI standards
        requirement).</p>
    */
    static int bitShift;
    
    /** The bitmask to be used when building <i>n</i>-mer encodings.
            
        This bitmask is used to drop higher order bits when building
        encodings for the <i>n</i>-mers constituting baton heads.
        This instance variable is actually passed on to the
        ESTCodec::NormalEncoder or ESTCodec::RevCompEncoder when
        computing encodings.  Since this is value is passed in a
        template parameter, it is defined to be static (to ensure that
        it has external linkage as per the ISO/ANSI standards
        requirement).
    */
    static int bitMask;
};

/** \fn std::ostream& operator<<(std::ostream&, const BatonList&)

    Insertion operator to stream BatonList information to a given
    output stream.  This method provides a convenient mechanism to
    dump the complete BatonList information for debugging purposes.
    This method displays baton information for each <i>n</i>-mer on
    its own separate line.

    \param[out] os The output stream to which the baton list
    information is to be written.

    \param[in] batonList The object whose information is to be
    displayed.

    \return As the generalized API contract for insertion operators,
    this method returns os (first parameter)
*/
extern std::ostream& operator<<(std::ostream& os, const BatonList& batonList);

#endif
