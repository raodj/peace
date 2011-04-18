#ifndef PARAMETER_SET_MANAGER_H
#define PARAMETER_SET_MANAGER_H

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

#include "ESTListListener.h"
#include "ParameterSet.h"
#include "Utilities.h"
#include <vector>

// Forward declarations
class ESTList;

/** Class that manages a list of parameter sets.
    
    <p>This class represents a shared list of parameter sets that
    analyzers and heuristics can access to obtain appropriate
    parameters for analysis of sequences with varying characteristics.
    The distinguish feature of sequences (that determines the
    parameter set to be used) is the length of the sequences.  The
    class initializes the appropriate parameter sets and provides the
    necessary parameter values when requested. The actual parameter
    values are encapsulated into ParameterSet objects. </p>

    <p> This class is (and must be used) in the following manner:

    <ol>

    <li>First the ParameterSetManager constructor must be invoked with
    the approriate parameters. The constructor sets up parameters for
    three different classes (based on length) of cDNA fragments.<li>

    <li>Next, the fragments to be analyzed and clustered are loaded
    from an appropriate data file(s) into the shared
    RuntimeContext::estList.</li>

    <li>Prior to each sub-system commencing its run, the Component
    that is interested in using adaptive parameters must invoke
    initialize() method.  This method ensures that the internal data
    structures are updated to reflect the number of ESTs in the list.
    This method pre-computes look-up tables to provide dynamic
    adaptation parameters quickly.</li>

    <li>The ParameterSetManager::getParameterSet() method is
    repeatedly invoked to obtain appropriate parameter to be used.<li>
    
    </ol>
    
    </p>
*/
class ParameterSetManager : ESTListListener {
    friend class HeuristicChain;
public:
    /** Gets the appropriate parameter set for two sequences given
        their indexes.

        \note This method must be called only after the
        ParameterSetManager has been initialized via call to the
        initialize() method.

        \param[in] seq1Index The index of the first sequence.  No
        special checks are made to ensure that the index is valid (for
        performance).

        \param[in] seq2Index The index of the second sequence.  No
        special checks are made to ensure that the index is valid (for
        performance).
        
        \return The index of the parameter set that should be used to
        compare the two sequences; -1 if the sequence lengths are too
        different and should not be compared.
    */
    inline int getParameterSet(const int seq1Index,
                               const int seq2Index) const {
        ASSERT((seq1Index >= 0) && (seq1Index < (int) lengthToParamSetMap.size()));
        ASSERT((seq2Index >= 0) && (seq2Index < (int) lengthToParamSetMap.size()));
        return lookupTable[lengthToParamSetMap[seq1Index]][lengthToParamSetMap[seq2Index]];
    }

    /** Obtain the parameter set at a given index.

        This method must be used to obtain the actual parameter set
        that is associated with a given index in the parameter set
        table.

        \param[in] index The index in the parameter set for which the
        actual parameter set is to be returned.  Note that no checks
        are made to validate the index.  Consequently, calling this
        method with an invalid index will cause issues.
    */
    inline const ParameterSet* getParameterSet(const int index) const {
        return parameterSets[index];
    }

    /** Setup the parameter set manager for handling a hybrid
        (sanger/454/illumina) data set.
	
        This method is typically invoked from an clustering or
        ESTAnalysis class that is designed to handle adaptive
        operations (such as: AdaptiveTwoPassD2).  This method is used
        to setup parameters for analyzing a hybrid data
        sets. Specifically, it initializes the parameter set list to
        default values.  The default is for 3 parameter sets, one for
        each of the following "classes" of sequences:

        <table>
        <tr><th>Type</th><th>Min Len</th><th>Max Len</th><th>Notes</th></tr>

        <tr><td>Short</td><td>0</td><td>150</td><td>In practice range
        is 50 to 150, as the LengthFilter will filter out sequences of
        fewer than 50 bases long</td></tr>
        
        <tr><td>Medium</td><td>150</td><td>400</td><td></td></tr>

        <tr><td>Long</td><td>400</td><td>no limit</td><td>The
        practical upper limit is really 65K</td></tr>
        </table>

        \note The list of category entries must be organized in
        ascending order of fragment lengths.  That is, the shorter
        fragment entries must occur first in the parameters to this
        method.  The fragment length values (passed-in as parameters)
        may not overlap.
        
        In the future this class should incorporate command-line
        arguments and customizability for the parameter sets.
        However, for now the defaults work just fine for a broad range
        of data sets.        
    */
    void setupParameters(int t1 = 20, int u1 = 4, int ws1 = 4,
                         int t2 = 35, int u2 = 6, int ws2 = 8,
                         int t3 = 50, int u3 = 8, int ws3 = 8);

    /** Method to setup adaptive parameters to \i fast(er) settings.
        
        <p>This version of the overloaded method must be used to setup
        the parameters to a single entry.  The entry is chosen to be
        optimal for longer reads.  Conseuqently, this option should be
        used when it is known that the data is from sanger sequences
        that typically have longer reads.  However, the single entry
        enables the clustering to run faster.</p>
        
        <p><b>Note:</b> This method is primarily meant to be used to
        configure parameters for longer reads.  It does cause the
        clustering to run faster, but should not be used with short
        reads.</p>
        
        \param[in] dummy This is just a dummy parameter to
        disambiguate the two overloaded version.  This parameter is
        not used by this method.
    */
    void setupParameters(const bool dummy);
	
    /** Initialize parameter sets by pre-computing a rapid look up table.

        <p>This method must be invoked after the parameters have been
        setup (via call to ParameterSetManager::setupParameters
        method) and the fragments to be analyzed have been loaded.</p>
        
        <p>This method first populates the legthToParamSetMap for each
        fragment to be analyzed.  Next it creates and populates the
        lookupTable using the various entries in the parameterSets
        vector.</p>

        <p><b>Note:</b> This method will be called multiple times in the
        following manner: <ol>
        
        <li> First when the FilteringSubSystem starts running (after
        dummy ESTs have been added).</li>
        
        <li> Another time when the ClusteringSubSystem starts
        running.</li>
        
        </ol></p>
    */
    void initialize();
    
    /** Gets the maximum possible frame size over all parameter sets.

		\note Prior to calling this method the setupParameters()
		method must have been invoked followed by the initialize()
		method.  If the aforementioned two methods are not invoked,
		then this method will result in a runtime error as the
		ParameterSetManager::parameterSets vector will have zero
		entries (using access to -1<sub>th</sub> entry in the vector).
		
        \return the maximum possible frame size, an integer.
    */
    inline int getMaxFrameSize() const {
        return parameterSets[parameterSets.size() - 1]->frameSize;
    } 
    
    /** Adds a parameter set to the parameter set manager.

        \param[in] p The parameter set to be added to this
        manager. The manager takes ownership of the parameters object
        passed in and arranges to have it deleted in the destructor.
    */
    void addParameterSet(ParameterSet* p);

    /** Handle notifications regarding entries added to the ESTList.

        This method is invoked (by one of the components in a
        SubSystem) when a set of EST entries have been added.  This
        method appropriately adds entries to the lengthToParamSetMap
        to handle look-up of the newly added entries.  This method is
        called only after initialization of the ParameterSetManager is
        complete.

        \param[in] startIndex The zero-based starting index from where
        new EST entries were added.

        \param[in] endIndex The zero-based ending index (into the
        shared list of cDNA fragments pointed by
        RutimeContext::estList) where the last entry was added.

        \note It is guaranteed that startIndex <= endIndex.
    */
    void entriesAdded   (int startIndex, int endIndex);

    /** Handle notifications regarding entries removed from the
        ESTList.

        This method is invoked (by one of the components in a
        SubSystem) when a set of EST entries have been removed.  This
        method appropriately updates (or removes) entries from the
        lengthToParamSetMap vector corresponding to the entries that
        have been removed.

        \param[in] startIndex The zero-based starting index from where
        new EST entries were removed. Note that this index value would
        not be really valid as entries have already been removed.

        \param[in] endIndex The zero-based ending index (into the
        shared list of cDNA fragments pointed by
        RutimeContext::estList) where the last entry was removed.

        \note It is guaranteed that startIndex <= endIndex.

        @see ESTListListener
    */    
    void entriesRemoved (int startIndex, int endIndex);
    
    /** The destructor.
        
        The destructor frees up all the parameter sets added to this
        parameter set manager.
    */
    virtual ~ParameterSetManager();

    
protected:
    /** Helper method to return the index within the parameter set for
        a given length.

        This method searches the list of entries in the parameter set
        index to determine the appropriate parameter set entry to be
        used.  The search is based on the given length.
        
        \param[in] estLength The length of the fragment for which a
        suitable parameter set is to be determined.

        \return The index in the parameterSets vector that contains
        the parameter set entry that accommodates the given fragment
        length.
    */
    int getParameterSetIndex(const int estLength) const;

    /** Helper method to return the preferred parameter set index
        given two parameter sets.

        This helper method is called from the initialize() method to
        compute the preferred index (into the parameterSets vector)
        based on two parameter set ranges.  This method operates as
        follows:

        <ol>

        <li>It obtains the ParameterSet objects at \c index1 and \c
        index2.</li>

        <li>It uses the \c minLength value from the above two objects
        to determine the minimum of the two entries and gives
        preference to that entry.</li>

        <li>If \c index1 corresponds to the shortest EST while \c
        index2 corresponds to the longest EST (or vice versa) and
        there are more than two entries in the parameterSets vector
        then this method returns -1 (to indicate that the ESTs must
        not be compared).</li>

        </ol>

        \param[in] index1 The index of the first entry in the
        parameterSets vector to be used for comparison.

        \param[in] index2 The index of the second entry in the
        parameterSets vector to be used for comparison.

        \return The preferred of the two indexs if the fragments that
        fit into these entries are comparable.  Otherwise this method
        returns -1 to indicate the entries are not comparable.
    */
    int getPreferredSetIndex(const int index1, const int index2) const;

    /** Setup a convenient reference to the shared cDNA list.

         This method is called by the HeuristicChain just before
         initializing this class to setup a convenient reference to
         the shared RuntimeContext::estList.
         
         \param[in] list A pointer to the shared cDNA list.  This is
        the same as RuntimeContext::estList.
    */
    void setESTList(const ESTList* list);

	/** Helper method to reset lookup-table data structure.

		This helper method is invoked from the initialize() method and
		the destructor to free-up memory allocated for various
		internal data structures.
	*/
	void resetLookupTable();

	/** Helper method to reset ParameterSetmanager::parameterSets data
		structure.

		This helper method is invoked from the setupParameters()
		method and the destructor to free-up memory allocated for
		various internal data structures.
	*/
	void resetParameterSets();
	
private:
    /** The constructor.
	
        This is made private because the parameter set manager is
        currently meant to be instantiated only from the
        HeuristicChain.  This constructor is typically used to setup
        parameters for analyzing a hybrid data set. Specifically, it
        initializes the parameter set list to default values.  The
        default is for 3 parameter sets, one for each of the following
        "classes" of sequences:

        <table>
        <tr><th>Type</th><th>Min Len</th><th>Max Len</th><th>Notes</th></tr>

        <tr><td>Short</td><td>0</td><td>150</td><td>In practice range
        is 50 to 150, as the LengthFilter will filter out sequences of
        fewer than 50 bases long</td></tr>
        
        <tr><td>Medium</td><td>150</td><td>400</td><td></td></tr>

        <tr><td>Long</td><td>400</td><td>no limit</td><td>The
        practical upper limit is really 65K</td></tr>
        </table>

        \note The list of entries must be organized in ascending order
        of fragment lengths.  That is, the shorter fragment entries
        must occur first in the list.  The fragment lengths may not
        overlap.
        
        In the future this class should incorporate command-line
        arguments and customizability for the parameter sets.
        However, for now the defaults work just fine for a broad range
        of data sets.        
    */
    ParameterSetManager(const ESTList* estList,
                        int t1 = 20, int u1 = 4, int ws1 = 4,
                        int t2 = 35, int u2 = 6, int ws2 = 8,
                        int t3 = 50, int u3 = 8, int ws3 = 8);
    
    /** The constructor.
        
        This is made private because the parameter set manager is
        currently meant to be instantiated only from the
        HeuristicChain.  This constructor is typically used to setup
        parameters for analyzing a simple data set consisting purely
        of sanger data.

        <p>Specifically, this version of the overloaded constructor
        must be used to setup the parameters to a single entry.  The
        entry is chosen to be optimal for longer reads.  Conseuqently,
        this option should be used when it is known that the data is
        from sanger sequences that typically have longer reads.
        However, the single entry enables the clustering to run
        faster.</p>

		<p><b>Note:</b> This method is primarily meant to be used to
		configure parameters for longer reads.  It does cause the
		clustering to run faster, but should not be used with short
		reads.</p>
    */
    ParameterSetManager();
    
    /** The vector containing the list of parameter sets.

        This vector contains the list of parameter sets that have been
        added via the addParameterSet() method in this class.

        \note The list of entries must be organized in ascending order
        of fragment lengths.  That is, the shorter fragment entries
        must occur first in the list.  The fragment lengths may not
        overlap.
    */
    std::vector<ParameterSet*> parameterSets;

    /** A 2-D array to rapidly look-up the index of the
        parameter set to be used for given pair of fragments.

        <p>Assume that the parameterSet vector contains \e N entries.
        This array is initialized to a 2-D array of \e N+1 x \e N+1
        entries in the initialize() method.  Each value in this matrix
        indicates the entry in the parameterSets vector to be used for
        a given pair of ESTs.  The last row and last column are
        reserved for invalid entries that don't fit into any
        parameterSet entry.  The look-up is performed by first mapping
        the length of the individual ESTs to appropriate entries in
        the parameterSets array.</p>

        <p>To enable rapid mapping of fragment-length to
        parameterSets-index, the these mappings are pre-computed in
        the initialize() method and stored in the lenghtToParamSetMap
        vector.</p>

        <p>Once the aforementioned mappings are computed the look up
        into this array is performed (in the getParameterSet() method)
        as follows:

        \code

        return lookupTable[lengthToParamSetMap[estIdx1]][lengthToParamSetMap[estIdx2]];

        \endcode

        </p>
    */
    int **lookupTable;

    /** Vector to map index of fragment to corresponding entry in
        paramSets array.
        
        To enable rapid mapping of fragment-length to
        parameterSets-index, a set of mappings are pre-computed in the
        initialize() method and stored in this vector.
    */
    std::vector<int> lengthToParamSetMap;

    /** A convenient reference to the shared list of cDNA fragments.

        This pointer provides a convenient shortcut to the shared list
        of cDNA fragments being currently analyzed.  This pointer is
        the same as RuntimeContext::estList.  It is set by the
        HeuristicChain just before initializing this class.
    */
    const ESTList* estList;
};

#endif
