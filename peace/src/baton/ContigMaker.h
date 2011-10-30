#ifndef CONTING_MAKER_H
#define CONTIG_MAKER_H

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

#include "BatonAlignmentInfo.h"
#include "ESTCodec.h"
#include "ESTList.h"
#include "NtDistr.h"
#include "Contig.h"

/** \file ContigMaker.h

    \brief Class to encapsulate a set of related alignment information
    to construct a contig.

    <p>This class essentially contains the set of related cDNA
    fragments that are to be combined together to form a contig. This
    class is primarily used by the BatonAssembler to build both local
    (on each parallel process) and global (the final one is built by
    the BatonAssemblerManager) contigs.</p>

	<p>The contigs are built using the BatonAlignmentInfo objects
	created by a SequenceAligner (typically the
	DefaultSequenceAligner). Each BatonAlignmentInfo object contains a
	list of Segment objects. Each Segment object provides information
	about how a given cDNA fragment must be aligned with-respect-to
	the reference cDNA fragment set when this ContigMaker was
	created.</p>

	<p>The contig created by this ContigMaker grows (to the left and
	right of the reference sequence) as new BatonAlignmentInfo objects are
	added to the ContigMaker. Intermediate consensus sequences are
	created to aid growing the contig further to the
	left-and-right. All consensus sequences (including the final
	contig) are created using the most frequently occuring base-pairs
	for each nucleotide position in the contig.</p>

	<p>In a parallel assembly run, the ContigMaker also performs the
	task of interacting with each other to build the global
	contig. Specially, the global contig generation proceeds in the
	following manner:

	<ol>

	<li>The BatonAssemblerManager along with any BatonAssemblerWorker
	classes collectively choose a cDNA fragment (through a distributed
	voting scheme) as the reference fragment.</li>

	<li>Each process (including BatonAssemberManager and
	BatonAssemberWorker) create a ContigManer with the reference
	sequence.</li>

	<li>Next, each process uses a SequenceAligner class to
	determine matching cDNA fragments in its purview (each process
	operates on a sub-set of cDNA fragments).  Matching cDNA fragments
	are individually aligned (using batons) to the reference sequence
	and the resulting Segments (encapsulated in a BatonAlignmentInfo
	object) are added to the ContigMaker </li>

	<li>Once all the local-cDNA fragments have been checked for
	matches, the ContigMakers build a local-intermediate contig and
	distribute it to the ContigMaker on the BatonAssemblerManager.
	The manager then builds a global intermediate-consensus and
	distributes it back to all BatonAssemblerWorkers.</li>

	<li>The above three steps are repated until no new cDNA fragments
	are added to the global consensus sequence.</li>

	<li>The final contig is built by the BatonAssemblerManager (when a
	matching cDNA fragment is no longer found by any of the processes)
	and the contig information is appropriately written.</li>

	<li>The above steps are repeated until all the cDNA fragments have
	been processed.</li>
	
	</ol>

	</p>

    <p>Note that for each set of related (and aligned) sequences, a
    new contig maker object is expected to be created.  The contig
    maker is not responsible for analyzing and computing the necessary
    assembly information.  It only performs the final stages of contig
    formation by utilizing the assembly information provided by other
    components of the baton assembler.  Refer to the various methods
    in this class for further documentation and details on this class.
    </p>
*/
class ContigMaker {
public:
    /** Constructor to start with an actual cDNA fragment as root/seed
        for sub-contig construction.

        The constructor is relatively straightforward. It initializes
        internal instance variables to their default initial values.
        In addition, if addRoot flag is true, then it adds the given
        root fragment as the first entry in the set of fragments to be
        aligned together. The global, logical starting position for
        this contig is set to be zero.

        \param[in] rootESTidx The index of the root (aka first)
        fragment constituting the contig being assembled.  This value
        must be int the range: 0 \f$ \le \f$ rootESTidx &lt
        EST::getSequenceCount().

        \param[in] estList A pointer to the cDNA fragments that are
        being currently being processed.  This list is not modified by
        this class.  It is used purely to access the cDNA nucleotide
        sequences. This pointer cannot be NULL.

        \param[in] tallyRoot If this flag is true then the root EST
        nucleotides are tallied into the final consensuses sequence
        for this contig. This flag is typically set to true only in
        the manager (and not by worker processes) to avoid
        multiple-inclusion of the root cDNA fragment.
    */
    ContigMaker(const int rootESTidx, const ESTList* estList,
                const bool addRoot);

    /** Constructor to create contigs using intermediate consensus
        sequence to build a sub-contig to left and/or right of another
        sub-contig.

        The constructor is relatively straightforward. It initializes
        internal instance variables to their default initial values.

        \param[in] baseSegIdx The index of the cDNA fragment
        from where the rootSeq was extracted. This information is
        used to create a suitable base segment entry.
        
        \param[in] rootSeq The cDNA fragment to be used as the
        root/seed for sub-contig construction. It is assumed that
        unprocessed cDNA fragments are compared (using batons) against
        this cDNA fragment.

        \param[in] rootSeqParentPos This value indicates the logical
        position of the rootSeq (parameter) in the parent contig. The
        parent contig is the ContigMaker from where the rootSeq was
        extracted. This value indicates the position of the rootSeq in
        the parent-contig (contig from where the rootSeq was
        extracted). This value is typically negative for contigs to
        the left and positive for contigs to the right. This value is
        used to finally merge this sub-contigs back into the parent
        contig (see ContigMaker::merge(const ContigMaker&) method).
        
        \param[in] estList A pointer to the cDNA fragments that are
        being currently being processed.  This list is not modified by
        this class.  It is used purely to access the cDNA nucleotide
        sequences. This pointer cannot be NULL.

        \param[in] tallyRoot If this flag is true then the root EST
        nucleotides are tallied into the final consensuses sequence
        for this contig. This flag is typically set to true only in
        the manager (and not by worker processes) to avoid
        multiple-inclusion of the root cDNA fragment.        
    */
    ContigMaker(const int baseSegIdx, const std::string& rootSeq,
                const int rootSeqParentPos,
                const int startPos, const int endPos,
                const ESTList* estList, const bool tallyRoot);
    
    /** The destructor.

        The destructor does not have any specific tasks to perform (as
        this class currently does not \e directly use any dynamic
        memory) and is present merely to adhere to coding conventions.
    */
    ~ContigMaker();

    /** Add a new fragment and its alignment information to this
        contig.

        This method must be used to add a new fragment to be assembled
        as a part of this contig.  The alignment information
        associated with the contig is the only parameter to this
        method.  This method performs the following tasks:

        <ol>

        <li>This method invokes the ensureCapacity() method to verify
        that there is sufficient room in the ntDistributions vector to
        hold the newly added alignment information. The necessary
        space is computed using the first and last segment information
        in the alignment information.</li>
        
        <li>The newly provided alignment information is added to the
        list of fragments to be assembled.</li>

		<li>The ntDistributions information is updated using the
		Segments in the alignment information object passed to this
		method.</li>
        
        <li>The left-most fragment information is checked and
        updated.</li>

        <li>The right-most fragment information is checked and
        updated.</li>
        
        </ol>

        \param[in] info The alignment information about the new
        fragment to be added to the list of fragments that are to be
        fused together to form the final contig.
    */
    void add(const BatonAlignmentInfo& info);

    /** Method to trigger global contig formation and to indicate if
		further operations to grow the contig to left and/or right is
		needed.

        This method must be used to initiate the final phase of contig
        formation (on all parallel processes) once all the pertinent
        fragments have been added to this contig maker.  This method
        performs the following tasks:

        <ol>

		<li>Process with MPI-rank zero is assumed to be the manager
		and other processes operate as workers (see
		MANAGER_RANK).</li>

		<li>Each worker-process sends its contig information to the
		manager-process. The contig information includes: lcal
		nucleotide occurrence frequencies (namely data in
		ntDistributions vector) along with associated counters.</li>

		<li>The manager-process receives the local-contig information
		from various worker-processes and merges the contig
		information into its local-contig information.  The operation
		of combining contig information parallel processes performed
		via the merge() method. Once data from all workers has been
		merged, the global contig is ready at the manager-process</li>

		<li>Next, the manager-process sends the global contig
		information to all the workers. The sending and receiving of
		data is performed via the sendContigData() and
		recvContigData() helper methods.</li>
		
        </ol>

		\note The alinedESTs data remains unaltered and remains
		local. We don't send and receive alined-segments as they are
		not really needed to form global consensus. However, as a
		side-effect, whenever segment information is needed (for
		instance writing detailed information to a file), then all
		parallel processes need to cooperatively generate the
		information.

		\return This method returns the total number of ESTs that were
		added by summing up the fragments added by each parallel
		process. If this method returns zero, then none of the
		parallel processes added any fragments. A zero return value
		also signifies that this contig cannot be grown further.
    */
    int formGlobalContig();

    /** Obtain the EST that is currently being used as the reference
        cDNA fragment by this contig maker.

        This method is a convenience method that can be used to
        determine the cDNA fragment this is currently being used by
        this Contig Maker to create/extend the contig.  The reference
        EST may change after every cycle of contig growing.

        \return A pointer to the reference cDNA that is being used by
        this contig maker to create/extend the contig.  If an actual
        EST is not the root (that is rootESTidx is -1 because an
        intermediate reference sequence is being used) then this method
        returns NULL.
    */
    inline const EST* getReferenceEST() const
    { return (rootESTidx != -1 ? estList->get(rootESTidx) : NULL);  }

    /** Obtain the root sequence being used by this contig maker.

        The reference nucleotide sequence being used by this contig
        maker to identify matching cDNA fragments for assembly. The
        returned value is always the root cDNA fragments immaterial of
        whether it is from a real cDNA fragment or from an
        intermediate reference sequence.

        \return The nucleotide string sequence being currently used as
        the reference sequence.
    */
    inline std::string getReferenceSeq() const {
        return rootSeq;
    }

    /** Method to detect if this contig has a left-overhang.

        This method can be used to detect of this contig has a
        left-overhang.

        \note The return value from this method is meaningful only
        after the formGlobalContig() method has been called.

        \return If this contig has a left-overhang then this method
        returns true. Otherwise this method returns false.
    */
	inline bool hasLeftOverhang() const {
        return !leftOverhang.empty();
    }

    /** Method to detect if this contig has a right-overhang.

        This method can be used to detect of this contig has a
        right-overhang.

        \note The return value from this method is meaningful only
        after the formGlobalContig() method has been called.

        \return If this contig has a right-overhang then this method
        returns true. Otherwise this method returns false.
    */
	inline bool hasRightOverhang() const {
        return !rightOverhang.empty();
    }


    /** Convenience method to construct a ContigMaker object that can
        be used to grow this contig to the left or right.

        This is a convenience method that is used to create a
        ContigMaker object that is pre-seeded with the left or the
        right overhang.  The returned object can be used to grow the
        contig to the left or right. Note that extra consensus
        nucleotides are appropriately added (to the left or right)
        based on the following rules:

		<ol>

		<li>If overhang len is shorter than 1 window size add
		sufficient consensus nucleotides to make rootNtSeq to be 1 window
		long.</li>

		<li>Otherwise add 1/2 window long consensus nucleotides to the
		overhang.</li>
		
		</ol>
		
        \param[in] leftOverhang If this flag is true, then this method
        returns a ContigMaker object seeded with the left-overhang. If
        this flag is false, then the returned object is seeded with
        the right-overhang.

		\param[in] windowSize The average window size to be used for
		ensuring that the root EST (seeded in the returned
		ContigMaker) is sufficiently long.  If the overhang is shorter
		than this window size then the overhang is appropriately
		extended to be at least as long as the window size. On the
		other hang, if the overhang is longer than the window Size,
		then half-window-long additional nucleotides from the
		overlapping regions are included to obtain a good root cDNA
		fragment.

        \param[in] tallyRoot If this flag is true then the root EST
        nucleotides are tallied into the final consensuses sequence
        for this contig. This flag is typically set to true only in
        the manager (and not by worker processes) to avoid
        multiple-inclusion of the root cDNA fragment.        
        
        \return A contig maker object seeded with the left or the
        right overhang.
    */
    ContigMaker getOverhangContigMaker(const bool leftOverhang,
									   const int windowSize,
                                       const  bool tallyRoot) const;

	/** Method to merge a sub-contig (created from this ContigMaker)
        back into this contig.

        This method is used by the BatonAssembler::localAssembly()
        method to merge a sub-contig back into this contig. There are
        two pre-requisites for this merge operation to work
        successfully:

        <ul>

        <li>The sub-contig must have been created via a successful
        call to getOverhangContigMaker() method in this (parent)
        object.</li>

        <li>After creating a overhang contig, no further changes to
        this (parent) object must have occurred.</li>

        <li>The sub-contig must have not been merged into this
        class. Merging the same sub-contig more than once will have
        undesired side-effects.</li>
        
        </ul>

        After the merge this (parent) contig will have all the
        necessary information from the sub-contig, including aligned
        EST Segment information. Furthermore, the left and right
        overhangs are appropriately updated.  Nucleotide frequencies
        are appropriately merged.
        
        \param[in] subContig The sub-contig (created via the
        getOverhangContigMaker() method in this class) to be merged
        into this (parent) contig.
    */
	void merge(const ContigMaker& subContig);	

    /** Generates a formatted contents of this contig for convenient
        analysis.

        This method provides a textual view of the alignment
        information for this complete contig.  The display is handy
        for analyzing alignments and troubleshooting any issues.

		\note The alignment information printed by this method applies
		only to the alignments available on the local process.
		
        \param[out] os The output stream to which the contents of this
        contig are to be written.
    */
    void prettyPrint(std::ostream& os) const;

    /** Method to obtain the current consensus sequence.

        This method utilizes the nucleotide distributions in the
        ntDistributions array to create a consensus sequence for the
        contig being formed.  The operation of this method is
        relatively straightforward.  For each nucleotide position in
        the given range, this method determines the base that has the
        highest frequency (ties in frequencies are randomy resolved)
        and uses that nucleotide character as the consensus base.

        \note This method internally calls the overloaded
        getConsensus(const int, const int) const method with suitable
        parameters to build the consensus sequence for the full contig.
        
        \return This method returns a nucleotide consensus sequence that
        represents the contig assembled.
    */
	std::string getConsensus() const;

    /** Method to obtain the current nucleotide occurrence frequencies
        for this contig.

        This method is a convenience emthod that can be used to obtain
        the raw occurrence frequencies (the values used to create a
        consensus sequence) in the ntDistributions array for this
        contig.  The distributions are returned in the form of a
        string that can be readily displayed for troubleshooting or
        other purposes.

        \param[in] separator The string to be used to separate
        consecutive occurrence-frequency values in the returned
        string.
        
        \return This method returns a string with nucleotide
        occurrence frequencies separated by the given seperator
        string.
    */
    std::string getNtDistribution(const std::string& separator = ",") const;
    
	/** Convenience method to populate data from this contig maker
		into a generic Contig object.

		This is a helper method (that is typically invoked after
		contig creation is successful) that is used by the
		BatonAssembler to populate a generic (assembler-neutral)
		Contig object with various information.

		\note This method does not set the contig ID.
		
		\param[out] contig The contig into which the information (from
		this contig maker) is to be populated.

		\param[in] setConsensus If this flag is \c true, then this
		method pouplates the consensus sequence in the contig.

		\param[in] setNtDistributions If this flag is \c true, then
		this method populates the nucleotide distribution data into
		the contig.
	*/
	void populateContig(Contig& contig, const bool setConsensus,
						const bool setNtDistributions) const;
	
protected:
    /** Helper method to form a consensus sequence.

        This method utilizes the nucleotide distributions in the
        ntDistributions array to create a consensus sequence for the
        contig being formed.  The operation of this method is
        relatively straightforward.  For each nucleotide position in
        the given range, this method determines the base that has the
        highest frequency and uses that nucleotide character as the
        consensus base.

        \param[in] startPos The starting position in the
        ntDistributions vector from where the consensus sequence is to
        be constructed.

        \param[in] endPos The ending position in the ntDistributions
        vector up to where the consensus sequence is to be
        constructed.

        \return This method returns a nucleotide consensus sequence that
        represents the contig assembled.
    */
	std::string getConsensus(const int startPos, const int endPos) const;

    /** Add a new fragment and its alignment information to this
        contig.

        This method must be used to add a new fragment to be assembled
        as a part of this contig.  The alignment information
        associated with the contig is the only parameter to this
        method.  This method performs the following tasks:

        <ol>

        <li>This method invokes the ensureCapacity() method to verify
        that there is sufficient room in the ntDistributions vector to
        hold the newly added alignment information. The necessary
        space is computed using the first and last segment information
        in the alignment information.</li>
        
        <li>The newly provided alignment information is added to the
        list of fragments to be assembled.</li>

		<li>If tallyBases (parameter) is true, the ntDistributions
		information is updated using the Segments in the alignment
		information object passed to this method.</li>
        
        <li>The left-most fragment information is checked and
        updated.</li>

        <li>The right-most fragment information is checked and
        updated.</li>
        
        </ol>

        \param[in] info The alignment information about the new
        fragment to be added to the list of fragments that are to be
        fused together to form the final contig.

        \param[in] tallyBases If this flag is set to true, then the
        nucleotide sequence is used to update the occurrence
        frequencies in the ntDistributions vector.
    */
    void add(const BatonAlignmentInfo& info, const bool tallyBases);

    /** Helper method to perform multi-process collective operations
        required to form the global contig on a manager process.

        This method is invoked form the formGlobalContig() method to
        perform the necessary operations to form a global contig after
        a round of local-contig formation has completed.  This method
        performs the following tasks:

        <ol>

        <li>It gathers local nucleotide occurrence frequencies (namely
        data in ntDistributions vector) along with associated counters
        from all the workers into this contig maker.  The operation of
        combining contig information from workers is performed via the
        merge() method. The merge operations also tracks overhangs to
        the left and right</li>

        <li>Once all the data from various workers has been merged
        then this method broadcasts the complete global contig back to
        all the workers for their local use.</li>
        
        </ol>

        \note This method assumes that the manager process's rank is
        zero.
        
		\return This method returns the total number of ESTs that were
		added by summing up the fragments added by each parallel
		process. If this method returns zero, then none of the
		parallel processes added any fragments. A zero return value
		also signifies that this contig cannot be grown further.
    */
	int managerFormGlobalContig();

    /** Helper method to perform multi-process collective operations
        required to form the global contig on worker processes.

        This method is invoked form the formGlobalContig() method to
        perform the necessary operations to form a global contig after
        a round of local-contig formation has completed.  This method
        performs the following tasks:

        <ol>

        <li>It sends local nucleotide occurrence frequencies (namely
        data in ntDistributions vector) along with associated counters
        to the manager.</li>

        <li>It then gathers the final global contig information from
        the manager into this contig manager.</li>
        
        </ol>

        \note This method assumes that the manager process's rank is
        zero.
        
		\return This method returns the total number of ESTs that were
		added by summing up the fragments added by each parallel
		process. If this method returns zero, then none of the
		parallel processes added any fragments. A zero return value
		also signifies that this contig cannot be grown further.
    */    
    int workerFormGlobalContig(const int managerMPIrank = 0);
	
	/** Helper method to send contig data from this contig maker to
		another process.

		This method is used in conjunction with the recvContigData()
		method to exchange contig information between a manager
		process (process with MPI rank 0) and worker processes
		involved in parallel, baton assembly.  This method performs
		the following tasks:

		<ol>

		<il>First it send out counter information to the designated
		(procID parameter) process.</li>

		<li>Next it sends the overhang nucleotide sequences.</li>
		
		<li>Finally, it dispatches the data from the ntDistributions
		vector.</li>

		</ol>

		\param[in] procID The MPI-rank of the destination process to
		which the contig information is to be sent.

		\param[in] alignedESTsCount The count of number of aligned
		ESTs to be sent. This value is typically alignedESTs.size() on
		workers and the sum-of-aligned-ESTs on the manager.

		@see recvContigData
	 */
    void sendContigData(const int procID, const int alignedESTsCount);

	/** Helper method to receive contig data from this contig maker to
		another process.

		This method is used in conjunction with the sendContigData()
		method to exchange contig information between a manager
		process (process with MPI rank 0) and worker processes
		involved in parallel, baton assembly.  This method performs
		the following tasks (coordinated with the sequence in
		sendContigData() method):

		<ol>

		<il>First it receives counter information from the designated
		(procID parameter) process. This information is used in the
		next two steps to pre-allocate sufficient buffers to receive
		data.</li>

		<li>Next, it receives the overhang nucleotide sequences.</li>

		<li>It then reads the data for the ntDistributions vector.</li>
		
		<li>Finally, if the merge flag (second parameter) is \c true
		then the received data is merged (used by manager process to
		accumulate data from workers) into this contig. Otherwise
		(i.e., when merge parameter is \c false) the contig data is
		replaced with the received information (this happens at
		workers that receive full contig information from the manager
		and replace their local information).  </li>
		
		</ol>

		\note The alinedESTs data remains unaltered and remains
		local. We don't send and receive alined-segments as they are
		not really needed to form global consensus. However, as a
		side-effect, whenever segment information is needed (for
		instance writing detailed information to a file), then all
		parallel processes need to cooperatively generate the
		information.
		
		\param[in] procID The MPI-rank of the destination process from
		which the contig information is to be read.

		\param[in] merge Flag to indicate if the received data must be
		merged with existing data or if existing data must be
		completely replaced with recieved data.

		\return This method returns the value of number of aligned
		ESTs as reported in the received counter values.
		
		@see recvContigData
	*/
    int recvContigData(const int procID, const bool merge);
	
protected:
	/** Helper method to ensure ntDistributions vector has sufficient
		capacity to hold nucleotide occurrence frequencies.

		This method is a helper method that is used to ensure that the
		ntDistributions vector has sufficient capacity to hold
		nucleotide occurrence frequencies.  The minimum required
		capacity is specified as the parameters. However, to minimize
		repeated growth of ntDistributions vector (which can be
		detrimental to performance) this method grows the vector to
		accomodate additional nucleotides. The extra growth is
		controlled by the padding value.

		\note Whenever the ntDistributions vector is grown to the left
		(happens when leftNtPos <= 0), then this method adjusts the
		following internal instance variables: leftMostNtPos,
		rightMostNtPos, rootESTNtPos, leftOverhangEnd,
		rightOverhangStart.
		
		\param[in] leftNtPos The left-most nucleotide position with
		respect to the current rootESTNtPos that is needed. This value
		can be negative.

		\param[in] rightNtPos The right-most nucleotide position with
		respect to the current rootESTNtPos.

		\param[in] padding The extra padding to be included whenever
		the the ntDistributions array is grown. This value cannot be
		negative (but can be zero).  The padding is applied to the
		right and left growth independently (that is, if the padding
		is 1024 and ntDistributions is grown to the left-and-right,
		then the net padding applied is 2048).
	*/
    void ensureCapacity(const int leftNtPos, const int rightNtPos,
						const int padding = 1024);

	/** Helper method to update occurrence frequencies (in
		ntDistributions) of nucleotides for a given segment.

		This method is a helper method (that is invoked from the add()
		method) to update the occurrence frequencies of nucleotides
		for a given segment.  The position of the nucleotides is
		determined using the position with-respect-to the position of
		the root fragment (namely rootESTNtPos). 

		\note This method assumes sufficient capacity of
		ntDistributions vector has been verified via call to
		ensureCapacity() method.
		
		\param[in] seg The segment whose nucleotide occurrence
		frequencies are to be added to the ntDistributions vector.

		\param[in] othSeq The sequence of nucleotides for the EST
		whose occurrence frequencies are to be upadated (in the
		ntDistributions vector).
	*/
	void tally(const Segment& seg, const std::string& othSeq);

	/** Helper method to merge contig nucleotide information from a
		remote ContigMaker into this class.

		Recollect that in a parallel run, each parallel process
		operates on a sub-set of ESTs. Consequently, although all
		processes start with the same root EST, the contigs built by
		each one will have a different sub-set of ESTs.  Consequently,
		at the end of processes all the local-ESTs, this method is
		used to builda global contig from the sub-contigs created by
		various parallel processes. Since all the ContigMakers operate
		with the same root EST, this method uses the position of the
		root EST to appropriately merge the distribution data.

		\note This method appropriately updates the following internal
		instance variables: leftMostNtPos, rightMostNtPos,
		rootESTNtPos, leftOverhangEnd, rightOverhangStart.
		
		\note The ntDistributions arrray does not include the
		nucleotide occurrences for the root fragment (until the end)
		to avoid multiple-recounting of the root fragment's
		nucleotides when these merges happen.

		\param[in] otherNtDistributions The nucleotide distributions
		from another process whose information is to be appropriately
		merged into this object.

		\param[in] otherLeftMostNtPos The left-most index postion in
		the otherNtDistributions vector that has the first
		nucleotide. All values to the left of otherLeftMostNtPos in
		otherNtDistributions vector are ignored.

		\param[in] otherRightMostNtPos The right-most index position
		in the otherNtDistributions vector that has the last
		nucleotide. All values to the right of otherRightMostNtPos in
		otherNtDistributions vector are ignored.

		\param[in] otherRootESTNtPos The location of the common root
		cDNA fragment in the otherNtDistributions. This value serves
		as the anchor to merge otherNtDistributions data with
		ntDistributions data.

		\param[in] otherLeftOverhangEnd The left-overhang end position
		in the otherNtDistributions. This value is used to update the
		leftOverhangEnd value in this class after the nucleotide
		frequency data is merged.

		\param[in] otherLeftOverhangESTidx The index of the cDNA
        fragment from which the left overhang was extracted. This
        information is used to track and compute base segment
        information.
        
		\param[in] otherRightOverhangStart The right-overhang start
		position in the otherNtDistributions. This value is used to
		update the rightOverhangStart value in this class after the
		nucleotide frequency data is merged.

		\param[in] otherRightOverhangESTidx The index of the cDNA
        fragment from which the right overhang was extracted. This
        information is used to track and compute base segment
        information.
	*/
	void merge(const NtDistrList& otherNtDistributions,
			   const int otherLeftMostNtPos, const int otherRightMostNtPos,
			   const int otherRootESTNtPos,
               const std::string& otherLeftOverhang,
			   const int otherLeftOverhangESTidx,
			   const std::string& otherRightOverhang,
			   const int otherRightOverhangESTidx);

private:
    /** The list of fragments that have been added to form the contig.

        This array maintains the list of alignment information
        associated with the fragments to be fused together to form the
        final consensus sequence.  The first entry to this list is
        added in the constructor.  Subsequent entries are added each
        time the add() method in this class is called.
    */
    BatonAlignmentInfoList alignedESTs;

    /** The list of base segment information for the contig being
        formed.

        This vector maintains the list of base segments associated
        with this contig. This information is used to create the final
        assembler-neutral contig for generating results in various
        file formats. New entries are added and updated when contigs
        are merged.
    */
    BaseSegmentInfoList baseSegments;
    
    /** Occurrence frequencies of base pairs in each column of the
        contig.

        <p>This array is used to maintain the occurence frequencies of
        nucleotides for each column in the contig. This vector is
        typically grown (dynamically as needed) to accomodate extra
        500 nt (to the left and right) each time it needs to be grown.
        This done to keep the memory allocations to a minimum.  The
        instance variables leftMostNtPos and rightMostNtPos indicate
        the index-locations in this vector that are actually used.
        The overlapStartPos and overlapRightPos indicate the left-most
        and right-most positions where overlaps begin and end
        respectively.</p>

		Here is an example of the data that is stored in the
		ntDistributions given the various ESTs in the alignedESTs
		list.  In the following case, we assume that alignedESTs
		contains 4 ESTs (namely, \c EST-A, \c EST-B, \c EST-C, and \c
		EST-D whith \c EST-A being the root EST). The data in
		ntDistributions list shown in the table below:

		<table cellpadding="3" cellspacing="0">
        <tr><th>Fragment</th><th>Nucleotides</th></tr>

        <tr><td>EST-A (root):      </td><td>\c AACTAGGTACG</td></tr>
        <tr><td>EST-B:             </td><td>\c ---TAGGTA--</td></tr>
        <tr><td>EST-C:             </td><td>\c AACTCAGG----</td></tr>
        <tr><td>EST-D:             </td><td>\c ----TGGTACG</td></tr>
        
        <tr><td colspan="2"><b>Data in ntDistributions</b></td></tr>
        
        <tr><td>Nt-A (\cfreq[0]) : </td><td>\c 22002140300</td></tr>
        <tr><td>Nt-T (\cfreq[0]) : </td><td>\c 00031003000</td></tr>
        <tr><td>Nt-C (\cfreq[0]) : </td><td>\c 00201000020</td></tr>
        <tr><td>Nt-G (\cfreq[0]) : </td><td>\c 00000341002</td></tr>
        <tr><td>Nt-N (\cfreq[0]) : </td><td>\c 00000000000</td></tr>		
        </table>
    */
    NtDistrList ntDistributions;    

	/** Variable to track the index-position in the ntDistributions
	    list where the current, left-most nucleotide information is
	    stored.

		<p>The ntDistributions vector is dynamically grown (as new cDNA
		are added) to accommodate growth of the contig to the
		left-and-right of the reference cDNA. However, the
		ntDistributions vector is grown extra entries (to the left and
		right) to keep memory allocations to a minimum.  Consequently,
		index position zero does not indicate the left-most
		entry. Instead, this instance variable must be used to
		determine the logical index of the first nucleotide
		information in the ntDistributions vector.</p>

		<p>This instance variable is updated in conjunction with
		rightMostNtPos and rootESTNtPos whenever the ntDistributions
		array is grown by the ensureCapacity() method in this
		class.</p>

		\note The following fact always holds: 0 &#8804; leftMostNtPos
		< rightMostNtPos
	*/
	int leftMostNtPos;
	
	/** Variable to track the index-position in the ntDistributions
	    list where the current, right-most nucleotide information is
	    stored.

		<p>The ntDistributions vector is dynamically grown (as new
		cDNA are added) to accommodate growth of the contig to the
		left-and-right of the reference cDNA. However, the
		ntDistributions vector is grown with extra entries (to the
		left and right) to keep memory allocations to a minimum.
		Consequently, index position ntDistributions.size() does not
		indicate the usable right-most entry. Instead, this instance
		variable must be used to determine the logical index of the
		first nucleotide information in the ntDistributions
		vector.</p>

		<p>This instance variable is updated in conjunction with
		leftMostNtPos and rootESTNtPos whenever the ntDistributions
		array is grown by the ensureCapacity() method in this
		class.</p>

		\note The following fact always holds: leftMostNtPos &lt;
		rightMostNtPos &#8804; ntDistribution.size()
	*/	
	int rightMostNtPos;
	
	/** Variable to track the index-position in the ntDistributions
	    list to track the starting index of the root cDNA.

		<p>The ntDistributions vector is dynamically grown (as new
		cDNA are added) to accommodate growth of the contig to the
		left-and-right of the reference cDNA.  Consequently, the index
		position where the nucleotide information for the root cDNA
		fragment is stored changes. This variable tracks the logical
		starting position of the root cDNA fragment to ensure
		information from Segments are correctly merged to produce the
		firnal contig.</p>

		<p>This instance variable is updated in conjunction with
		leftMostNtPos and rightMostNtPos whenever the ntDistributions
		array is grown by the ensureCapacity() method in this
		class.</p>

		\note The following fact always holds: leftMostNtPos &#8804;
		rootESTNtPos &lt; rightMostNtPos
	*/	
	int rootESTNtPos;

	/** The length of the root cDNA fragment with which other cDNA
		fragments are currently being aligned.

		This instance variable is used to track the length of the cDNA
		fragment which is currently serving as the root fragment with
		which other fragments are being aligned.
	*/
	int rootESTLen;

    /** The index (in estList) of the current cDNA fragment that is
        serving as the seed for contig extension.

        This instance variable is used to track the current cDNA
        fragment that is being used for contig extension. The index is
        into the estList set in this class.  This value initialized in
        the constructor and may change when sub-contigs are merged.

        \note If a valid rootESTidx is not available then this value
        will be -1.
    */
    int rootESTidx;

    /** The nucleotide sequence for the root cDNA fragment for this
        contig maker.

        This instance variable tracks the nucleotide sequence for the
        root cDNA fragment.  This is the cDNA against which other cDNA
        constituting this contig have been aligned. This value can
        either be the sequence for the rootESTidx (if it is valid) or
        an intermediate sub-sequence (typically a left or a right
        overhang) from another contig.
    */
    std::string rootSeq;

    /** This value indicates the logical position of the root sequence
        (for this contig) in the parent contig. The parent contig is
        the contig from which this contig was created to explore the
        parent's left/right overhang.  This value indicates the
        position of the rootSeq in the parent-contig (contig from
        where the rootSeq was extracted). The parent contig starts at
        value zero.  Therefore, this value is typically negative for
        contigs to the left and positive for contigs to the
        right. This value is used to merge sub-contigs back into their
        parent contig.
    */
	int rootSeqParentPos;
    
	/** The left overhang nucleotide sequence for this contig
		maker.

		The left-overhang is a sub-fragment (or a complete fragment)
		that has the maximum overhang to the left of the root EST.
		The left-overhang begins at leftMostNtPos and ends at
		leftOverhangEnd.  This value is initialized to an empty
		sequence (\c "") just before a round of contig formulation
		commences.  This value remains an empty string if a left
		overhang is not found during contig extension. This value is
		tracked in by the add() method.
	*/
	std::string leftOverhang;

    /** The index (in estList) of the current cDNA fragment that is
        contributing to the leftOverhang for contig extension.

        This instance variable is used to track the current cDNA
        fragment that is contributing to the left overhang. The index
        is into the estList set in this class.  This value changed
        each time a new left overhang is detected or when sub-contigs
        are merged. If a left overhang is not available, then this
        instance variable has the value of -1. This variable is
        primarily used to track base segments for the contig.
    */    
    int leftOverhangESTidx;
    
	/** The right overhang nucleotide sequence for this contig maker.

		The right-overhang is a sub-fragment (or a complete fragment)
		that has a overhang (i.e., it has nucleotides that do not
		overlap with the root cDNA fragment) to the right of the root
		EST.  This value is initialized to an empty string (\c "")
		just before a round of contig formulation commences.  This
		value remains an empty string if a right overhang is not found
		during contig extension.  This value is tracked in by the
		add() method.
	*/
	std::string rightOverhang;

    /** The index (in estList) of the current cDNA fragment that is
        contributing to the rightOverhang for contig extension.

        This instance variable is used to track the current cDNA
        fragment that is contributing to the right overhang. The index
        is into the estList set in this class.  This value changed
        each time a new right overhang is detected or when sub-contigs
        are merged. If a right overhang is not available, then this
        instance variable has the value of -1. This variable is
        primarily used to track base segments for the contig.
    */
    int rightOverhangESTidx;
    
    /** A convenience pointer to the list of cDNA fragments being
        processed.

        This pointer is passed in via the constructor whenever a
        contig maker object is instantiated.  This pointer is used to
        refer to the cDNA fragments during contig formulation.
    */
    const ESTList* estList;
};

#endif
