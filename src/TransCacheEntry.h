#ifndef TRANS_CACHE_ENTRY_H
#define TRANS_CACHE_ENTRY_H

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

#include "HashMap.h"
#include "MSTCache.h"

/** \def HashMap<int, CachedESTInfo> TransCacheMap

    \brief A shortcut to refer to a hash map of CachedESTInfo.

    This typedef provides a convenient shortcut to refer to a hash map
    containing the information to compute metrics using
    conditional-transitivity.  The key into this hash map is the index
    of the EST to which a CachedESTInfo corresponds to.  This index is
    the value in CachedESTInfo::estIdx.
*/
typedef HashMap<int, CachedESTInfo> TransCacheMap;

/** Class to encapsulate information needed to apply
    conditional-transitivity.

    This class has been introduced to provide a convenient intreface
    to encapsulate information required to apply
    conditional-transitivity between a given pair of ESTs.  This class
    contains metric information pertaining to a single EST whose index
    is identified by the estIdx instance variable in this class.  The
    metric information is stored as a set of CachedESTInfo objects
    that are organized into the following two hash maps:

    <ol>

    <li>The parentInfo hash map contains CachedESTInfo objects
    corresponding to the parent(s) of estIdx. These entries form one
    set of relationships necessary to apply conditional-transitivity
    with various entries in the peerInfo hash map.</li>

    <li>peerInfo hash map contains CachedESTInfo objects corresponding
    to other ESTs with which estIdx could \i potentially have
    conditional-transivity relationships.  However, the applicability
    of transitivity is established only when the need arises. </li>

    </ol>

    The TransCacheEntry is created by the TransMSTClusterMaker after
    an EST has been added to the MST. Once an EST has been added to
    the MST, the top subset of CachedESTInfo is broadcasted to all the
    processes.  Each process then creates a set of TransCaceEntry
    objects for all the newly added ESTs and updates the list of
    related ESTs in parentInfo and peerInfo hash maps in this class.
*/
class TransCacheEntry {
    // Friend declaration to enable TransMSTClusterMaker to
    // instantiate objects of this type.
    friend class TransMSTClusterMaker;
public: 
    /** Constructor.

        The one and only constructor for this class. The constructor
        has been made public to enable this class to be used in
        conjunction STL data structures.  Currently this constructor
        is called from the TransMSTClusterMaker::populateCaches()
        method.

        \param[in] estIdx The index of the EST to which this
        TransCacheEntry is associated with.  This value is simply
        copied into the corresponding instance variable.
    */
    TransCacheEntry(const int estIdx = -1);
    
    /** The destructor.

        Currently, this class does not use any dynamic memory and all
        the memory is automatically managed by the encapsulated
        objects. Consequently, the destructor does not have any
        specific tasks to perform.
    */
    virtual ~TransCacheEntry() {}

    /** Extract and add metric entries in this cache entry.

        This method extracts information pertaining to ESTs in the
        range startIndex <= otherESTidx < endIndex and the parentInfo
        hash map must not have an existing entry for
        reference.refESTidx.
        
        \param[in] reference The reference information with respect to
        which entries in the metrics list have to be processed.  Note
        that reference.estIdx must be equal to this.estIdx.
        
        \param[in] metrics The list of metrics from which necessary
        entries for this cache must be extracted and stored.  Note
        that in a distributed scenario, it is not necessary to store
        all the metrics because not all comparisons are performed by
        all the processes.
        
        \param[in] startIndex The index of the EST starting from which
        the owner process is responsible for computing metrics.
        
        \param[in] startIndex The index of the last EST before which
        the owner process is responsible for computing metrics.       
    */
    void addEntries(const CachedESTInfo& reference, const SMList& metrics,
                    const int startIndex, const int endIndex);

    /** Obtain an existing metric from this cache.

        This method may be used to obtain an estimated metric value
        between this.estIdx and otherESTidx, assuming a suitable entry
        is available in this cache.

        \param[in] otherESTidx The index of EST to which an estimated
        metric value is expected.  If otherESTidx == this.estIdx this
        method always returns \c true and sets metric (parameter) to 0
        (zero).

        \param[out] metric If a valid entry is found in this cache
        then this parameter is updated with the estimated value.
        Otherwise this parameter is unaltered by this method.

        \return This method returns \c true if an estimated entry
        could be determined using transitivity.  Otherwise this method
        returns false.
    */
    bool getMetric(const int otherESTidx, float& metric) const;
    
protected:   
    /** The index of the EST to which this trans cache entry pertains
        to.

        This instance variable maintains the index of the EST to which
        this trans cache entry pertains to.  This value is set after
        the object is instantiated and is never changed during the
        life time of this object.
    */
    int estIdx;

private:    
   /** The parent metrics related to estIdx EST.

        This hash map contains a set of CachedESTInfo that are used to
        determine one of the metrics needed for computing metrics
        using conditional-transitivity.  These are the reference
        entries with respect to which entries are mainted in the
        peerInfo hash map. Entries in this cache are added by the
        addEntries() method in this class.
   */
    TransCacheMap parentInfo;
    
    /** The actual cache of metrics related to estIdx EST.

        This hash map contains a set of CachedESTInfo that are used to
        compute metrics using conditional-transitivity.  Entries in
        this cache are added by the addEntries() method in this class.
    */
    TransCacheMap peerInfo;
};

#endif
