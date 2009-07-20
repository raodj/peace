#ifndef PMST_HEAP_CACHE_H
#define PMST_HEAP_CACHE_H

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
//          James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "MSTHeapCache.h"

/** A version of the MSTHeapCache that is specialized to work with
    the Partitioned MST Cluster Maker.  It provides some special
    functionality needed to construct the full MST using Kruskal's
    algorithm after merging MSTs that were constructed in parallel.

    This class is not designed for use by any subclasses and should
    not be extended.
*/
class PMSTHeapCache : public MSTHeapCache {
public:
    /** The constructor.
        
        This MST cache requires information about the number of ESTs
        in the system in order to effectively track the ESTs already
        added to the MST.  These values must be passed in as the
        various parameters to this class.  The constructor essentially
        passes off the parameters to the base class that saves the
        information (for future reference) in suitable instance
        variables.

        \note This object is typically instantiated in the
        PMSTClusterMaker::mergeManager() method.

        \param[in] totalESTCount The total number of ESTs to be
        processed (that is, to be added to the MST).

        \param[in] startOwnESTidx The starting index of the EST that
        is owned by the process that is using this cache.  This value
        is used to internally normalize est index values to 0 for the
        first EST owned by this process.

        \param[in] numOwnESTs The number of ESTs owned by the process
        that is using this cache.  This information is used to reserve
        necessary space to hold SMLists for each EST owned by the
        manager/worker process.

        \param[in] analyzer The analyzer to be used for EST
        comparisons, comparing metrics, and ordering entries in this
        cache depending on the type of EST analyzer (whether it is
        similarity based or distance metric based).
    */
    PMSTHeapCache(const int totalESTCount, const int startOwnESTidx,
		  const int numOwnESTs, const ESTAnalyzer *analyzer);
    
    /** The destructor.
        
        The STL data structures used by this class manage all the
        memory operations for this class.  Consequently, the
        destructor does not have any special tasks to perform.  It
        here present merely to adhere to coding conventions.
    */
    virtual ~PMSTHeapCache() {}

    /** Obtains and removes the top-most similar entry from the MSTCache.

        This method simply uses the top-most entry in the heap and
        populates the parameters with the appropriate value.  Note
        that the parameters are initialized to -1, -1, -1.0f, and -1
        respectively.

        \param[out] srcESTidx The source EST index from where the
        similarity metric is being measured.  The srcESTidx is already
        present in the MST.

        \param[out] destESTidx The destination EST index that is the
        best choice to be added to the MST (based on the local
        information).

        \param[out] metric The similarity/distance metric between the
        srcESTidx and the destESTidx.

        \param[out] alignmentData The alignment data associated with
        the srcESTidx and the destESTidx.
    */
    virtual void popBestEntry(int& srcESTidx, int& destESTidx,
			      float& metric, int& alignmentData);
    
protected:
    // Currently this class does not have any protected instance
    // variables or methods for use.
    
private:
    // Currently this class does not have any private instance
    // variables or methods for use.
};

#endif
