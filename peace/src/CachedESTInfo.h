#ifndef CACHED_EST_INFO_H
#define CACHED_EST_INFO_H

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

/** A simple class to encapsulate data to be cached about ESTs.

    This class represents the final (atomic) entry in the MSTCache
    data structure.  This class is meant to serve just as an
    encapsulation class. Consequently, it is very straightfoward and
    exposes all the fields in it directly.

    \note This class is meant to be used only by the MSTCache.
*/
class CachedESTInfo {
public:
    /** The default constructor.

        This class is heavily used in conjunction with STL containers
        (such as std::vector). Consequently, the default constructor
        is called very frequently when entries are added to the
        MSTCache. Therefore the constructor is very straightfoward to
        the point that it does absolutely nothing and is merely
        present to adhere to coding conventions.
    */
    CachedESTInfo() {}

    /** A convenience all parameter constructor.

        This is a convenience constructor that is used to set up all
        the instance variables to corresponding values specified by
        the parameters.

        \param[in] estIdxValue The EST index value to be stored in
        this class.

        \param[in] metricValue The similarity/distance metric value to
        be stored in this object.

        \param[in] alignmentDataValue Additional alignment information
        to be stored in this class.
    */
    CachedESTInfo(const int estIdxValue, const float metricValue,
                  const int alignmentDataValue) :
      estIdx(estIdxValue), metric(metricValue),
      alignmentData(alignmentDataValue) {}
    
    /** The destructor.

        This class does not contain any members that utilize dynamic
        memory to hold data.  Consequently, the destructor has no
        special tasks to perform and is merely present to adhere to
        coding conventions.
    */
    ~CachedESTInfo() {}

    /** The index of the EST for which this data is relevant.

        This instance variable is used to store the index of one of
        the ESTs for which the remainder of the information has been
        stored.
    */
    int estIdx;

    /** The distance/similarity metric value.

        This instance variable is used to maintain the
        similarity/distance metric corresponding to the EST referred
        by estIdx.  Note that this class does not contain the index of
        the reference EST as this information is maintained in the
        MSTCache directly.
    */
    float metric;

    /** The alignment information for this EST.

        This instance variable is used to maintain the alignment data
        (if any) associated with the EST indicated by estIdx.
    */
    int alignmentData;

protected:
    // Currently this class does not have any protected members.
    
private:
    // Currently this class does not have any private members.
};

#endif
