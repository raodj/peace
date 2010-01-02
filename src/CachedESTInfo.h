#ifndef CACHED_EST_INFO_H
#define CACHED_EST_INFO_H

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

        \param[in] refESTidxValue The zero-based index of the EST that
        was used as the reference EST which which EST at index \c
        estIdxValue was compared/analyzed to generate the \c
        metricValue.
        
        \param[in] estIdxValue The zero-based index of the EST that
        was compared against the EST at index \c refESTidx to generate
        \c metricValue.
        
        \param[in] metricValue The similarity/distance metric value to
        be stored in this object.
        
        \param[in] alignmentDataValue Additional alignment information
        to be stored in this class.
    */
    CachedESTInfo(const int refESTidxValue, const int estIdxValue,
		  const float metricValue, const int alignmentDataValue) :
        refESTidx(refESTidxValue), estIdx(estIdxValue), metric(metricValue),
        alignmentData(alignmentDataValue) {}
    
    /** The destructor.

        This class does not contain any members that utilize dynamic
        memory to hold data.  Consequently, the destructor has no
        special tasks to perform and is merely present to adhere to
        coding conventions.
    */
    ~CachedESTInfo() {}

    /** The zero-based index of the reference EST.

        This instance variable is used to store the zero-based index
        of the reference EST against which the EST at \c estIdx was
        analyzed to generate the \c metric and \c alginmentData.
    */
    int refESTidx;
    
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
