#ifndef CACHED_EST_INFO_HELPER_H
#define CACHED_EST_INFO_HELPER_H

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

#include "CachedESTInfo.h"

/** \file CachedESTInfoHelper.h

    \brief A collection of functors and static utility methods for
    varions operations on CachedESTInfo objects.

    This file was primarily introduced to try and centralized some of
    the utility methods and definitions associated with manaing
    CachedESTInfo objects.  These classes and methods were moved here
    as they technicaly don't belong to any specific class and they
    were unnecessarily cluttering other files.
*/

/** Functor for CachedESTInfo sorting.
    
    This Functor is used when sorting ESTs based on
    similarity/distance metric by the sort method defined in this
    class.  This functor sorts methods in ascending order with the
    smallest values at top.
*/
class LessCachedESTInfo :
    public std::binary_function<CachedESTInfo, CachedESTInfo, bool> {
public:
    /** Constructor.
        
        The constructor requires a pointer to the ESTAnalyzer that
        is being used for analysis.  The analyzer is used to
        compare the metrics for sorting CachedESTInfo objects.
        This constructor is called in the MSTCache class for
        sorting cache entries.
        
        \param[in] analyzer The analyzer to be used for comparing
        the metric values associated with two CachedESTInfo
        objects.
    */
    LessCachedESTInfo(const ESTAnalyzer *analyzer) :
        comparator(analyzer) {}
    
    /** \fn operator()
        
        \brief operator() for CachedESTInfo.
        
        The following operator provides a convenient mechanism for
        comparing CachedESTInfo objects for sorting.  This
        operator overrides the default STL operator() by comparing
        only the metrics of two CachedESTInfo objects (ignoring
        the EST indexes and any other information).  This method
        uses the metric comparison functionality provided by all
        EST analyzers.
        
        \param[in] estInfo1 The first CachedESTInfo object to be
        used for comparison.
        
        \param[in] estInfo2 The second CachedESTInfo object to be
        used for comparison.
        
        \return This method returns \c true if \c estInfo1.metric
        is better than \c estInfo2.metric.  Otherwise it returns
        \c false.
    */
    inline bool operator()(const CachedESTInfo& estInfo1,
                           const CachedESTInfo& estInfo2)
    { return comparator->compareMetrics(estInfo1.metric, estInfo2.metric); }
    
private:
    /** The functor for comparing.
        
        This comparator is set based on the template parameter
        associated with the enclosing class.
    */
    const ESTAnalyzer *const comparator;
};

/** Functor for CachedESTInfo sorting.
    
    This Functor is used when sorting ESTs based on
    similarity/distance metric by the sort method defined in this
    class.  This functor sorts methods in \b descending order with the
    smallest values at the bottom.  This functor is useful for
    maintaing heaps with the best metric at the top.
*/
class GreaterCachedESTInfo :
    public std::binary_function<CachedESTInfo, CachedESTInfo, bool> {
public:
    /** Constructor.
        
        The constructor requires a pointer to the ESTAnalyzer that is
        being used for analysis.  The analyzer is used to compare the
        metrics for sorting CachedESTInfo objects.  This constructor
        is called in the MSTCache derived classes for sorting cache
        entries.
        
        \param[in] analyzer The analyzer to be used for comparing
        the metric values associated with two CachedESTInfo
        objects.
    */
    GreaterCachedESTInfo(const ESTAnalyzer *analyzer) :
        comparator(analyzer) {}
    
    /** \fn operator()
        
        \brief operator() for CachedESTInfo.
        
        The following operator provides a convenient mechanism for
        comparing CachedESTInfo objects for sorting.  This
        operator overrides the default STL operator() by comparing
        only the metrics of two CachedESTInfo objects (ignoring
        the EST indexes and any other information).  This method
        uses the metric comparison functionality provided by all
        EST analyzers.
        
        \param[in] estInfo1 The first CachedESTInfo object to be
        used for comparison.
        
        \param[in] estInfo2 The second CachedESTInfo object to be
        used for comparison.
        
        \return This method returns \c true if \c estInfo2.metric is
        better than \c estInfo1.metric.  Otherwise it returns \c
        false.
    */
    inline bool operator()(const CachedESTInfo& estInfo1,
                           const CachedESTInfo& estInfo2)
    { return comparator->compareMetrics(estInfo2.metric, estInfo1.metric); }
    
private:
    /** The functor for comparing.
        
        This comparator is set based on the template parameter
        associated with the enclosing class.
    */
    const ESTAnalyzer *const comparator;
};

#endif
