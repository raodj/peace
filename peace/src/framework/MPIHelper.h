#ifndef MPI_HELPER_H
#define MPI_HELPER_H

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

/** \file MPIHelper.h

    \brief File and macros to enable conditional compilation of MPI
    code so that PEACE can be built with and without MPI.

    This header file has been introduced to permit PEACE to compile
    with and without the presence of MPI. The code is designed in such
    a way to enable the use of MPI without hurting the readability and
    development of PEACE.
*/
#ifndef _WINDOWS
#include "config.h"
#endif

#include "Utilities.h"

#ifdef HAVE_LIBMPI

// We have MPI. In this case the standard/default operation of the
// system is good to go. Supress unused parameter warnings under GCC
// just for this header.

#if defined (__GNUG__) && !defined(ICC)
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <mpi.h>

#if defined (__GNUG__) && !defined(ICC)
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif

#endif

/** \def MPI_GET_RANK

    \brief A simple, convenience macro to conditionally call
    MPI::COMM_WORLD.Get_rank()

    <p>This macro provides a convenient wrapper around MPI's
    MPI::COMM_WORD.Get_rank() method to conditionally use it.  If MPI
    is available/enabled, this macro reduces to calling the actual
    Get_rank() method. However, if MPI is not available, this macro
    reduces to the constant 0 (one) making it appear as if there is
    only one process to work with. </p>
   
    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
        // ... some code goes here ..
	int myRank = MPI_GET_RANK();
        // ... more code goes here ..
    }

    \endcode
*/
#ifdef HAVE_LIBMPI
// We have MPI enabled
int  MPI_GET_RANK();
#else
// We don't have MPI
#define MPI_GET_RANK() 0
#endif

/** \def MPI_GET_SIZE

    \brief A simple, convenience macro to conditionally call
    MPI::COMM_WORLD.Get_size()

    <p>This macro provides a convenient wrapper around MPI's
    MPI::COMM_WORD.Get_size() method to conditionally use it.  If MPI
    is available/enabled, this macro reduces to calling the actual
    Get_size() method. However, if MPI is not available, this macro
    reduces to the constant 1 (one) making it appear as if there is
    only one process to work with. </p>
   
    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
        // ... some code goes here ..
	int WorkerCount = MPI_GET_SIZE();
        // ... more code goes here ..
    }

    \endcode
*/
#ifdef HAVE_LIBMPI
// We have MPI enabled
int MPI_GET_SIZE();
#else
// We don't have MPI
#define MPI_GET_SIZE() 1
#endif

#ifndef HAVE_LIBMPI
/** \def MPI_ANY_SOURCE

    \brief If MPI is not available, this macro defines a dummy
    MPI_ANY_SOURCE.
*/
#define MPI_ANY_SOURCE -1
#endif

/** \def MPI_STATUS

    \brief A simple, convenience macro to conditionally provide a
    suitable definition for MPI::Status data structure depending on
    whether MPI is available (or not).

    <p>This macro provides a convenient mechanism to work with
    MPI::Status data structure.  If MPI is available, then it defaults
    to MPI::Status.  On the other hand, if MPI is disabled then this
    macro provides a suitable definition so that the core code base
    does not get cluttered with unnecessary logic.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
        // ... some code goes here ..
		MPI_STATUS msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        // ... more code goes here ..
    }

    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_STATUS MPI_Status
#else
// We don't have MPI. So provide a suitable definition for MPI_Status
class MPI_STATUS {
public:
    int MPI_SOURCE = 0;
    int MPI_TAG    = 0;
};
#endif

/** Convenience method to extract the number of "top level" items.

    This is a convenience method that must be used to obtain the count
    of the number of items of a specific data type as indicated by the
    MPI_STATUS structure.

    \param[in] status The status structure from where the count is to
    be obtained.

    \param[in] mpiType The data type of the top-level item to be
    retrieved.

    \return The count of number of elements.  If MPI is not available
    then this method just returns 0 as convenience.
 */
#ifdef HAVE_LIBMPI
int MPI_GET_COUNT(MPI_STATUS& status, MPI_Datatype mpiType);
#else
int MPI_GET_COUNT(MPI_STATUS& status, int mpiType);
#endif

/** \def MPI_TYPE_INT

    \brief Macro to map MPI_TYPE_INT to MPI::INT (if MPI is enabled)
    or 0 if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::INT enumerated constant. If MPI is available,
    then MPI_TYPE_INT defaults to MPI::INT.  On the other hand, if MPI
    is disabled then this macro simply reduces to 0.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
        // ... some code goes here ..
	MPI_STATUS msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        const int dataSize = msgInfo.Get_count(MPI_TYPE_INT);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_TYPE_INT MPI_INT
#else
// MPI is not available
#define MPI_TYPE_INT 0
#endif

/** \def MPI_TYPE_CHAR

    \brief Macro to map MPI_TYPE_CHAR to MPI::CHAR (if MPI is enabled)
    or 0 if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::CHAR enumerated constant. If MPI is available,
    then MPI_TYPE_CHAR defaults to MPI::CHAR.  On the other hand, if
    MPI is disabled then this macro simply reduces to 0.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
        // ... some code goes here ..
	MPI_STATUS msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        const int dataSize = msgInfo.Get_count(MPI_TYPE_CHAR);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_TYPE_CHAR MPI_CHAR
#else
// MPI is not available
#define MPI_TYPE_CHAR 0
#endif

/** \def MPI_TYPE_2INT

    \brief Macro to map MPI_TYPE_2INT to MPI_2INT (if MPI is enabled)
    or 0 if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI_2INT enumerated constant. If MPI is available,
    then MPI_TYPE_2INT defaults to MPI_2TWO.  On the other hand, if
    MPI is disabled then this macro simply reduces to 0.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
        // ... some code goes here ..
	MPI_STATUS msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        const int dataSize = msgInfo.Get_count(MPI_TYPE_2INT);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_TYPE_2INT MPI_2INT
#else
// MPI is not available
#define MPI_TYPE_2INT 0
#endif

/** \def MPI_TYPE_FLOAT

    \brief Macro to map MPI_TYPE_FLOAT to MPI::FLOAT (if MPI is
    enabled) or 0 if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::FLOAT enumerated constant. If MPI is available,
    then MPI_TYPE_FLOAT defaults to MPI::FLOAT.  On the other hand, if
    MPI is disabled then this macro simply reduces to 0.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
        // ... some code goes here ..
	MPI_STATUS msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        const int dataSize = msgInfo.Get_count(MPI_TYPE_FLOAT);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_TYPE_FLOAT MPI_FLOAT
#else
// MPI is not available
#define MPI_TYPE_FLOAT 0
#endif

/** \def MPI_OP_SUM

    \brief Macro to map MPI_OP_SUM to MPI::SUM (if MPI is enabled)
    or 0 if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::SUM enumerated constant. If MPI is available,
    then MPI_OP_SUM defaults to MPI::SUM.  On the other hand, if
    MPI is disabled then this macro simply reduces to 0.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
	    // ... some code goes here ..
		int localCount = smList.size();
		int totalCount = 0;
		MPI_ALL_REDUCE(&localCount, &totalCount, 1, MPI_TYPE_INT, MPI_OP_SUM);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_OP_SUM MPI_SUM
#else
// MPI is not available
#define MPI_OP_SUM 0
#endif

/** \def MPI_OP_MAXLOC

    \brief Macro to map MPI_OP_MAXLOC to MPI::MAXLOC (if MPI is enabled)
    or 0 if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::MAXLOC enumerated constant. If MPI is available,
    then MPI_OP_MAXLOC defaults to MPI::MAXLOC.  On the other hand, if
    MPI is disabled then this macro simply reduces to 0.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    void someMethod() {
	    // ... some code goes here ..
		int localCount    = smList.size();
		int totalCount[2] = {0, 0};
		MPI_ALL_REDUCE(&localCount, &totalCount, 1, MPI_TYPE_2INT, MPI_OP_MAXLOC);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_OP_MAXLOC MPI_MAXLOC
#else
// MPI is not available
#define MPI_OP_MAXLOC 0
#endif

/** \def MPI_INIT

    \brief Macro to map MPI_INIT to MPI::Init (if MPI is enabled) or
    an empty method call if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::Init method. If MPI is available, then MPI_INIT
    defaults to MPI::Init.  On the other hand, if MPI is disabled then
    this macro simply reduces to a blank method.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..
	MPI_INIT(argc, argv);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_INIT(argc, argv) MPI_Init(argc, argv)
#else
// MPI is not available
#define MPI_INIT(argc, argv)
#endif

/** \def MPI_FINALIZE

    \brief Macro to map MPI_FINALIZE to MPI_Finalize (if MPI is
    enabled) or an empty method call if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::Finalize method. If MPI is available, then
    MPI_FINALIZE defaults to MPI::Init.  On the other hand, if MPI is
    disabled then this macro simply reduces to a blank method.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..
	MPI_INIT(argc, argv);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_FINALIZE() MPI_Finalize()
#else
// MPI is not available
#define MPI_FINALIZE()
#endif

/** \def MPI_WTIME

    \brief Macro to map MPI_WTIME to MPI::Wtime (if MPI is enabled) or
    to a suitable OS-independent implementation for Wtime.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::Wtime method. If MPI is available, then MPI_WTIME
    defaults to MPI::Wtime.  On the other hand, if MPI is disabled
    then this macro simply reduces to a suitable OS-independent method
    call.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..
	const double startTime = MPI::Wtime();
        // ... more code goes here ..
	const double elapsedTime = (MPI::Wtime() - startTime) * 1000.0;
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_WTIME MPI_Wtime
#else
// MPI is not available
extern double MPI_WTIME();
#endif

/** \def MPI_CODE

    \brief Macro to work around warnings about unused MPI variables or
    conditionally compile-out MPI code.

    <p>This macro provides a convenient, conditionally defined macro
    to conditionally compile-out MPI related code when MPI is
    disabled.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..
	MPI_CODE({
	    MPI_STATUS msgInfo;
            MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        });
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_CODE(x) x
#else
// MPI is not available
#define MPI_CODE(x)
#endif

/** \def MANAGER_RANK

    \brief Macro to refer to the rank of the manager process.


    <p>This macro provides a convenient mechanism to refer to the
	manager process while avoiding magic numbers being spread out
	in the code. Manager process is typically the process with
	rank 0 and is used to coordinate various distributed tasks
	and generate outputs.</p>

    This macro can be used as shown below:

    \code

    #include "MPIHelper.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..

		if (MPI_GET_RANK == MANAGER_RANK) {
		   // Do something special
		}
        // ... more code goes here ..
    }
    \endcode
*/
const int MANAGER_RANK = 0;

/** Determine the owner process Rank for a given estIdx.
    
    This method is a convenience method to determine the Rank of
    the process that logically owns a given EST.  The owning
    process is responsible for maintaining the cache for a given
    EST.  The owners are assigned in a simple fashion and ESTs are
    evenly divided up amongst all the processes.

    \param[in] estListSize The number of entries currently being
    processed.
    
    \param[in] estIdx The index of the EST whose owner process's rank
    is requested.  It is assumed that the estIdx is valid -- i.e., 0
    <= estIdx < estListSize.  If invalid EST index values are supplied
    then the operation of this method is undefined.
    
    \return The rank of the owner process for the given estIdx.
*/
int getOwnerProcess(const int estListSize, const int estIdx);

/** Helper method to compute the start and ending indexes of the
    EST that this process owns.
    
    This method was introduced to keep the math and logic clutter
    involved in computing the list of owned ESTs out of the
    methods that use the information.  This method returns the
    range, such that: \c startIndex <= \em ownedESTidx < \c
    endIndex.

    \param[in] estListSize The number of entries to be evely
    subdivided to determine the start and end index.
    
    \param[out] startIndex The starting (zero-based) index value
    of the contiguous range of ESTs that this process owns.
    
    \param[out] endIndex The ending (zero-based) index value of
    the contiguous range ESTs that this process owns.  The value
    returned in this parameter is \b not included in the range of
    values.
*/
void getLocallyOwnedESTidx(const int estListSize,
                           int& startIndex, int& endIndex);

#endif
