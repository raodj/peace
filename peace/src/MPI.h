#ifndef MPI_H
#define MPI_H

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

/** \file MPI.h

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
// system is good to go.
#include <mpi.h>

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

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
	int myRank = MPI_GET_RANK();
        // ... more code goes here ..
    }

    \endcode
*/
#ifdef HAVE_LIBMPI
// We have MPI enabled
#define MPI_GET_RANK() MPI::COMM_WORLD.Get_rank()
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

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
	int WorkerCount = MPI_GET_SIZE();
        // ... more code goes here ..
    }

    \endcode
*/
#ifdef HAVE_LIBMPI
// We have MPI enabled
#define MPI_GET_SIZE() MPI::COMM_WORLD.Get_size()
#else
// We don't have MPI
#define MPI_GET_SIZE() 1
#endif

/** \def MPI_ANY_SOURCE

    \brief If MPI is not available, this macro defines a dummy
    MPI_ANY_SOURCE.
*/
#ifndef HAVE_LIBMPI
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

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
	MPI_STATUS msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        // ... more code goes here ..
    }

    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_STATUS MPI::Status
#else
// We don't have MPI. So provide a suitable definition for MPI_Status
class MPI_STATUS {
public:
    inline int Get_source() const { return 0; }
    inline int Get_count(const int UNREFERENCED_PARAMETER(type) = 0) const { return 0; }
    inline int Get_tag() const { return 0; }
};
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

    #include "MPI.h"

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
#define MPI_TYPE_INT MPI::INT
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

    #include "MPI.h"

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
#define MPI_TYPE_CHAR MPI::CHAR
#else
// MPI is not available
#define MPI_TYPE_CHAR 0
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

    #include "MPI.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..
	MPI_INIT(argc, argv);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_INIT(argc, argv) MPI::Init(argc, argv)
#else
// MPI is not available
#define MPI_INIT(argc, argv)
#endif

/** \def MPI_FINALIZE

    \brief Macro to map MPI_FINALIZE to MPI::Finalize (if MPI is
    enabled) or an empty method call if MPI is unavailable.

    <p>This macro provides a convenient, conditionally defined macro
    to refer to MPI::Finalize method. If MPI is available, then
    MPI_FINALIZE defaults to MPI::Init.  On the other hand, if MPI is
    disabled then this macro simply reduces to a blank method.</p>

    This macro can be used as shown below:

    \code

    #include "MPI.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..
	MPI_INIT(argc, argv);
        // ... more code goes here ..
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_FINALIZE() MPI::Finalize
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

    #include "MPI.h"

    int main(int argc, char *argv[]) {
        // ... some code goes here ..
	const double startTime = MPI::Wtime();
        // ... more code goes here ..
	const double elapsedTime = (MPI::Wtime() - startTime) * 1000.0;
    }
    \endcode
*/
#ifdef HAVE_LIBMPI
#define MPI_WTIME MPI::Wtime
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

    #include "MPI.h"

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

#endif
