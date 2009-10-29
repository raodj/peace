#ifndef MPI_STATS_H
#define MPI_STATS_H

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

#include <iostream>
#include <string>

// Include our non-native MPI header to enable/diable MPI usage.
#include "MPI.h"

/** A non-instantiable wrapper class to track MPI-related statistics.

    This class is a simple wrapper class that is used to encapsulate
    the statistic gathered during program execution.  The statistics
    collated by this class related to MPI calls.  In order to
    facilitate statistics collection without cluttering MPI code, this
    class provides a set of handy macros for various MPI macro calls,
    namely: MPI_PROBE, MPI_SEND, MPI_RECV, MPI_BCAST, and
    TRACK_IDLE_TIME.  Use the aforementioned macros to streamline your
    code.  Finally, at the end of MPI processing use
    MPIStats::displayStats() method do display the statistics gathered
    by this method.
*/
class MPIStats {
public:
    /** Cumulative time spent idling.

        This \b static instance variable (initialized to 0) maintains
        the cumulative time spent by this process idling, waiting for
        messages to arrive from another process.  This value is
        tracked by wrapping each MPI probe call (or MPI recv as the
        case may be) with the TRACK_IDLE_TIME macro.  The time is
        tracked using the Wtime() method defined by MPI.
    */
    static double idleTime;

    /** Number of messages sent by this process.

        <p> This \b static instance variable (initialized to 0 (zero))
        is used to essentially track the number of calls to MPI send
        method.  This variable must be incremented each time MPI's
        send method is invoked.  This variable is finally printed in
        the displayStats() method.<p>

        \note Using the MPI_SEND macro call streamlines this process
        as it increments the send count as well.
    */
    static int sendCount;

    /** Number of messages received by this process.

        <p>This \b static instance variable (initialized to 0 (zero))
        is used to essentially track the number of calls to MPI recv
        method.  This variable must be incremented each time MPI's
        recv method is invoked.  This variable is finally printed in
        the displayStats() method.</p>

        <p>If you would also like to track the time spent to recv a
        certain message as idle time then you can easily achieve this
        objective as shown in the code fragement below:</p>

        \code
        #include "MPI.h"

        void someMethod() {
            TRACK_IDLE_TIME(MPI_SEND(b, 1, MPI_TYPE_INT, rank, tag));
        }
        \endcode

        \note Using the MPI_RECV macro call streamlines this process
        as it increments the recvCount as well.
    */
    static int recvCount;

    /** Number of messages broad casted by this process.

        This \b static instance variable (initialized to 0 (zero)) is
        used to essentially track the number of calls to MPI bcast
        method.  This variable must be incremented each time MPI's
        bcast method is invoked.  This variable is finally printed in
        the displayStats() method.

        \note Using the MPI_BCAST macro call streamlines this process
        as it increments the bcastCount as well.
    */
    static int bcastCount;

    /** Number of times reduction operation was performed by this process.

        This \b static instance variable (initialized to 0 (zero)) is
        used to essentially track the number of calls to MPI reduce
        method.  This variable must be incremented each time MPI's
        reduce method is invoked.  This variable is finally printed in
        the displayStats() method.

        \note Using the MPI_REDUCE macro call streamlines this process
        as it increments the reduceCount as well.
    */
    static int reduceCount;
	
    /** Number of time this process probed for messages.

        This \b static instance variable (initialized to 0 (zero)) is
        used to essentially track the number of calls to MPI probe
        method.  This variable must be incremented each time MPI's
        bcast method is invoked.  This variable is finally printed in
        the displayStats() method.

        \note Using the MPI_PROBE macro call streamlines this process
        as it increments the probeCount as well.
    */
    static int probeCount;

    /** Displays statistics collated thus far.

        This method is invoked at the end of the makeClusters() method
        to display performance related statistics for this class.
        This method displays the following information:

        <ul>

        <li> The time spent idling waiting for messages to arrive (in
        seconds).  This information is tracked in the idleTime
        instance variable. </li>

        <li> The number of messages that were sent (essentially calls
        to MPI Send method).  This information is tracked in the
        sendCount instance variable. </li>

        <li> The number of messages that were received (essentially
        calls to MPI Recv method).  This information is tracked in the
        recvCount instance variable. </li>

        <li> The number of calls to broadcast (essentially calls to
        MPI Bcast method).  This information is tracked in the
        bcastCount instance variable. </li>

        <li> The number of reduction operations (essentially calls to
        MPI Reduce method).  This information is tracked in the
        reduceCount instance variable. </li>

       <li> The number of calls to proble (essentially calls to MPI
        probe method).  This information is tracked in the probeCount
        instance variable. </li>
        
        </ul>

        \param[out] os The output stream to which all the information
        must be written.
    */
    static void displayStats(std::ostream& os,
                             const std::string& indent = "    ");    
protected:
    // Currently this class has no protected members.

private:
    /** The default constructor.

        The default constructor has is private in order to ensure that
        this class is never instantiated.  Instead, the static methods
        and static instance variables in this class are directly used.
    */    
    MPIStats() {}

    /** The destructor.

        The destructor is private to ensure that objects of this class
        can never be deleted (as they cannot be created either).
    */    
    ~MPIStats() {}
};

/** \def TRACK_IDLE_TIME

    \brief Convenience macro to track time spent in calling a MPI
    method (only if MPI is available).

    This macro provides a convenient wrapper around a given MPI method
    call to track the time spent in calling a MPI method.  This macro
    assumes that the MPI call being made must be accounted as idle
    time for the process and uses the Wtime() method to appropriately
    time the method call and adds the elapsed time to
    MPIStats::idleTime variable.  This macro must be used as shown
    below:

    \code

    #include "MPI.h"

    void someMethod() {
        MPI::Status statusInfo;
        TRACK_IDLE_TIME(MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, MPI_ANY_TAG,
                                              statusInfo));
    }

    \endcode

    \note If you are tracking idle time for the Probe method, use
    MPI_PROBE macro directly as it automatically tracks the idle time
    too.
*/
#ifdef HAVE_LIBMPI
#define TRACK_IDLE_TIME(MPIMethodCall) {                        \
        const double startTime = MPI::Wtime();                  \
        MPIMethodCall;                                          \
        MPIStats::idleTime += (MPI::Wtime() - startTime);       \
    }
#else
// If we don't have MPI this TRACK_IDLE_TIME reduces to nothing.
#define TRACK_IDLE_TIME(MPIMethodCall)
#endif

/** \def MPI_PROBE

    \brief Convenience macro to track time spent in calling a
    MPI::COMM_WORLD.Probe() method only if MPI is enabled. If MPI is
    not enabled, then this macro does nothing and call to MPI_Probe is
    never made.

    This macro provides a convenient wrapper around MPI's Probe method
    call to update MPIStats::probeCount variable and update the
    MPIStats::idleTime variable with the spent in the blocking Probe
    method.  This macro assumes that the Probe call being made must be
    accounted as idle time for the process and uses the Wtime() method
    to appropriately time the method call and adds the elapsed time to
    MPIStats::idleTime variable.  This macro must be used as shown
    below:

    \code

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
        MPI::Status statusInfo;
        MPI_PROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, statusInfo));
        // ... more code goes here ..
    }

    \endcode

    \note This macro call checks the statsInfo after the Probe call to
    ensure that the MPI call completed successfully.  If errors are
    found, then this method throws an exception.
*/
#ifdef HAVE_LIBMPI
#define MPI_PROBE(rank, tag, statusInfo)  {                \
        const double startTime = MPI::Wtime();             \
        MPI::COMM_WORLD.Probe(rank, tag, statusInfo);      \
        if ((statusInfo.Get_error() != MPI::SUCCESS) ||    \
            (statusInfo.Is_cancelled())) {                 \
            throw MPI::Exception(statusInfo.Get_error());  \
        }                                                  \
        MPIStats::idleTime += (MPI::Wtime() - startTime);  \
        MPIStats::probeCount++;                            \
    }
#else
// If we don't have MPI this MPI_PROBE reduces to nothing.
#define MPI_PROBE(rank, tag, statusInfo)
#endif

/** \def MPI_SEND

    \brief Convenience macro to track number of calls to
    MPI::COMM_WORLD.Send() method.

    This macro provides a convenient wrapper around MPI's Send method
    call to update MPIStats::sendCount variable.  This macro can be
    used as shown below:

    \code

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
        MPI_SEND(data, count, mpiDataType, rank, tag);
        // ... more code goes here ..
    }

    \endcode

    \note This macro calls MPI::COMM_WORD.Send only if MPI is
    enabled. If MPI is not enabled (or is unavailable) then this macro
    reduces to nothing and the actual send call is never made.
*/
#ifdef HAVE_LIBMPI
#define MPI_SEND(data, count, dataType, rank, tag)  {           \
        MPI::COMM_WORLD.Send(data, count, dataType, rank, tag); \
        MPIStats::sendCount++;                                  \
    }
#else
// If we don't have MPI then MPI_SEND reduces to nothing.
#define MPI_SEND(data, count, dataType, rank, tag)
#endif

/** \def MPI_RECV

    \brief Convenience macro to track number of calls to
    MPI::COMM_WORLD.Recv() method.

    This macro provides a convenient wrapper around MPI's Recv method
    call to update MPIStats::recvCount variable.  This macro can be
    used as shown below:

    \code

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
        MPI_RECV(data, count, mpiDataType, rank, tag);
        // ... more code goes here ..
    }

    \endcode

    \note This macro calls MPI::COMM_WORD.Send only if MPI is
    enabled. If MPI is not enabled (or is unavailable) then this macro
    reduces to nothing and the actual receive call is never made.
*/
#ifdef HAVE_LIBMPI
#define MPI_RECV(data, count, dataType, rank, tag)  {           \
        MPI::COMM_WORLD.Recv(data, count, dataType, rank, tag); \
        MPIStats::recvCount++;                                  \
    }
#else
#define MPI_RECV(data, count, dataType, rank, tag)
#endif

/** \def MPI_BCAST

    \brief A simple, convenience macro to track number of calls to
    MPI::COMM_WORLD.Bcast() method.

    This macro provides a convenient wrapper around MPI's Bcast method
    call to update MPIStats::bcastCount variable.  This macro can be
    used as shown below:

    \code

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
        MPI_BCAST(data, count, mpiDataType, senderRank);
        // ... more code goes here ..
    }

    \endcode

    \note This macro calls MPI::COMM_WORD.Bcast only if MPI is
    enabled. If MPI is not enabled (or is unavailable) then this macro
    reduces to nothing and the actual broadcast call is never made.
*/
#ifdef HAVE_LIBMPI
#define MPI_BCAST(data, count, dataType, senderRank)  {                 \
        MPI::COMM_WORLD.Bcast(data, count, dataType, senderRank);       \
        MPIStats::bcastCount++;                                         \
    }
#else
// If we don't have MPI then MPI_BCAST reduces to nothing.
#define MPI_BCAST(data, count, dataType, senderRank)
#endif

/** \def MPI_REDUCE

    \brief A simple, convenience macro to track number of calls to
    MPI::COMM_WORLD.Reduce() method.

    This macro provides a convenient wrapper around MPI's Reduce
    method call to update MPIStats::reduceCount variable.  This macro
    can be used as shown below:

    \code

    #include "MPI.h"

    void someMethod() {
        // ... some code goes here ..
        int localSuccess = 10, totalSuccess = 0;
        MPI_REDUCE(&localSuccess, &totalSuccess, 1, MPI::INT, MPI::SUM, 0);
        // ... more code goes here ..
    }

    \endcode

    \note This macro calls MPI::COMM_WORD.Reduce only if MPI is
    enabled. If MPI is not enabled (or is unavailable) then this macro
    reduces to nothing and the actual broadcast call is never made.
*/
#ifdef HAVE_LIBMPI
#define MPI_REDUCE(srcBufr, destBufr, count, dataType, op, destRank)  { \
   	    MPI::COMM_WORLD.Reduce(srcBufr, destBufr, count, dataType,		\
							   op, destRank);							\
        MPIStats::reduceCount++;                                        \
    }
#else
// If we don't have MPI then MPI_BCAST reduces to nothing.
#define MPI_REDUCE(srcBufr, destBufr, count, dataType, op, destRank)
#endif

#endif
