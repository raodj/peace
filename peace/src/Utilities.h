#ifndef UTILITIES_H
#define UTILITIES_H

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OH.
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

#ifdef _WINDOWS

/** A global not-equal-to function (not an member function).

    This function provides a plug-in replacement for the corresponding
	string comparison function present in Linux version of STL 
	implementation but is missing in Windows version. Note that this
	method is defined only under Windows. On Linux the default 
	operator!=() provided by basic_string is used instead.
*/
inline bool
operator!=(const std::string& s1, const std::string& s2) {
	return (s1.compare(s2) != 0);
}

/** A global insertion operator for string (not an member function).

    This function provides a plug-in replacement for the corresponding
	insertion operator for std::string present in Linux version of STL 
	implementation but is missing in Windows version. Note that this
	method is defined only under Windows. On Linux the default 
	operator<<() provided by std::basic_string is used instead.
*/
inline std::ostream&
operator<<(std::ostream& os, const std::string& str) {
	return (os << str.c_str());
}

#endif // _WINDOWS

/** \def UNREFERENCED_PARAMETER Workaround warning C4100 in VS 2005.

    This macro is a simple work around to supress the C4100 warning
	(unreferenced parameter).  This warning is generated at Level 4
	under visual studio 2005.  This warning has not been observed 
	under gcc (4.2) and therefore on Linux/gcc this macro resolves to
	nothing.
*/
#ifdef _WINDOWS

#ifndef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(param) param
#endif

#else
#define UNREFERENCED_PARAMETER(param)
#endif

/** \def fmax Macro to return the maximum of 2 values.

    This macro provides a replacement for the fmax() function
	defined in math.h under Linux.  However, VS 2005 does not
	have this method and this macro serves as a replacement.
*/
#ifdef _WINDOWS
#define fmax(x, y) ((x < y) ? y : x)
#endif

/** \def fmin Macro to return the minimum of 2 values.

    This macro provides a replacement for the fmin() function
	defined in math.h under Linux.  However, VS 2005 does not
	have this method and this macro serves as a replacement.
*/
#ifdef _WINDOWS
#define fmin(x, y) ((x < y) ? x : y)
#endif

/** \def vsnprintf_s Macro to define vsnprintf_s if not defined.

    This macro provides a replacement for the vsnprintf_s function
	defined in Windows but is absent in Unix/Linux. This macro 
    simply defines vsnprintf_s as vsnprinf in Unix and Linux.
*/
#if (!defined(_WINDOWS) && !defined(vsnprintf_s))
#define vsnprintf_s(buffer, bufSize, count, format, argptr) vsnprintf(buffer, count, format, argptr)
#endif

/** \def _fileno Macro to define _fileno if not defined.

    This macro provides a replacement for the \c _fileno function
	defined in Windows but is absent in Unix/Linux. This macro 
    simply defines \c _fileno as \c fileno in Unix and Linux.
*/
#if (!defined(_WINDOWS) && !defined(_fileno))
#define _fileno fileno
#endif

/** \def ctime_s Macro to define ctime_s if not defined.

    This macro provides a replacement for the \c ctime_s function
	defined in Windows but is absent in Unix/Linux. This macro 
    simply defines \c ctime_s as \c ctime_r in Unix and Linux.
*/
#if (!defined(_WINDOWS) && !defined(ctime_s))
#define ctime_s(buffer, size, time) ctime_r(time, buffer);
#endif

/** \def ASSERT(x)

    \brief Define a convenient macro for using c asserts.

    Define a custom macro ASSERT (note the all caps) method to be used
    instead of the c \c assert method to enable compiling ASSERT
    usages in and out of the source code.  When the compiler flag \c
    DEVELOPER_ASSERTIONS is specified then the ASSERT call defaults
    to the normal \c assert method.  For example ASSERT (\<<em>boolean
    expression</em>\>)} will be mapped to \c assert (\<<em>boolean
    expression</em>\>).  If the compiler flag \c DEVELOPER_ASSERTIONS
    is not specified then the ASSERT simply gets compiled out.
*/
#ifndef ASSERT
#ifdef DEVELOPER_ASSERTIONS
#include <assert.h>

#define ASSERT(x) assert(x)

#else // !DEVELOPER_ASSERTIONS

#define ASSERT(x)

#endif
#endif

/** \def VALIDATE(x)

    \brief Define a convenient macro for extra checks that can be
    easily compiled out to improve performance.

    Define a custom macro VALIDATE (note the all caps) to be used
    around additional extra checks that can be easily compiled out to
    improve performance. When the compiler flag \c ENABLE_VALIDATES is
    specified then the VALIDATE call defaults to the normal operation
    and are included in the source code.  If the compiler flag \c
    ENABLE_VALIDATES is not specified then the VALIDATE code simply
    gets compiled out.
*/
#ifndef VALIDATE
#ifdef ENABLE_VALIDATES

#define VALIDATE(x) x

#else // !ENABLE_VALIDATES

#define VALIDATE(x)

#endif
#endif


/** \def getTimeStamp

    \brief Get's the file modification timestamp for a given file
    name.

    This method provides a portable (to Windows and Linux/Unix)
    implementation for a helper method to obtain the modification
    timestamp for a given file.

    \param[in] fileName The file name (with full path) for which the
    modification time stamp is desired.  If the fileName is NULL then
    this method simply returns the buffer without any modifications.

    \param[out] buffer The buffer into which the time stamp is to be
    written.

    \return A pointer to the buffer.
*/
char* getTimeStamp(const char *fileName, char *buffer);

/** \def getTime

    \brief Returns the string representation of the supplied time data
    structure.

    This method provides a portable (to Windows and Linux/Unix)
    implementation for a helper method to obtain the string
    representation of a given encoded time_t
    datastructure.

    \param[out] buffer The buffer into which the string representation
    of the supplied time is to be written.  This pointer must be
    capable of holding at least 128 characters.  If this pointer is
    NULL, then this method exits immediately.

    \param[in] codedTime The encoded time_t data structure to be
    converted to a string representation.  If this parameter is NULL,
    then the current system time is converted to a string and filled
    into the supplied buffer.

    \return The pointer to the buffer passed in.
*/
char* getTime(char *buffer, const time_t *codedTime = NULL);

#endif
