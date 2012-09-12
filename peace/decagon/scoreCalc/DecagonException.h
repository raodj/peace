#ifndef DECAGON_EXCEPTION_H
#define DECAGON_EXCEPTION_H

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

/** \file DecagonException.h

    \brief Class and macros for generating a DecagonException

    This file provides the definition for the DecagonException class
    that is used by the various methods and classes in a general
    purpose way to throw suitable exceptions when error conditions
    occur.  The DecagonException has been designed to include as much
    information about the exception as possible.  However, some of the
    information may not be useful/suitable for reporting exceptions to
    end users or should not be placed in release versions.
    Consequently, such information can be effectively compiled out.
    To ease uniform generation of exceptions this file also defines a
    set of EXP macros that must be used for throwing an exception.
    The working of the macros are controlled by compiler flags which
    automatically include or exclude additional debugging information.
    Refer to the documentation on the macros and the class for further
    details.
*/

#include <ostream>

/** A general purpose DecagonException class that may be readily reused.

    This class defines a generic exception that may be thrown by
    different classes and method.  The \c DecagonException class has
    been defined to encapsulate as much information as possible
    regarding an exception thrown.  The objective is to provide a
    detailed error message to aid in troubleshooting and debugging the
    source code.  The class provides various member functions that can
    be used to access the various pieces of information stored in the
    exception class, such as:

    <UL> <li> The error code associated with the exception (must be
	non-zero). </li>
        
	<li> The error string for the exception (maybe NULL). </li>

	<li> Any suggestions to recover from the error (maybe NULL). </li>

	<li> Source file name (NULL is not available). </li>

	<li> Line number in source file (zero if not available). </li>
    </UL>

    <p>The exceptions thrown by all methods constituting the
    implementation of a class may extend this class for further
    customization.  Note that the \c DecagonException class is fully
    self sufficient and can be directly used for generating
    exceptions.  However, if need for further customization is
    desired, then a new exception class can be created with the \c
    DecagonException class as its parent. </p>

    <p>Depending on the type of environment being developed (\em i.e.,
    for release or internal debugging) the DecagonException can
    include additional information (such as source code file name and
    line number) to aid in debugging.  To ease generation of such
    exceptions, a EXP macro has been defined, that suitably creates a
    new \c DecagonException class using suitable parameters.</p>

    \note When throwing a \c DecagonException always use the \c EXP
    macro defined in this header file.  Here is a short example of how
    to use the \c BGE macro:
    
    \code
    throw EXP(1, "Out of memory", "Try increasing virtual memory size"); 
    \endcode
*/
class DecagonException {
    friend std::ostream& operator<<(std::ostream&, const DecagonException&);
public:
    /** The default and only constructor.

        <p>The \c DecagonException class has been defined to
        encapsulate as much information as possible regarding an
        exception thrown.  The objective is to provide a detailed
        error message to aid in troubleshooting and debugging the
        source code.</p>

        <p>All the information regarding an exception must be provided
        as a part of the constructor.  The primary objective is to
        eliminate any unwanted side effects that may arise when
        creating the exception (and the members of the exception are
        all constant objects) and to avoid uncontrolled changes to an
        exception once it has been created.</p>

        \param errorCode The error code set for this exception.

        \param errorMessage The error message string.

        \param suggestion A text message suggesting how to avoid this
        exception.

        \param fileName An \em optional source code fileName.

        \param lineNumber An \em optional source code line number.
    */
    DecagonException(const int errorCode, const char* errorMessage,
                     const char* suggestion = "",
                     const char* fileName   = "",
                     const int   lineNumber = 0);
    
    /** Obtain the error code for this exception.

        This method must be used to obtain the error code associated
        with the exception.  The error code is a non-zero integer
        value that indicates the error number associated with the
        exception.

        \return The error code associated with the exception.
    */
    int getErrorCode() const;

    /** Obtain the error message for this exception.

        This method must be used to obtain the error message
        associated with this exception.  The error message is
        represented as C string ('\0' terminted array of characters)
        containing an ASCII text description of the error.

        \return The error message associated with the exception (never
        NULL).
    */
    const char* getErrorMessage() const;

    /** Obtain the help or suggestion text for this exception.

        This method must be used to access any help or suggestion
        messages associated with the exception.  The error message is
        represented as C string ('\0' terminted array of characters)
        containing an ASCII text that provides tips or suggestions on
        how to avoid the exception being thrown.  Note that the
        suggestion is optional and not all exceptions may include a
        suggestion.

        \return The suggestion associated with the exception (never NULL).
    */
    const char* getSuggestion() const;

    /** Source file name for the exception (only Debug version).

        This method maybe used to determine the source code file name
        from where this exception was thrown.  The source file name is
        filled in only in the debug versions.  In the non-debug
        versions, the value returned by this method is an empty string
        (\i i.e., "\0") but not NULL.

        \return An optional source file name for the exception.
    */
    const char *getFileName() const;

    /** Obtain the source code line number (only Debug version).

        This method maybe used to determine the source code line
        number from where this exception was thrown.  The source code
        line number is filled in only for the dbeug version.  In the
        non-debug versions, this method returns 0 (zero).

        \return The optional source file line number.
    */
	int getLineNumber() const;
    
    /** The destructor.

        The DecagonException class does not allocate any memory or
        perform signficiant operations to avoid auxiliary exceptions
        from being generated.  Consequently, the destructor has no
        tasks to perform.  It is here just to adhere to the standard
        convention of having a constructor-destructor pair for every
        class.
    */
    ~DecagonException();
    
protected:
    // Currently this class has no protected members.
    
private:
    /// DecagonException error code
    /** This member contains the error code set for this exception.
        The error code cannot be 0 (zero).
    */
    const int errorCode;

    /// DecagonException error message
    /** This member contains a detailed string description for the
        given exception.  Some exceptions may not have an error
        message associated with them.  However, it is \em strongly
        recommended that every exception have a error message
        associated with it,

        \note A message can never be longer than 256 bytes!  This has
        been adopted to ensure that creation of an exception does not
        throw further exceptions (which could occur if memory for the
        message strings are dynamically allocated).
    */
    char message[256];

    /// A suggestion to avoid or recover from exception
    /** The objective of this member is to provide a brief text (ASCII
        string) message providing some feedback or suggestion on how
        to avoid (future exceptions) or recover from the exception.
        This is an optional piece of information that maybe omitted.
    */
    char suggestion[256];
    
    /// The source file which generated exception.
    /** This member contains the source code file name from which this
        exception was generated.  Note that this information is
        optional and may be omitted in certain versions (such as
        release versions).
    */
    char fileName[1024];
    
    /// The source code line number
    /** This member contains the line number in the sourceCode
        (specified by the fileName member) from which the exception
        was generated.  If no line number information was provided (as
        in the case of release binaries) this member is set to 0
        (zero).
    */
    const int lineNumber;

};

/** \def EXP(code, errorMessage, suggestion)

    \brief Macro to ease throwing a DecagonException.

    This macro provides a convenient mechanism to create a DecagonException,
    depending on wether the DEVELOPER_ASSERTIONS flag has been turned
    on or off.  In case the DEVELOPER_ASSERTIONS flag has been
    sepcified, then this macro automatically tags the source code file
    name and source code line number from where the exception was
    raised.  This makes debugging much faster.
*/

#ifdef DEVELOPER_ASSERTIONS

#define EXP(code, errorMessage, suggestion) DecagonException(code, errorMessage, suggestion, __FILE__, __LINE__)

#else

#define EXP(code, errorMessage, suggestion) DecagonException(code, errorMessage, suggestion)

#endif

/** \fn std::ostream& operator<<(std::ostream&, const DecagonException&)

    Insertion operator to stream exception information to a given
    output stream.  This method provides a convenient mechanism to
    dump the complete exception information for debugging purposes.

	\param[out] os The output stream to which the formatted exception
	information is to be written.

	\param[in] de The decagon exception to be written.

	\return This method returns os as per the API requirements.
*/
extern std::ostream& operator<<(std::ostream& os, const DecagonException& de);

#endif
