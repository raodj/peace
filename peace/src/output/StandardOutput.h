#ifndef STANDARD_OUTPUT_H
#define STANDARD_OUTPUT_H

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

#include "Component.h"
#include "LogLevel.h"
#include <ostream>

/** An output sub-system component to redirect standard output streams.

    This class is a component of the output sub-system. Its
    responsibilities include:

    <ul>

    <li>Redirection of the standard output stream std::cout to a file
	if the user has specified \c --stdcout command line argument.  The
	destination file name is obtained from the command line argument
	and suffixed with the MPI rank to make the file name unique for
	each MPI process.</li>

    <li>Redirection of the standard error stream std::clog to a file
	if the user has specified \c --stdlog command line argument.  The
	destination file name is obtained from the command line argument
	and suffixed with the MPI rank to make the file name unique for
	each MPI process.</li>

    <li>Setup the log level at which logs are to be generated to the
	standard log file.  This task is performed based on the \c
	--logLevel command line argument.  Valid log levels in increasing
	levels of verbosity are: \c none, \c info, and \c debug. The
	default log level is \c info.</li>
    
    </ul>
*/
class StandardOutput : public Component {
    friend class OutputSubSystem;
    friend std::ostream& operator<<(std::ostream&, const LogLevel&);
public:
    /** The destructor.

        The destructor is merely present to adhere to coding
        conventions.  Currently, this does not have any special tasks
        to perform.
    */
    virtual ~StandardOutput();

    /** Add the set of command line parameters for this component.

        This method is typically invoked by the OutputSubSystem when
        it is being requested for command line arguments.  This method
        adds the various command line arguments that can be used to
        further customize the operations of this component.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    void addCommandLineArguments(ArgParser& argParser);

    /** Initialize the operations of this component.

        This method is invoked when the OutputSubSystem is being
        initialized.  This method (along with helper methods) performs
        the necessary operations to redirect standard output streams
        as specified by the user.

        \return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    bool initialize();

    /** Wind-up the operations of this component.

        This method is invoked when the OutputSubSystem is being
        finalized.  This method closes any redirected output streams
        that was opened in the initialize() method.
    */
    void finalize();
    
protected: 
    /** The constructor.

        This is the only constructor for this class. The
        constructor(s) are not public to ensure that this class is not
        instantiated directly but only by the friends of this
        component.  Currently, only the OutputSubSystem class can
        instantiate this component.  Currently the constructor does
        not have any special tasks to perform and merely initializes
        instance variables to their default initial value.
    */
    StandardOutput();

    /** Redirect a given stream to a file.

        This is a helper method that is called from the initialize()
        method to redirect a given output stream to a given file.

        \param[in] os The output stream to be redirected. Typically
        this is std::cout or std::clog.

        \param[out] currStreamBuf A pointer to the current stream
        buffer that is being used by this stream.
        
        \param[in] fileName The target file name to which the stream
        is to be redirected.  If running under MPI, then this file
        name is prefixed with the rank of the process on which this
        method is being invoked.

        \return On successful redirection, this method returns a
        pointer to the stream to which the data is being
        redirected. On errors this method returns NULL.
    */
    std::ostream* redirect(std::ostream& os, std::streambuf* &currStreamBuf,
                           const std::string& fileName);
    
private:
    /** The file name to which the standard output stream is to be
        redirected.

        This instance variable is used to track the file name
        specified by the user (if any) along with the \c --stdout
        command line argument.
    */
    std::string coutFileName;

    /** The output stream to which the standard output is being
        redirected.

        This instance variable is used to maintain a pointer to the
        output stream to which the standard output is being
        redirected. This pointer is set in the initialize() method if
        the user has specified a valid file via the \c --stdout
        command-line argument.
    */
    std::ostream *coutStream;

    /** Pointer to the default/system stream buffer used by std::cout.
            
        This instance variable is used to hold a pointer to the
        default/system-provided stream buffer that is used by
        std::cout to generate log messages.  This pointer is setup by
        the initialize() method.
    */
    std::streambuf *coutSystemStreamBuf;
    
    /** The file name to which the standard log stream is to be
        redirected.

        This instance variable is used to track the file name
        specified by the user (if any) along with the \c --stdlog
        command line argument.
    */
    std::string clogFileName;

    /** Pointer to the default/system stream buffer used by std::clog.
            
        This instance variable is used to hold a pointer to the
        default/system-provided stream buffer that is used by
        std::clog to generate log messages.  This pointer is setup by
        the initialize() method.
    */
    std::streambuf *clogSystemStreamBuf;

    /** The output stream to which the standard log is being
        redirected.

        This instance variable is used to maintain a pointer to the
        output stream to which the standard log is being
        redirected. This pointer is set in the initialize() method if
        the user has specified a valid file via the \c --stdlog
        command-line argument.
    */
    std::ostream *clogStream;

    /** Pointer to a wrapper stream buffer that is used to turn
        logging on/off.

        This instance variable is used to hold a pointer to a custom
        stream buffer to enable/disable logging of any data written to
        the log stream.  This pointer is setup by the initialize()
        method and deleted by the finalize() method.
    */
    std::streambuf *clogWrapperStreamBuf;
    
    /** The log level that the user has selected.
        
        This instance variable tracks the preferred level of logging
        that the user has selected.  This value is set via the \c
        --logLevel command-line argument.
    */
    std::string logLevelStr;

    /** The current log level set by the user.

        This instance variable tracks the current log level that has
        been set by the user via a command line argument.  This value
        is set by the initialize() method. This value is used to
        switch between #clogOldStreamBuf and #clogDummyStreamBuf
        whenever some code fragment in PEACE writes a log level
        constant object to an output stream as shown in the example
        below:

        \code
        
        std::clog << LogLevel::INFO_LOG << "Number of ESTs loaded: "
                  << estCount << std::endl;

        \endcode

        The actual change in object occurs in the operator<<() method
        implemented by this class.
    */
    static LogLevel::LevelTypes logLevel;
};

/** \fn std::ostream& operator<<(std::ostream&, const LogLevel&)

    Global method to handle log redirection when a LogLevel object is
    written to an output stream.  This method switches between the
    default and dummy (that simply discards data) stream buffers so
    that logs are generated or discarded depending on the log level
    set by the user.

    \param[in] os The output stream to which the log level is being
    written.

    \param[in] ll The LogLevel object that is being used to flag the
    log level information of the data that is going to be written.

    \return This method returns os back.
*/
std::ostream& operator<<(std::ostream& os, const LogLevel& ll);

#endif
