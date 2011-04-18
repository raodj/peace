#ifndef LOG_LEVEL_H
#define LOG_LEVEL_H

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

/** Class to define various log levels supported by the
    OutputSubSystem.

    This class provides a set of strongly-typed log level constant
    objects that can be used to set/change the log-level of a set of
    logs generated. Here is an example:

    \code

    #include "OutputSubSystem.h"

    void doSomething() {
        //....

        std::clog << LogLevel::INFO_LOG << "Number of ESTs loaded: "
                  << estCount << std::endl;
        // Continue an informational log
        std::clog << "More information...";

        // Change logs to debug logs...
        std::clog << LogLevel::DEBUG_LOG << "Starting clustering...\n";
        
        //...
    }

    \endcode
*/
class LogLevel {
    friend class StandardOutput;

    /** Enumerations identifying the different log levels.

        This enumeration defines the various log levels that are
        supported by the OutputSubSystem of PEACE.
     */
    enum LevelTypes {
        LOG_LEVEL_NONE, /**< Generate no log information.            */
        LOG_LEVEL_INFO, /**< Generate logs at informational level.   */
        LOG_LEVEL_DEBUG /**< Generate all logs written to std::clog. */
    };

public:
    /** Convenience operator  to compare two log levels.

        This method provides a convenience operator to compare two
        LogLevel objects.

        \param[in] rhs The other LogLevel object with which the
        current object is to be compared.

        \return This method returns true if this object has a
        logically greater log level than the rhs (parameter). The
        levels of logging are defined by the LogLevel::LevelTypes
        enumeration.
    */
    inline bool operator>=(const LogLevel& rhs) const {
        return level >= rhs.level;
    }

    /** Convenience operator to compare two log levels.

        This method provides a convenience operator to compare two
        LogLevel objects.

        \param[in] rhs A enumeration constant from the LevelTypes
        enumeration.

        \return This method returns true if this object has a
        logically greater log level than the rhs (parameter). The
        levels of logging are defined by the LogLevel::LevelTypes
        enumeration.
    */
    inline bool operator>=(const LevelTypes& rhs) const {
        return level >= rhs;
    }

    /** \var const LogLevel INFO_LOG
        
        \brief Constant to set/change the log level to informational log.
        
        This external global constant object must be used to identify
        log messages being generated at the informational level.
        Informational log messages are the default log levels.  These
        log messages provide the user with meaningful information.
    */
    static const LogLevel INFO_LOG;

    /** \var const LogLevel DEBUG_LOG
        
        \brief Constant to set/change the log level to debug logs.
        
        This external global constant object must be used to identify log
        messages being generated at a debug level.  Debug log messages are
        used to generate temporary messages for troubleshooting and
        debugging PEACE.  These log messages are primarily meant for the
        a programmer and may not have meaningful information to an user.
    */
    static const LogLevel DEBUG_LOG;

    /** \var const LogLevel NONE_LOG
        
        \brief Constant to identify the log level where no logs are
        generated.
        
        This external global constant object is a convenience constant
        that is used to identify the 'none' log level where all log
        messages are suppressed.  This log level must <b>not</b> be used
        to generate logs.  It is just a constant that is used for
        comparisons.
    */
    static const LogLevel NONE_LOG;
    
protected:
    /** The log level at which this log class is operating.
        
        This instance variable contain the log level that is defined
        by this log level class.  This value is set by the constructor
        and is never changed during the life time of this class.
    */
    const LevelTypes level;
    
private:
    /** The only constructor.

        The constructor for this class. The constructor has been made
        private to ensure only a small subset of objects of LogLevel
        type are defined.

        \param[in] logLevel The log level value to be set for this
        object.
    */
    LogLevel(const LevelTypes logLevel) : level(logLevel) {};

    /** The destructor.

        The destructor is private to ensure that objects of LogLevel
        type are never deleted (even by accident).
    */
    ~LogLevel() {}

    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    LogLevel& operator=(const LogLevel& src);

};

#endif
