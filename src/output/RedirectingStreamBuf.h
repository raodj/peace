#ifndef REDIRECTING_STREAM_BUF_H
#define REDIRECTING_STREAM_BUF_H

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

#include <streambuf>

/** A custom stream buffer that redirects output to another stream buffer.

    This class is a custom stream buffer that redirects the data being
    written to it to another stream buffer. The redirection of data
    can be enabled and disabled dynmaically.  This class is very
    straightforward.
*/
class RedirectingStreamBuf : public std::streambuf {
public:
    /** The constructor.

        The constructor sets up the target stream buffer to which data
        is to be redirected. In addition, it enables redirection (by
        default).

        \param[in] destination The target stream buffer to which data
        is to be redirected. This pointer cannot be NULL.
    */
    RedirectingStreamBuf(std::streambuf* destination);

    /** The destructor.

        Currently, the destructor does not have any specific tasks to
        perform and is present here to adhere to coding conventions.
    */
    ~RedirectingStreamBuf() {}

    /** Enable/disable redirection of output.

        This method can be used to enable/disable redirection.

        \param[in] flag If the flag is \c true then redirection is
        enabled. If the flag is \c false redirection is disabled.
    */
    void redirect(bool flag) {
        redirectFlag = flag;
    }

    /** Determine if redirection is enabled/disabled.

        \return This method returns \c true if redirection is enabled
        and data is being written to the original stream buffer.
    */
    inline bool isRedirecting() const { return redirectFlag; }

    /** Override the default base-class method to redirect data.

        This method override the default implementation in the base
        class to appropriately redirect (or suppress) the data to the
        original stream buffer.

        \param[in] data The data to be written.

        \param[in] n The number of bytes to be written.

        \return The number of bytes actually return. If the data is
        suppressed this method always returns n pretending all the
        data was successfully written.
    */
    virtual std::streamsize sputn (const char *data, std::streamsize n) {
        return (redirectFlag ? target->sputn(data, n) : n);
    }

    /** Override the default base-class method to redirect data.

        This method override the default implementation in the base
        class to appropriately redirect (or suppress) the data to the
        original stream buffer.

        \param[in] c The character to be written.

        \return This method returns traits_type::eof() if data was
        supressed. Otherwise it returns the value returned by calling
        the corresponding method on the original string buffer.
    */
    virtual int sputc ( char c ) {
        return (redirectFlag ? target->sputc(c) : traits_type::eof());
    }

protected:
    /** A dummy operator=
            
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    RedirectingStreamBuf& operator=(const RedirectingStreamBuf& src);

private: 
    /** The stream buffer to which data is being written.

        This value is initialized in the constructor to point to the
        stream buffer to which data is to be redirected.
    */
    std::streambuf* const target;

    /** Flag to indicate if redirection is enabled/disabled.

        If this flag is set to \c true, then data is redirected to the
        #target stream buffer.
    */
    bool redirectFlag;
};


#endif
