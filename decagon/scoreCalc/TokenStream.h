#ifndef TOKEN_STREAM_H
#define TOKEN_STREAM_H

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

/** \file TokenStream.h

    \brief Class for reading an textual input stream as a series of
    tokens or words.

    This file provides the definition for the TokenStream class. This
    class provides a convenient mechanism to read data from an input
    stream (or file), one word (aka token) at a time.  The class also
    handles undo of tokens to enable predicting processing of files.
    Refer to the documentation for the various methodds in this class
    for further details.
*/

#include <sstream>

/** A general purpose class to read data from a text file as tokens or
	words.

	This class provides a generic mechanism to read data from a text
    file, one token (aka word) at a time. Each token is delimited by a
    white space (\c ' ', '\r', '\n', '\t'). Tokens are processed
    one-line at a time.  Blank lines result in an empty string being
    returned.

    \code
    throw EXP(1, "Out of memory", "Try increasing virtual memory size"); 
    \endcode
*/
class TokenStream {
public:
    /** Contructor to create token stream to generate tokens from a
        given input stream.

        This constructor must be used to create a token stream object
        that will generate tokens using the data from a given input
        stream. The data from the input stream is read, on demand, one
        line at a time and tokenized as needed.

        \note This class maintains an internal reference to the
        supplied input stream parameter.  Consequently, the stream
        should not be closed or read external to this class.

        \param[in] is The input stream from which tokens are to be
        generated on-demand.
    */
    TokenStream(std::istream& is);

    /** Obtain the next token to be processed.

        This method must be used to obtain the next token to be
        processed. This method returns the next token from either the
        source input stream or the tokens that were pushed back to be
        reprocessed via the undo() method in this class.

        \note Empty lines in the input source are returned as an empty
        string (\c "").
        
        \param[out] token The string which should contain the next
        token. If a next token is not available (that is the return
        value is \c false) then the contents of the string is
        undefined.

        \return This method returns \c true if the token parameter
        contains a valid token to be processed.  Otherwise this method
        returns \c false.
    */
    bool getToken(std::string& token);
    
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
};

#endif
