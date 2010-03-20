#ifndef INTERACTIVE_CONSOLE_H
#define INTERACTIVE_CONSOLE_H

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

#ifndef _WINDOWS
#include "config.h"
#endif

#include "Utilities.h"
#include <vector>
#include <string>

// Make compiler happy and fast
class ESTAnalyzer;

/** An interactive console for manually studying ESTs.

    This class provides an text-based interactive console for studying
    and comparing ESTs loaded from a given FASTA file. This class
    provides a prompt for users to enter various commands.  This class
    processes the commands provided by the user and displays the
    results on the display. 
*/
class InteractiveConsole {
public:
    /** The constructor.

        The constructor is pretty straightforward and merely
        initializes all the instance variables to their default
        initial values.
    */
    InteractiveConsole(ESTAnalyzer* analyzer);

    /** The destructor.

        The destructor frees up memory allocated to instance
        variables.
    */
    ~InteractiveConsole();

    /** The main method of the console.

        This method is invoked from the main() function after the
        InteractiveConsole has been created.  This method performs the
        following tasks:

        <ol>

        <li>It first initializes the analyzer which causes the ESTs to
        be loaded from the supplied FASTA file.  If errors occur
        during initialization then this method reports an error and
        exits.</li>

        <li>It reads a command from the user and if the command is
        "exit" it exits from this method.</li>

        <li>Otherwise it tokenizes the input command line into
        multiple tokens (as a vector of strings) and then uses a
        dispatch table to invoke various methods to do the actual
        processing for a given command.</li>

        <li>
        
        </ol>
    */
    void processCommands();

protected:
    /** Helper method to tokenize a given string.

        This is a helper method that is used tokenize a given string
        based on a given set of delimiters.  This method uses a
        standard set of methods in std::string class to extract
        sub-strings from str and returns a std::vector containing the
        set of tokens.

        \param[in] str The string that must be broken to multiple
        tokens.

        \param[in] delims The set of characters to be used as
        delimiters. The default set of delimiters contains just the
        standard white space characters (namely: \c " \n\t\r").

        \return This method returns a std::vector containing the
        tokens in \c str broken using the set of \c delims specified.
    */
    static
    std::vector<std::string> tokenize(const std::string& str,
                                      const std::string& delims = " \n\t\r");
    
    /** A helper method to print statistics about ESTs.

        This method is invoked from the \c processCommands method
        (indirectly through a dispatch table).  This method computes
        and prints statistics about the list of ESTs loaded for
        analysis.

        \param[in] cmdWords The set of words/tokens from the command
        entered by the user.
    */
    void printStats(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords));

    /** Mehtod to windup the interactive console.

        This method is invoked from the \c processCommands method
        (indirectly through a dispatch table).  This method merely
        prints an exit message.
        
        \param[in] cmdWords The set of words/tokens from the command
        entered by the user. Currently, this list is ignored by this
        method.
    */
    void exit(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords));

    /** List the ESTs curently loaded.

        This method is invoked from the \c processCommands method
        (indirectly through a dispatch table).  This method prints the
        list of ESTs that are currently available for D2 analysis.
        
        \param[in] cmdWords The set of words/tokens from the command
        entered by the user. Currently, this list is ignored by this
        method.
    */
    void list(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords));

    /** Analyze a given pair of ESTs.

        This method is invoked from the \c processCommands method
        (indirectly through a dispatch table).  This method assumes
        that the user has specified the index of thw two ESTs to be
        analyzed and prints the result of analyzing the two ESTs.
        
        \param[in] cmdWords The set of words/tokens from the command
        entered by the user. This vector must have exactly 3 words,
        the first word being "analyze" and the other two words being
        either the index of ESTs or the fasta header for the ESTs.
    */
    void analyze(const std::vector<std::string>& cmdWords);

    /** Display a brief help regarding supported commands.

        This method is invoked whenever the user types the command
        "help" at the prompt. This method is invoked (from the \c
        processCommands method) via its method pointer stored in the
        \c cmdHandlerList. This method simply displays a static help
        text with the list of commands and some example usage.

        \param[in] cmdWords The set of words/tokens from the command
        entered by the user. Currently, this list is ignored by this
        method.
    */
    void help(const std::vector<std::string>& UNREFERENCED_PARAMETER(cmdWords));

    /** Display detailed information about a given EST

        This method is invoked whenever the user types the command
        "print" at the prompt. This method is invoked (from the \c
        processCommands method) via its method pointer stored in the
        \c cmdHandlerList. This method simply displays a static help
        text with the list of commands and some example usage.

        \param[in] cmdWords The set of words/tokens from the command
        entered by the user. This method expects exactly one parameter
        and process it assuming it is an EST index or an FASTA
        identifier corresponding to the EST whose information is to be
        displayed.
    */
    void print(const std::vector<std::string>& cmdWords);

#ifndef HAVE_LIBREADLINE
    /** A helper method used only under Windows.

	This is a replacement for the wonderful and interactive
	readline functionality that is available under Linux. This
	method displays the prompt and reads a line from the standard
	input from the user.

	\param[in] prompt The prompt string to be displayed to the user.

	\return A std::string containing the line of input entered
	by the user. Note that the returned pointer must be free'd
	by the caller.
    */
    static char* readline(const char *prompt);
#endif

private:    
    /** Helper method to convert index or FASTA identifier to index.

        This method is a helper method that is invoked from the \c
        analyze method to convert either an index value or a FASTA
        header to a index value.

        \param[in] id The identifier (either an number) or the FASTA
        identifier to be converted to an index value.

        \return The index value of the EST. If the index is invalid or
        the FASTA identifier was not found, this method prints a
        suitable message and returns -1.
    */
    int getESTIndex(const std::string& id) const;

    /** Helper method to perform initialization.

        This is a refactored helper method that was introduced to keep
        the code clutter in the processCommands method to a minimum.
        This method is invoked only once, right when the
        processCommands method is invoked from the \c main function.
        This method displays a standard startup message and then
        initializes the analyzer which causes the ESTs to be loaded
        from the supplied FASTA file.  This method also tracks and
        reports the time taken for loading the ESTs.

        \return This method returns \c true if the initialization was
        successful.  On errors it returns \c false.
    */
    bool initialize();
    
    /** The analyzer to be used for analysis.

        The analyzer to be used for performing the analysis to compute
        similarity or distance metrics when the user requests it. This
        pointer is set in the constructor and is never changed during
        the life time of this object.
    */
    ESTAnalyzer* const analyzer;

    /** \struct CmdEntry

        \brief A structure to ease delegating processing of commands
        to various methods in this class.
    */
    typedef struct Entry {
        /// The command with which this handler is associated.
        const char* cmd;
        /// The method to be invoked for processing a command.
        void (InteractiveConsole::*handler)(const std::vector<std::string>&);
    } CmdEntry;

    /** The list of command handlers associated in this class.

        This array contains the list of command handlers associated
        with this class. This list is a static list of handlers
        associated with various methods in this class.  This list
        enables invoking various methods from the \c processCommands
        method.
    */
    static CmdEntry cmdHandlerList[];

    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    InteractiveConsole& operator=(const InteractiveConsole& src);
};


#endif
