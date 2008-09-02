#ifndef RESULT_LOG_H
#define RESULT_LOG_H

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

#include <stdio.h>
#include <string>

class ResultLog {
public:    
    /** Constructor.

        The constructor establishes a ResultLog.  The result outputs
        are written to standard output or a file in plain text or HTML
        format.

        \param[in] fileName The output/destination file to which logs
        are to be cut.  If this parameter is NULL then data is written
        to standard output.

        \param[in] html If this parameter is \c true then logs are
        written in HTML format so that they can be easily viewed in
        via a browser.
    */
    ResultLog(const std::string& fileName = std::string(""),
	      const bool html = false);

    /** Destructor.

        Winds up the process of generating a log.  If HTML logs are
        being generated then it generates trailers to wind up HTML
        formats.  If a file was opened then the destructor closes the
        destination file.
    */
    ~ResultLog();

    /** Dumps a line of result to the log.

        This method provides a convenient mechanism to dump a line of
        output to the specified destination.  If no destination is
        specified it dumps the line to standard output.  This method
        also handles the task of including necessary html tags if HTML
        output is requested.
        
        \param[in] line The line of text to be displayed.  This line
        can take the standard format as the printf() method thereby
        facilitating dumping of complex output in a simple manner.
    */
    void reportLine(const char *line, ...);

    /** Dumps a line of result using a 3-column table format.

        This method provides a convenient mechanism to dump a line of
        output in 3-column format to the specified destination.  If no
        destination is specified it dumps the line to standard output.
        This method also handles the task of including necessary html
        tags if HTML output is requested.
        
        \param[in] col1 The data to be displayed in the first column.
        This text can take the standard format as the printf() method
        thereby facilitating dumping of complex output in a simple
        manner.

        \param[in] col2 The data to be displayed in the second
        column. This text can take the standard format as the printf()
        method thereby facilitating dumping of complex output in a
        simple manner.

        \param[in] col3 The data to be displayed in the third
        column. This text can take the standard format as the printf()
        method thereby facilitating dumping of complex output in a
        simple manner.

        \note This method accepts a variable number of arguments
        corresponding to the printf() style format specifiers in the
        three text entries.  The order of arguments must match the
        order of format specifiers in col1, col2, and col3
        respectively.  Here is an example:

        \code
        
        ResultLog log("test.html", true);
        log.report("%s", "%c %d", "%ld", "test", '*', 123, 3.142);

        \endcode
    */
    void report(const char* col1, const char*col2, const char *col3, ...);

    /** Starts an HTML table if necessary.

        This method generates HTML code to start a table.  This method
        peforms all the checks for handling HTML mode and generates
        headers only if it is needed.  Consequently, it is safe to
        call this method repeatedly without much performance penalty.

        \param[in] titles The titles for the rows to be used when
        creating the table.  If this parameter is NULL, then no titles
        are included for the table.
        
    */
    void startTable(const char* titles[] = NULL);

    /** Ends an HTML table if necessary.

        This method generates HTML code to end a HTML table.  This
        method peforms all the checks for handling HTML mode and
        generates headers only if it is needed.  Consequently, it is
        safe to call this method repeatedly without much performance
        penalty.
    */
    void endTable();
    
protected:
    /** Helper method to dump HTML headers.

        This method is a helper method that is invoked from the
        constructor to dump HTML headers.
    */
    void dumpHTMLHeaders();  
    
private:
    /** The file to which results are to be logged.

        This member object holds the file to which results are to be
        logged.  This pointer is initialized in the constructor to
        either point to a result file or standard output.
    */
    FILE *output;

    /** Flag to indicate if HTML format is being used.

        This instance variable is used to determine if the logs are to
        be generated in HTML format or not.  This value is set in the
        constructor and is never changed during the life time of this
        class.
    */
    const bool html;

    /** Flag to indicate if a HTML table is currently open.

        This flag is initialized to false.  It is set to true each
        time a HTML table is started.  It is reset to false when the
        table is ended.  This flag is used in the polymorphic report()
        methods in this class to appropriately handle HTML tables.
    */
    bool insideTable;

    /** A dummy operator=

        The operator=() is supressed for this class as it has constant members
        whose value is set when the object is created.  These values cannot be
        changed during the lifetime of this object.

        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.

        \return Reference to this.
    */
    ResultLog& operator=(const ResultLog& src);

};

#endif
