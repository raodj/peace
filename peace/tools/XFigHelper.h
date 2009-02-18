#ifndef XFIG_HELPER_H
#define XFIG_HELPER_H

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
#include <exception>
#include <string>

class XFigHelper {
public:
    /** The constructor.
            
        The constructor saves off the reference to a output stream for
        further use in this class.

        \param[out] os The output stream to which data has to be
        written.

	\param[in] genCustomColors If this flag is set to \c true then
	this method adds an user-defined color table to generated
	XFig. The user-defined color table guarantees a set of visibly
	different colors for use in applications.	
    */
    XFigHelper(std::ostream &os, const bool genCustomColors = false);

    /** The destructor.

        Currently the destructor does not have any specific tasks to
        perform.  It is present merely to adhere to coding convetions.
    */
    ~XFigHelper();

    /** Draw a line in the XFig code.

        This method can be used to generate the code necessary to draw
        a line.

        \param[in] x1 The x-coordinate of the first point.

        \param[in] y1 The y-coordinate of the first point.

        \param[in] x2 The x-coordinate of the second point.

        \param[in] y2 The y-coordinate of the second point.

        \param[in] colorCode An optional color code to specify color
        of the line.
    */
    void drawLine(const int x1, const int y1,
                  const int x2, const int y2,
                  const int colorCode = 0, const int level = 50);

    /** Dump code for displaying a text

        This method can be used to generate code for printing a text
        (or string) at a specific location.

        \return Return height of text in xfig units!
    */
    int drawText(const std::string& text,
                 const int x, const int y,
                 const int fontCode = 14,
                 const int fontSize = 12,
                 const int colorCode = 0,
                 const int level = 50);


    void drawRect(int x, int y, int width, int height, int colorCode = 0);
    
protected:
    /** Method to dump initial XFig header to the file.
        
        This is a helper method that must be used to dump a standard
        XFig header to the supplied output file.

	\param[in] genCustomColors If this flag is set to \c true then
	this method adds an user-defined color table to generated
	XFig. The user-defined color table guarantees a set of visibly
	different colors for use in applications.
    */
    void dumpHeader(const bool genCustomColors) const;

    /** Helper method generate color table.

        This method is invoked from the dumpHeader() method to
        actually generate the table of user-defined colors at the end
        of the header.
    */
    void generateColorTable() const;
    
private:
    /** The output stream to which data is to be written.

        This instance variable holds a reference to the output stream
        to which this class must dump all the XFig information.  The
        refernece is instialized in the constructor and used by
        various methods in this class.
    */
    std::ostream& os;
};

#endif
