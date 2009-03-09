#ifndef DRAW_MATRIX_H
#define DRAW_MATRIX_H

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
#include <arg_parser.h>
#include "XFigHelper.h"

// Forward declaration to keep compiler happy
class EST;

class DrawMatrix {
public:
    /** The main method to perform tasks for this class.

        This method is invoked from the ::main() method associated
        with PEACE tools once it has detected that the tool to be used
        is this DrawMatrix tool.  Any pending (unprocessed) command
        line parameters are passed to this method for its use.

        \param[in] argc The number of remaining command line arguments
        for use by this class/method.

        \param[in,out] argv The actual command line arguments for use
        by this class/method.
    */
    static int main(int argc, char *argv[]);
    
    /** Dump the matrix information to the output file.

	This method reads line by line from the input matrix file and
	draws graphical squares for each entry for each row.

	\param[in,out] input The input stream from where the matrix
	data is to be read.
	
	\param[in] estCount The number of ESTs in the matrix.

	\param[in] sqSize The size of the square to be drawn for each
	entry in the matrix.
    */
    bool drawMatrix(std::istream& input, const int estCount,
                    const int sqSize);
    
protected:
    /** Obtain a row of matrix information to draw.

	This method is used to obtain a row of matrix information to
	draw. This method first reads a line from the file an then
	reads pairs of information (d2-distance EST index) from the
	file and stores the color codes corresponding to distance
	values into the supplied array of codes.  Unspecified entry
	values are set to WHITE.

	\param[in] estCount The number of ESTs in each row.

	\param[out] codes The color codes associated with each entry
	in the row.
	
	\return This method returns \c true if a line of data was
	read. Otherwise this method returns \c false.
    */
    bool getRow(const int estCount, std::istream& inFile,char *codes);

 private:
    /** The helper class to generate XFig information.
        
	This object is used to generate a XFig information to display
	the aligned information.
    */
    XFigHelper xfig;
    
    /** The default constructor for this class.

        This class is meant to be instantiated only from the public
        static main() method.  Consequently, the constructor has been
        made private to force users to use the main() method instead.

	\param[in] os The output stream to which the Xfig information
	is to be written.
    */
    DrawMatrix();
    
    /** The destructor.
        
        The destructor frees any dynamic memory allocated to store the
        data members in this class.
    */
    ~DrawMatrix();
    
    /** Helper method to show usage information.

        This is a helper method that was introduced to reudce
        redundant code to display usage information.

        \param[in] ap The argument parser that is used to show some of
        the usage information.
    */
    static void showUsage(const arg_parser& ap);
};

#endif
