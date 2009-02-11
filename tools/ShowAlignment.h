#ifndef SHOW_ALIGNMENT_H
#define SHOW_ALIGNMENT_H

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
#include <vector>

// Forward declaration to keep compiler happy
class EST;

class ShowAlignment {
public:
    /** The main method to perform tasks for this class.

        This method is invoked from the ::main() method associated
        with PEACE tools once it has detected that the tool to be used
        is this ShowAlignment tool.  Any pending (unprocessed) command
        line parameters are passed to this method for its use.

        \param[in] argc The number of remaining command line arguments
        for use by this class/method.

        \param[in,out] argv The actual command line arguments for use
        by this class/method.
    */
    static int main(int argc, char *argv[]);
    
    /** Method to load data for processing from corresponding fasta files.

        This method must be used initially to load the necessary
        information from two different fasta files.

        \param[in] srcDataIndex The index of the source
        genome/transcript in the FASTA-compatible srcFileName.

        \param[in] srcFileName The FASTA-compatible data file that
        contains the source genes/transcripts from which ESTs were
        generated.  This method reads ESTs from this file until the
        necessary EST is read.

        \param[in] estFileName The name of the file that contains the
        actual sequence of ESTs to be aligned.

        \return This method returns \c true if all the data was read
        successfully.  On errors this method generates suitable error
        messages and returns \c false.
    */
    bool loadData(const int srcDataIndex, const char* srcFileName,
                  const char *estFileName);

    /** Dump the alignment to the output file.

	This method merely calls the drawEST() method to draw all the
	ESTs.  It first draws the refernece EST and then draws the
	remaining ESTs.

	\param[in] rectHeight The height of the rectangle (in XFig
	units) to be drawn for each EST.  If this value is -1, then
	the actual sequence is displayed.
    */
    bool drawAlignment(const int rectHeight);

    /** Obtain a row to draw an EST.

	This method may be used to obtain a row to draw an EST.  This
	method iteratively searches through the rowUsage vector to
	find a row which is not filled up to the start column. If a
	row is not found then a new row is added. The row usage is
	updated to the end column specified.

	\param[in] start The starting column from where the EST is to
	be drawn.

	\param[in] end The ending column where the EST ends.

	\return The logical row in the figure where the EST is to be
	drawn.
    */
    int getRow(const int start, const int end);
    
protected:
    /** Instance variable to hold the reference EST sequence.

        This instance variable is used to hold the reference sequence
        (typically a gene or transcript) using which the ESTs were
        generated.  This instance variable is pouplated in the
        loadData() method.
    */
    EST *refEST;

    /** Tracks the columns on each row used thus far.

	This vector tracks the maximum column (on each row) on which
	data has already been written.  This array is used to ensure
	that EST information do not overlap.
    */
    std::vector<int> rowUsage;
    
    /** The helper class to generate XFig information.

	This object is used to generate a XFig information to display
	the aligned information.
    */
    XFigHelper xfig;
    
    /** Helper method to load sequences from a FASTA file.

        This method is a helper method that is used to load data from
        a given FASTA file.

        \param[in] fileName The FASTA file name from which the FASTA
        sequences are to be loaded.
       
        \return This method returns \c true if all the data from the
        fasta file was successfully read.  On errors this method
        generates suitable error messages and returns \c false.
    */
    bool loadFastaFile(const char* fileName);

    /** Helper method to draw a given est.

	This helper method is used to draw a single EST as an XFig
	entity. This method operates as follows:

	<ol>

	<li>If the index is not -1, then it assumes the EST is not a
	reference EST and parses the header assuming it is in the
	format: g001_000001_000001. If the gene/transcript number (the
	first 3-digts) is not a match to the gene being processed then
	this method returns immediately. Otherwise it extracts the
	start and end positions for the EST (the second and third 6
	digit numbers respectively) for further use below.</li>

	<li>Next, it converts the index to a string to ease display.</li>

	<li>Next it determines the logical row in which the EST must
	be displayed by calling the getRow() method.</li>

	<li>First the index is rendered in Courier font using a fixed
	font size (currently 6 pt).</li>

	<li>If rectHeight is not -1, then this method draws a
	rectangle to fill the space for the EST. Otherwise it draws
	the actual sequence for the EST.</li>

	</ol>
     */
    void drawEST(const int index, const EST* est, const int rectHeight = -1);

private:
    /** The default constructor for this class.

        This class is meant to be instantiated only from the public
        static main() method.  Consequently, the constructor has been
        made private to force users to use the main() method instead.

	\param[in] os The output stream to which the Xfig information
	is to be written.
    */
    ShowAlignment(std::ostream &os);

    /** The destructor.

        The destructor frees any dynamic memory allocated to store the
        data members in this class.
    */
    ~ShowAlignment();

    /** Helper method to show usage information.

        This is a helper method that was introduced to reudce
        redundant code to display usage information.

        \param[in] ap The argument parser that is used to show some of
        the usage information.
    */
    static void showUsage(const arg_parser& ap);
};

#endif
