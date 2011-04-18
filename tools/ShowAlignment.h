#ifndef SHOW_ALIGNMENT_H
#define SHOW_ALIGNMENT_H

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

#include "Tool.h"

#include <stdio.h>
#include <vector>

class ShowAlignment : public Tool {
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
        generated.  This method reads ESTs from this file and extracts
        the cDNA fragment at the given srcDataIndex.  If this file
        name is an empty string (\c "") then this file is ignored.

        \param[in] estFileName The name of the file that contains the
        actual sequence of ESTs to be aligned.

        \return This method returns \c true if all the data was read
        successfully.  On errors this method generates suitable error
        messages and returns \c false.
    */
    bool loadData(const int srcDataIndex, const std::string& srcFileName,
                  const std::string& estFileName);

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

    /** \typedef std::pair<int, int> UseEntry

        \brief A typedef for a pair of values to store starting and
        ending column usages in a given row.

        This typedef provides a shortcut to refer to entries
        containing the space usage in each row.  \c UsageEntry.first
        indicates the smallest (or first) column (on a given row)
        where EST information has already been written while \c
        UsageEntry.second contains the last column where a EST entry
        has been written.
    */
    typedef std::pair<int, int> UseEntry;
    
    /** Tracks the columns on each row used thus far.

        This vector tracks the maximum column (on each row) on which
        data has already been written.  This array is used to ensure
        that EST information do not overlap.
    */
    std::vector<UseEntry> rowUsage;
          
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
    */
    ShowAlignment();

    /** The destructor.

        The destructor frees any dynamic memory allocated to store the
        data members in this class.
    */
    ~ShowAlignment();

    /** Determine if this EST of interest for showing alignment.

        This method is used to determine if an EST (given its index)
        is of interest for displaying in the alignment figure. This
        method takes into account any reference sequences that may
        have been sepcified and compares its gene squence with the
        gene information (assuming that the EST sequence is in the
        form \c g001_ or \c g002 etc.).

        \param[in] estIndex The index of EST that must be tested to
        verify if it is an EST of interest to be displayed in the
        alignment information.
    */
    bool isOfInterest(const int estIndex) const;

    /** Obtain start and end column from EST's FASTA header.

        This is a helper method that is used to obtain the alignment
        information from an EST's FASTA header.  This method
        consistently trys and extracts information from two different
        formats of ESTs: "g001_<stcol>_<endCol>" or
        "...|<genStCol><genEndCol>".

        \param[in] est The EST whose starting and ending column are to
        be extracted.

        \param[out] startCol This parameter is set with the starting
        column where the EST is to be aligned.

        \param[out] endCol This parameter is set to the ending column
        where the EST ends.

        \return This method returns \c true if the data was
        successfully read.  If data was not read successfully, then
        this method generates a suitable error message and returns \c
        false.
    */
    bool getStartEnd(const EST* est, int &startCol, int& endCol) const;

    /** Obtain color coding information (if any) for a given EST.

        This method is a helper method that is used to determine
        additional color coding information for a given EST.  If a
        color code is available (based on clustering) then this method
        returns that color. Otherwise this method returns black as the
        default color.

        \note This method consistently trys and extracts original gene
        name from two different formats of ESTs:
        "g001_<stcol>_<endCol>" or "...|<genStCol><genEndCol>" to
        provide consistent color coding based on clustering
        information.

        \param[in] est The EST for which color coding information is
        to be computed.

        \return This method returns a color code extracted from the
        colorMap.  If a color code is not found, then this method
        returns black as the default color code.
    */
    int getColor(const EST* est) const;
};

#endif
