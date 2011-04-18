#ifndef WINDOW_PAIR_H
#define WINDOW_PAIR_H

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

#include <iostream>

/** \file WindowPair.h

    \brief A simple class to encapsulate score and window index
    information about a pair of windows that exceed a given threshold.

    This file contains the declaration for the WindowPair class.

	@see BatonList::getWindowPairs()
	@see WindowPairList
*/

/** Encapsulate information about window pairs with significant number
    of identical batons
    
    This class encapsulates a 3-tuple that contains information about
    a pair of windows that exceed a given threshold.  Each window
    identifies a sub-fragment of a given cDNA fragment.  The score
    indicates the number of identical batons that were found between
    the given pair of windows These objects are created by
    BatonList::getWindowPairs() method.
*/
class WindowPair {
    friend std::ostream& operator<<(std::ostream&, const WindowPair&);
public:
    /** The one and only constructor.

        This is a convenience constructor used to build a WindowPair
        object.

        \param[in] win1 The index of the window in the first (or
        reference) cDNA fragment that is being analyzed to determine
        the number of identical batons.

        \param[in] win2 The index of the window in the second (or
        other) cDNA fragment that is being analyzed to determine the
        number of identical batons.

        \param[in] scoreVal The number of identical batons between the
        pair of windows.  Identical batons have the same \i n -mer
        heads and have the same length.

        @see BatonList::tallyBatons()
    */
    inline WindowPair(int win1, int win2, int scoreVal) :
        window1(win1), window2(win2), score(scoreVal) {}

    /** The destructor.

        The destructor does not have any special task to perform. It
        is here merely to adhere to coding conventions.
    */
    inline ~WindowPair() {}

    /** Default comparison operator.

        This method provides a custom implementation of operator<()
        (less-than operator) to compare two WindowPair objects.  This
        operator is handy in maintaining a sorted list of WindowPair
        objects. 

        \param[in] rhs The other WindowPair object with which this
        object is being compared.

        \return This method returns \c true if the score of the rhs
        WindowPair is greater than the score of this WindowPair.
    */
    inline bool operator<(const WindowPair& rhs) const {
        return (score < rhs.score);
    }

    /** Obtain the window in the first (or reference) cDNA.

        This method must be used to obtain the window in the first (or
        reference) cDNA for which this object contains the data.  The
        value returned by this method was the value set via the
        constructor (when the object was instantiated).

        \return The window in the first (or reference) cDNA.
    */
    inline int getWindow1() const { return window1; }

    /** Obtain the window in the second (or other) cDNA.

        This method must be used to obtain the window in the second
        (or other) cDNA for which this object contains the data.  The
        value returned by this method was the value set via the
        constructor (when the object was instantiated).

        \return The window in the second (or other) cDNA.
    */
    inline int getWindow2() const { return window2; }

    /** Obtain the score between the pair of windows associated with
        this object.

        This method returns the score between the pair of windows
        associated with this object.  This score indicates the number
        of identical batons in the given pair of windows.  The value
        returned by this method was the value set via the constructor
        (when the object was instantiated).

        \return The score between the pair of windows associated with
        this object.
    */
    inline int getScore() const { return score; }
    
private:
    /** The window in the first cDNA fragment.

        This instance variable maintains the index of the window in
        the first cDNA fragment.  This is the reference fragment used
        when computing identical batons via the BatonAnalyzer.
    */
    int window1;

    /** The window in the second cDNA fragment.

        This instance variable maintains the index of the window in
        the first cDNA fragment.  This is the reference fragment used
        when computing identical batons via the BatonAnalyzer.
    */    
    int window2;

    /** The number of identical batons in the given window pair.

        The score indicates the number of identical batons in the
        given window pair.  Identical batons have the name \i n -mer
        heads and same length.
    */
    int score;
};

#endif
