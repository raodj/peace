#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H

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

/** Simple class to encapsulate a set of parameters for analyzers
    and/or heuristics.

    <p>The primary information associated with a ParameterSet object
    is the length of the cDNA fragment with which the object is
    associated.  The length of the cDNA fragment is represented by the
    minLength and maxLength instance variables.  If the length of a
    cDNA fragment falls within this range, then the remainder of the
    entries are used as parameters to various algorithms and
    heuristics.  Currently, this information is used by the adaptive
    TwoPassD2 analyzer.</p>
    
    <p>This class is created and managed via the ParameterSetManager
    class.  Consequently, the constructor for this class is
    private.</p>
*/
class ParameterSet {
    friend class ParameterSetManager;
public:
    /** The length of the shortest fragment which which this
        ParameterSet is to be associated.

        This value is set when then the object is instantiated and is
        never changed during the lifetime of this object.
    */
    const int minLength;

    /** The length of the longest fragment which which this
        ParameterSet is to be associated.

        This value is set when then the object is instantiated and is
        never changed during the lifetime of this object.
    */
    const int maxLength;

    /** The window (aka frame) size (in bp) to be used by the
        analyzer.
        
        This value is set when then the object is instantiated and is
        never changed during the lifetime of this object.  This value
        overrides the default values set in the analyzer.
    */
    const int frameSize;

    /** The shift (number of bp to skip) when running heuristics or
        the actual distance/similarity metric generating analyzer.
        
        This value is set when then the object is instantiated and is
        never changed during the lifetime of this object.
    */       
    const int frameShift;

    /** Threshold value in the analyzer that determines when it is
        safe to report that two fragments are sufficiently similar.

        The threshold value is used to permit the analyzer algorithm
        to return earlier (without checking all windows) when a window
        with sufficient similarty is found.  This value (may it be a
        distance or similarity metric) must be setup to be appropriate
        to the analyzer being used.
    */ 
    const int threshold;

    /** Threshold value in the analyzer that determines when it is
        safe to report that two fragments are sufficiently \e
        dissimilar.

        The threshold value is used to decide if the score from the
        faster asymmetric D2 algorithm is sufficiently large, thereby
        permitting bypassing the symmetric D2 phase.  This permits the
        TwoPassD2 analysis to return earlier (without checking all
        windows).  This value (may it be a distance or similarity
        metric) must be setup to be appropriate to the analyzer being
        used.
    */
    const int maxThreshold;

    /** The number of minumum number of common words, a TV Heuristic
        parameter (also uses frameSize).

        This instance variable contains the minimum number of words
        (that are close without having to be idential) that have
        matching values in pairs of cDNA fragments.  However, this
        value can be overridden by a command line argument.  This is
        used to adaptively set parameters for the <i>t/v</i> heuristic.
    */
    const int t;
    
    /** The threshold for number of common words, used by the
        <i>u/v</i> heuristic.

        This parameter contains the minimum number of common words
        that two cDNA fragments must contain in order to "pass" the
        <i>u/i</i> heuristic.
    */
    const int u;

    /** Number of base pairs to skip when comparing words in
        NewUVHeuristic.
        
		The number of base pairs that are to be skipped (in each pass)
		when building hashes and checking for common words.
	*/
    const int wordShift;
    
    /** The destructor.
        
        The destructor for the parameter set.  Currently, this class
        is just encapsulation class.  Consequently, the constructor
        does not have much work to do.
    */
    virtual ~ParameterSet() {}

protected:
    /** The constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead the class
        must be instantiated via the ParameterSetManager.

        \param[in] min  The length of the shortest fragment which which this
        ParameterSet is to be associated.

        \param[in] max The length of the longest fragment which which
        this ParameterSet is to be associated.

        \param[in] fsz The window (aka frame) size (in bp) to be used
        by the analyzer.

        \param[in] fsh The shift (number of bp to skip) when running
        heuristics or the actual distance/similarity metric generating
        analyzer.

        \param[in] thresh Threshold value in the analyzer that
        determines when it is safe to report that two fragments are
        sufficiently similar.

        \param[in] maxThresh Threshold value in the analyzer that
        determines when it is safe to report that two fragments are
        sufficiently \e dissimilar.

        \param[in] tIn The number of minumum number of common words, a
        TV Heuristic parameter (also uses frameSize).

        \param[in] uIn The threshold for number of common words, used
        by the <i>u/v</i> heuristic.

        \param[in] ws Number of base pairs to skip when comparing
        words in NewUVHeuristic.
    */
    ParameterSet(const int min, const int max, const int fsz, const int fsh,
                 const int thresh, const int maxThresh, const int tIn,
                 const int uIn,  const int ws) :
        minLength(min), maxLength(max), frameSize(fsz), frameShift(fsh),
        threshold(thresh), maxThreshold(maxThresh), t(tIn), u(uIn), 
        wordShift(ws) {
        // Nothing else to be done in the constructor for now.
    }
    
private:
    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.

        \return Reference to this.
    */
    ParameterSet& operator=(const ParameterSet& src);		
};

#endif
