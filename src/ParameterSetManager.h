#ifndef PARAMETER_SET_MANAGER_H
#define PARAMETER_SET_MANAGER_H

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

#include <vector>
#include "ParameterSet.h"

/** Class that manages a list of parameter sets.

    <p>This class represents a list of parameter sets that analyzers
    and heuristics can access to obtain appropriate parameters for
    analysis of sequences with varying characteristics, typically
    length.  The class initializes the appropriate parameter sets
    and provides the necessary parameter values when requested. </p>

*/
class ParameterSetManager {
public:
    /** Gets the appropriate parameter set for the given sequence lengths.

        \return The parameter set that should be used; NULL if the sequence
        lengths are too different to be compared.
    */
    ParameterSet* getParameterSet(const int seq1Len, const int seq2Len);

    /** Gets the maximum possible frame size over all parameter sets.

        \return the maximum possible frame size, an integer
    */
    int getMaxFrameSize();

    /** Get a pointer to the instance of the parameter set manager.
        
        Since this class is a singleton, the constructor is private
        and the only way to obtain an instance of the class is through
        this method.  The parameter set manager is available only after the
        \c setupParameters method (that is invoked from main right after
        command line arguments are validated) has successfully
        completed its operation.  Until such time this method simply
        returns \c NULL.

        \return The process-wide unique pointer to the parameter set manager.
    */
    static inline ParameterSetManager* getParameterSetManager() {
        return ptrInstance;
    }

    /** Set up the parameter set manager.

        Initializes the parameter set list to default values.
	The default is for 3 parameter sets, one for each of the following
	"classes" of sequences:
	Short  - 0 to 150 bases (but in practice it is 50 to 150, as the
	LengthFilter will filter out sequences of fewer than 50 bases length)
	Medium - 150 to 400 bases
	Long   - 400 to (no limit) bases
	
	In the future this class should incorporate command-line arguments
	and customizability for the parameter sets.
    */
    static void setupParameters();
    
    /** The destructor.

        The destructor frees up all the parameter sets added to this
        parameter set manager.
    */
    virtual ~ParameterSetManager();
    
protected:
    
private:
    /** The constructor.
	
        This is made private because the parameter set manager is
	a singleton, and should only be instantiated from the
	getParameterSetManager() static method.
    */
    ParameterSetManager();

    
    /** The vector containing the list of parameter sets.
	
    */
    std::vector<ParameterSet*> parameterSets;
    
    /** The pointer to the singleton instance of this class.

        Again, this is made private so that only methods of this class
        can access it. The getParameterSetManager() method in this class
        must be used to obtain an instance of this class.
    */
    static ParameterSetManager* ptrInstance;
};

#endif
