#ifndef EST_ANALYZER_FACTORY_H
#define EST_ANALYZER_FACTORY_H

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

// For declaration to make compiler happy.
class ESTAnalyzer;

class ESTAnalyzerFactory {
public:
    /** Method to instantiate a suitable EST Analyzer.

        This method must be used to instantiate a suitable EST
        analyzer instance.  This method uses the name (parameter) to
        suitably instantiate an EST analyzer.  If the name is not
        valid, then this method returns NULL.

        \param[in] name The name of the EST analyzer to be
        instantiated.

        \param[in] refESTidx The reference EST index value to be used
        when instantiating an EST analyzer.  Index values start from 0
        (zero).  If this value is negative then this method returns
        NULL.

        \param[in] outputFileName The target file to which the
        analysis report is to be written (if any).  Note that this
        parameter may be ignored if this EST analyzer is used to
        generate clusters.  If the value of outputFileName is ""
        (empty string) then the outputs are streamed to standard out.
    */
    static ESTAnalyzer* create(const char* name, const int refESTidx,
                               const std::string& outputFileName);

    /** Method to display the list of EST analyzers available.

        This method is typically used to display the list of EST
        analyzers currently available.  This method is typically used
        in the main() method when displaying usage information.

        \param[out] os The output stream to which the list of EST
        names must be written.
    */
    static void displayList(std::ostream& os);
    
protected:
    // Currently this class has no protected members.
    
private:
    /** The default constructor.

        The default constructor has is private in order to ensure that
        this class is never instantiated.  Instead, the static methods
        in this class must be directly used to instantiate a EST
        analyzer.
    */
    ESTAnalyzerFactory() {}

    /** The destructor.

        The destructor is private to ensure that objects of this class
        are never deleted.
    */
    ~ESTAnalyzerFactory() {}
};


#endif
