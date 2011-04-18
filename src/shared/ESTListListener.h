#ifndef EST_LIST_LISTENER_H
#define EST_LIST_LISTENER_H

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

/** Interface to intercept dynamic changes to ESTList.

    <p>This interface is used to dispatch updates to various
    components notifying them about changes in the list of cDNA
    fragments to be processed.  This interface was primarily
    introduced to provide a clean API for updating adaptive parameters
    (and other similar data structures) whenever EST entries are added
    or removed (typically by filters).</p>

    <p>Note that this class is meant to be a pure interface class and
    is not meant to be directly instantiated.  Instead, classes that
    are interested in receiving notifications about changes to the
    shared list of cDNA fragments must implement this interface and
    add themselves to the shared list of listeners in
    RuntimeContext::estListenerList.</p>

    @see ParameterSetManager
*/
class ESTListListener {
public:
    /** Notifications regarding entries added to the ESTList.

        This method is invoked (by one of the components in a
        SubSystem) when a set of EST entries have been added.

        \param[in] startIndex The zero-based starting index from where
        new EST entries were added.

        \param[in] endIndex The zero-based ending index (into the
        shared list of cDNA fragments pointed by
        RutimeContext::estList) where the last entry was added.

        \note It is guaranteed that startIndex <= endIndex.
    */
    virtual void entriesAdded   (int startIndex, int endIndex) = 0;

    /** Notifications regarding entries removed from the ESTList.

        This method is invoked (by one of the components in a
        SubSystem) when a set of EST entries have been removed.

        \param[in] startIndex The zero-based starting index from where
        new EST entries were removed. Note that this index value would
        not be really valid as entries have already been removed.

        \param[in] endIndex The zero-based ending index (into the
        shared list of cDNA fragments pointed by
        RutimeContext::estList) where the last entry was removed.

        \note It is guaranteed that startIndex <= endIndex.
    */    
    virtual void entriesRemoved (int startIndex, int endIndex) = 0;

    /** The destructor.

        The destructor has no specific operation to perform and is
        present merely to adhere to coding conventions.
    */
    virtual ~ESTListListener() {}
    
protected:
    /** A default dummy constructor.

        The constructor has been made protected to ensure that this
        class is never directly instantiated.
    */
    ESTListListener() {}
    
private:
};

#endif
