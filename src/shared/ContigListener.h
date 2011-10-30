#ifndef CONTIG_LISTENER_H
#define CONTIG_LISTENER_H

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

// Forward declaration to keep compiler happy.
class Contig;
class Assembler;

/** \file ContigListener.h

	\brief A event-listener interferace class to notify interested
	listeners that a genome-assembler has finished forming a contig.

	<p>This file contains the definition for a pure-interface class
	that provides a convenient mechanism to dispatch notifcations
	whenever a contig has been formed. Interested listeners (classes
	that implement the pure-virtual method in this class) can perform
	any necessary operations needed for their operations.</p>

    <p>This file also defines a convenience vector to operate on a
    list of contig listeners. Refer to the documentation on the class
    and interface method for additional details.</p>
*/

/** A simple interface class for dispatching notifications once a
    contig has been formed.

    This class provides a convenient mechanism to notify various
    independent sub-systems constituting PEACE about formation of
    contigs.  This class contains only a single pure-virtual method
    that interested classes can implement and register themselves with
    the Assembler (see Assembler::addContigListener() method).
    Whenever a contig formulation has been completed, the Assembler
    notifies (calls the pure-virtual method) the listeners so they can
    process the contig suitably.  For example, see: SAMWriter.
*/
class ContigListener {
public:
    /** Interface method to dispatch notification about successful
        contig formation.

        This method is called on all parallel processes once a contig
        has been successfully formed. The newly formed contig is
        passed in as the parameter. Note that not all
        genomic-assemblers collate all the information on all parallel
        processes.  Consequently, it must not be assumed that the
        contig contains full information.  Instead the contig may
        contain only partial information that is available on the
        parallel-process on which this method is being invoked.

        \param[in] assembler The genomic-assembler that has formulated
        the contig.
        
        \param[in] contig The newly formed contig with partial
        information that is available on the process on which this
        method is invoked.

        \param[in] fullContig If this flag is \c true then the contig
        has full information.  If this flag is \c false then this
        contig has partial information and data from other
        parallel-processes must be fused together to obtain full
        contig information.
        
        \return This method is expected to return \c true if the
        contig has been processed successfully.  On errors this method
        must return \c false.
    */
    virtual bool contigFormed(const Assembler& assembler,
                              const Contig& contig, const bool fullContig) = 0;
    
protected:
    /** Default constructor.

        This is a simple default constructor there is present more to
        adhere to standards.  The constructor is protected to ensure
        that this class is never directly instantiated (as it is not
        meant to be instantiated directly). Instead a derived class
        that implements the pure-virtual method in this class must be
        instantiated and used.
    */
    ContigListener();

    /** Destructor.

        The destructor has no specific operations to perform. It is
        here merely to adhere to coding conventions.
     */
    virtual ~ContigListener();

private:
    // This class does not have any private members
};

/** \typedef std::vector<ContigListener*> ContigListenerList

    \brief A vector that contains a list of ContigListener objects.

    This typedef provides a convenient short cut to refer to a vector
    that is used to hold a list of ContigListener objects.  This
    vector is used in more than one spot in the code. Consequently, it
    has been defined here to serve as a common data structure.
*/
typedef std::vector<ContigListener*> ContigListenerList;

#endif
