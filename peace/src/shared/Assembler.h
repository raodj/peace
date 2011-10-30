#ifndef ASSEMBLER_H
#define ASSEMBLER_H

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

#include "Component.h"
#include "ArgParser.h"
#include "Contig.h"
#include "ContigListener.h"

#include <fstream>

// Some forward declarations to keep compiler fast and happy
class ESTAnalyzer;
class ESTList;

/** \file Assembler.h
    
    \brief The common base class and primary Application Program
    Interface (API) for all gene assemblers.

    This file contains the class declaration for the Assembler class.
    This class is the primary API for all assemblers in PEACE.
    Additional documentation regarding the API provided by this class
    is available with the class and method documentation.
*/

/** The base class of all assemblers in PEACE.

    <p>This class must be the base class of all assemblers in the
    system. This class provides some default functionality that can be
    readily used by each assembler.  However, this class is very
    generic, thereby providing the necessary flexibility for
    implementing different types of assemblers.</p>

	<p>This class cannot be directly instantiated.  Instead a suitable
    derived assembler can be created via the AssemberFactory.  Refer
    to the documentation on the AssemblerFactory for additional
    details.</p>

    @see AssemblerFactory
    @see BatonAssembler
*/
class Assembler : public Component {
public:
   /** Add valid command line arguments for this assembler.

        This method must be used to add all valid command line options
        that are supported by this assembler.  Note that derived
        classes may override this method to add additional command
        line options that are applicable to it.  This method is
        invoked when the clustering sub-system is initialized.

        \note Derived assembler classes <b>must</b> override this
        method to display help for their custom command line
        arguments.  When this method is overridden don't forget to
        call the corresponding base class implementation to add common
        options.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);
	
    /** Primary method to perform cDNA assembly.

        This method is the primary interface method that is invoked
        (on each and every MPI process) to perform the complete
        assembly and output generation process.  This method is
        invoked once all the command line parameters are processed by
        the various subsystems constituting PEACE.  This method is a
        pure-virtual method.  Therefore all derived assembler classes
        must override this method to perform all the necessary
        operations.

        \note This method must be invoked only after the initialize()
        method has been invoked.

        \return This method returns zero if the assembly process was
        successfully completed.  Otherwise this method returns a
        non-zero error code.
    */
    virtual int assemble() = 0;

    /** Set the ESTList to be used by this assembler.

        This method must be used to set the list of cDNA fragments
        being processed by a given run of PEACE.  This list is
        typically a shared (between all components associated with a
        single instance of PEACE class) list of cDNA fragments that is
        used and processed by this class.

        \note This class maintains a pointer to this list internally
        until the finalize() method is called.  The methods in this
        class (and its child classes) do not add/remove entries in
        this list. They only update various attributes associated with
        each cDNA entry in the list.
        
        \param[in] estList The ESTList to be used by this analyzer.
        This pointer is used until the finalize() method is invoked.
     */
    void setESTList(ESTList* estList);
	
    /** A method to handle initialization tasks for the Assembler.
		
        This method is called after an assembler object has been
        instantiated but \e before the cDNAs (to be assembled) have
        been loaded.  This method performs the following tasks:

        <ol>

        <li> First it loads all the ESTs and other cDNA fragments to
        be processed from the specified input file.</li>

        <li>Next it creates the output file (if specified, otherwise
        output is dumped to standard out) to write assembler
        data.</li>
        
        </ol>
        
        \return This method returns zero on success. On errors, this
        method returns a non-zero value.
    */
    virtual bool initialize();

    /** Obtain a mutable pointer to the list of ESTs associated with
        this assembler.

        This method can be used to obtain a mutable list of ESTs that
        this assembler is operating with.  This method essentially
        returns a pointer to the shared RuntimeContext::estList
        object.

        \note This method returns a valid pointer only after the
        initialize() method has been invoked.
        
        \return A pointer to the list of ESTs associated with this run
        of PEACE.
    */
    inline ESTList* getESTList() { return estList; }

    /** Obtain a immutable pointer to the list of ESTs associated
        with this assembler.

        This method can be used to obtain a reference to the list of
        ESTs that this assembler is operating with.  This method
        essentially returns a pointer to the shared
        RuntimeContext::estList object.

        \return A constant pointer to the list of ESTs associated with
        this run of PEACE.
    */
    inline const ESTList* getESTList() const { return estList; }
    
    /** The destructor.

        The destructor frees memory allocated for holding any data in
        this base class.
    */
    virtual ~Assembler();

    /** Obtain the list of contigs that have already been built by
        this assembler.

        <p>This method can be used to obtain a immutable reference to
        the contigs that have been created by this
        assembler. Typically this method is used to obtain the final
        set of contigs that have been created at the end of
        gnomic-assembly process.</p>

        <p>It is important to note that, in parallel assembly
        scenarios, the contig information in this list typically has
        only the local information (pertaining to the process) and
        does not have global information.  This contig information
        remains distributed across the various parallel processes to
        minimize memory overheads.  Consequently, to get full contig
        information, the contig data from all the parallel processes
        must be appropriately fused together.</p>

        \return The list of contigs that have been created by this
        assembler.
    */
    inline const ContigList& getContigs() const { return contigList; }

    /** Add a contig listener to be notified whenever a new contig is
        successfully created by this assembler.

        This method can be used to register a contig listener with
        this assembler.  The contig listener is notified (by this
        assembler) whenever a new contig is successfully created by
        this assembler.  The assembler maintains a pointer to the
        contig listener until assembly operation is successfully
        completed.

        \param[in] cl The contig listener to be notified whenever this
        assembler successfully creates a new contig.
    */
    void addContigListener(ContigListener* cl);
    
protected:
    /** The default constructor.

        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived assembler classes must be instantiated via the
        AssemblerFactory API methods.

        \param[in] name The human readable name for this assembler.
        This name is used when generating errors, warnings, and other
        output messages for this object.
    */
    Assembler(const std::string& name);

    /** Helper method to write progress information to a given data
		file.

        This method is a helper method that is invoked from the main
        assemble method (in the derived class) to report progress logs
        as ESTs are analyzed and assembled.  This method cuts logs
        only if the progressFile is valid.  It logs the number of ESTs
        processed and the total number of fragments using the values
        returned by ESTList::getESTCount() and the number of entries
        reported as having been processed.

		\param[in] estsProcessed The number of cDNA fragments that
		have been processed.  This value is used to report the
		progress.
    */
    void reportProgress(const int estsProcessed);

    /** Helper method that is used (by derived classes) to notify
        contig listeners that a contig has been successfully created.

        <p>This method must be invoked by the derived class methods to
        notify all registered contig listeners that a contig has been
        successfully formed. This method must be invoked from all
        parallel-processes immaterial of whether the contig has full
        or partial information.</p>

        <p>This method iterates over the contig listeners in this
        class and invokes the interface method on each one of the
        registered contig listeners.</p>
        
        \param[in] contig The contig that has just been successfully
        created by this assembler.
        
        \param[in] fullContig If this flag is \c true then the contig
        has full information.  If this flag is \c false then this
        contig has partial information and data from other
        parallel-processes must be fused together to obtain full
        contig information.

        \return This method returns \c true if all the listeners
        reported \c true (indicating the contig has been successfully
        processed). Otherwise this method returns \c false.
    */
    bool notifyContigListeners(const Contig& contig, 
                               const bool fullContig) const;
    
    /** A shortcut reference to the shared list of cDNA fragments
        being analyzed.
		
        This list holds a pointer to the shared list of cDNA fragments
        currently being assembled.  This pointer is initialized to
        NULL in the constructor.  A valid pointer is filled in via the
        setESTList() method or when the initialize() method is invoked
        (the pointer to the ESTList is obtained from the
        SubSystem::runtimeContext via the Component::subSystem member
        in the base class).
    */
    ESTList *estList;

    /** The list of contigs that have been built by this assembler.

        <p>This instance variable contains the set of contigs that
        have been created by this assembler. Derived classes directly
        add new contigs to this list whenever a contig has been
        successfully created.</p>

        <p>It is important to note that, in parallel assembly
        scenarios, the contig information in this list typically has
        only the local information (pertaining to the process) and
        does not have global information.  This contig information
        remains distributed across the various parallel processes to
        minimize memory overheads.  Consequently, to get full contig
        information, the contig data from all the parallel processes
        must be appropriately fused together.</p>

        \return The list of contigs that have been created by this
        assembler.
    */
    ContigList contigList;
    
    /** Name of file to report progress in during assembly.
		
        This command line argument provides the name of the log file
		where progress information is to be written. The progress
		information is in the form: \#estsProcessed, \#ests. The file
		name is specified via a command line argument \c --progress.
		This value is used only by the BatonAssemblerManager class.
		Possibly this parameter must be moved down to the child class.
    */	
	std::string progFileName;

    /** A list of liteners to be notified whenever a contig has been
        successfully formed.

        This list contains a set of classes to be notified whenever a
        contig has been formed.  The listeners utlize the contig
        information to perform suitable operations. 

        @see SAMWriter
    */
    ContigListenerList contigListeners;
    
private:
    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \note This method is intentionally has a dummy implementation.
        
        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.
        
        \return Reference to this.
    */
    Assembler& operator=(const Assembler& src);

    /** File stream to log progress information.
        
        This output stream is created in the initialize() method if
        the progFileName is specified.  The progress information is
        generated by the updateProgress method.
    */
    std::ofstream progressFile;	
};

#endif
