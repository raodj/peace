#ifndef BATON_LIST_CACHE_H
#define BATON_LIST_CACHE_H

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

#include <map>

// Forward declarations to keep compiler fast and happy
class BatonList;
class EST;
class ArgParser;

/** A simple class to maintain an in-memory cache of a given set of
    BatonList objects.

    This class provides a convenient interface to manage a set of
    BatonList objects.  For each EST being analyzed, this class
    maintains two BatonList entries.  One BatonList corresponds to the
    forward (or normal) nucleotide sequence.  The other entry
    corresponds to the reverse complementary representation.  Entries
    are created on-demand --- that is, when a specific BatonList is
    requested and the cache does not have it, it is computed.
*/
class BatonListCache {
public:
    /** The constructor.

        The constructor does not have any special tasks to perform and
        merely present to adhere to coding conventions.
    */
    BatonListCache();

    /** The destructor.

        The destructor frees up the BatonList objects maintained in
        the internal cache.
    */
    ~BatonListCache();

    /** Add the set of command line parameters for this component.
        
        This method is invoked by the BatonAssembler (that logically
        owns the cache) when it is being requested for command line
        arguments.  This method is expected to add the various command
        line arguments that can be used to further customize the
        operations of the cache.  This method adds the \c --nMers and
        \c --window command line parameters.

        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    void addCommandLineArguments(ArgParser& argParser);
    
    /** Get the size (in nucleotides) of the overlapping windows used
		for searching identical batons.

        This method can be used to get the windowSize that has been
        set when a BatonList object is created by this cache.  This
        method is used by the BatonAnalyzer and the SequenceAligner
        class hierarchies.

        \return The window size to be used when this cache creates
        BatonList objects.
    */
    inline int getWindowSize() const { return windowSize; }

    /** Get the size (in nucleotides) of baton heads/ends.

        This method can be used to get the number of nucleotides that
        constitute the baton heads/ends.  This size is used to create
        the batons in the BatonList objects cached by this class.
        This method is used by the BatonAnalyzer and the
        SequenceAligner class hierarchies.

        \return The size (in number of nucleotide bases) of the baton
        heads.
    */
    inline int getBatonHeadSize() const { return nMerSize; }
    
    /** Helper method to create and obtain normal and reverse
        complement baton lists for future reference.

        This method is an internal helper method that is used to
        create the necessary baton list(s) in a lazy-manner.  This
        method is invoked from different classes that share the cache
        (such as: BatonAnalyzer::analyze(),
        BatonAnalyzer::getMetric(), SequenceAligner::align()) methods
        to obtain the list of batons for processing.  This method
        first checks to see if the requested baton list is already
        available in the blCache (read as: baton list cache).  If not,
        this method builds the required baton list, updates the
        blCache, and returns the requested baton list.

        \param[in] est A pointer to an immutable EST object whose
        baton list is to be returned by this method.  This parameter
        cannot be NULL.
        
        \param[in] getRC If this flag is \c true, then this method
        returns the baton list corresponding to the reverce
        complementary sequence of the specified est.

        \return This method returns the requested baton list for
        further use. Note that the returned object must \b not be
        mutated in any manner as the data is cached.
    */
    BatonList* getBatonList(const EST* est, const bool getRC);
    
protected:
	/** The local cache of Baton Lists.

		This instance variable blCache (read as: baton list cache) is
		used to maintain a cache of baton list objects that have been
		pre-computed via an on-demand approach.  The entries in this
		list consist of baton lists for both normal and reverse
		complement (RC) representation of a given cDNA fragment.  The
		normal and RC baton lists for a cDNA fragment with ID \i k are
		stored in \c blCache with key value <i>k*2</i> and
		<i>k*2+1</i> respectively.
	*/
	std::map<int, BatonList*> blCache;    

    /** Setup the window size value to be used.

        This instance variable tracks the window size to be used when
        generating the BatonList objects.  This member is used to hold
        the window size (in number of nucleotides) into which a cDNA
        fragment must be subdivided to search for identical batons.
        Note that this value is a suggested average value.  Each baton
        list uses a slightly different value (around this average) as
        its window size. The window size plays an important role in
        effective identification of related cDNA fragments for
        clustering and assmebly.  The default value is 100.  The
        threshold value is a command line argument that can be set by
        the user via the \c --window command line argument.
    */
    int windowSize;

    /** Instance variable to track the number of nucleotides to be
		used for baton ends/heads.
		
        This member is used to hold the size (in nt) of the baton
        heads.  If this value is 3, then the baton head consists of
        3-mers of the form \c ATC, \c AAA, \c CGT etc.  The default
        value is 3.  This value can be changed via the \c --nMers
        command line parameter.
    */
    int nMerSize;
    
private:
    /** A dummy operator=
        
        The operator=() is supressed for this class as we don't have
        the cache to be copied or duplicated.

        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.

        \return Reference to this.
    */
    BatonListCache& operator=(const BatonListCache& src);

    /** A dummy copy constructor.

        The copy constructor has been made private to ensure that
        BatonListCache objects are never duplicated (even
        accidentally) to ensure cache coherence and minimal memory
        footprint.

        \param[in] src The source object from where the data is to be
        copied.
    */
    BatonListCache(const BatonListCache& src);    
};

#endif
