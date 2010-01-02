#ifndef CLU_H
#define CLU_H

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

#include "FWAnalyzer.h"
#include "ESTCustomData.h"
#include <vector>

// Forward declaration to keep compiler happy
class EST;
class ResultLog;

/** CLU: An implementation of the CLU's analysis algorithm, originally
    developed by Andrey Ptitsyn and Winston Hide.

    <p>This analyzer implements the similarity metric generation part
    of CLU, a clustering algorithm developed by Andrey Ptitsyn and
    Winston Hide.  The Reference to clu is: <br><br>


    "CLU: A new algorithm for EST clustering", A. Ptitsyn and W. Hide,
    BMC Bioinformatics, 6(2), 2005. doi:
    10.1186/1471-2105-6-S2-S3. <br><br> </p>

    <p>The implementation in this class has been developed by suitably
    adapting parts of the source code for CLU available from
    http://lamar.colostate.edu/~ptitsyn/ </p>
    
    <p>This class has been implemented by extending the FWAnalyzer
    base class.  The FWAnalyzer base class provides most of the
    standard functionality involved in reading FASTA files and
    generating formatted output and processing some parameters.  This
    class adds functionality to compare EST's using CLU's similarity
    comparison metrics.</p>

    \note This class instantiated via the ESTAnalyzerFactory::create
    method.
*/
class CLU : public FWAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~CLU();

    /** Display valid command line arguments for this analyzer.

        This method must be used to display all valid command line
        options that are supported by this analyzer.  Note that
        derived classes may override this method to display additional
        command line options that are applicable to it.  This method
        is typically used in the main() method when displaying usage
        information.

        \note Derived EST analyzer classes <b>must</b> override this
        method to display help for their custom command line
        arguments.  When this method is overridden don't forget to
        call the corresponding base class implementation to display
        common options.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.
        
        This method is used to process command line arguments specific
        to this EST analyzer.  This method is typically used from the
        main method just after the EST analyzer has been instantiated.
        This method consumes all valid command line arguments.  If the
        command line arguments were valid and successfully processed,
        then this method returns \c true.

        \note Derived EST analyzer classes <b>must</b> override this
        method to process any command line arguments that are custom
        to their operation.  When this method is overridden don't
        forget to call the corresponding base class implementation to
        display common options.
        
        \param[in,out] argc The number of command line arguments to be
        processed.

        \param[in,out] argv The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.  This method checks to ensure that a valid
        frame size and a valid word size have been specified.
    */
    virtual bool parseArguments(int& argc, char **argv);

    /** Method to begin EST analysis.

        This method is invoked just before commencement of EST
        analysis.  This method first invokes the base class method
        that loads the list of ESTs from a given input multi-FASTA
        file and pouplates the list of ESTs.  If the ESTs were
        successfully loaded, then this method initializes the custom
        data for each EST (with empty hash maps).

        \return If the ESTs were successfully loaded from the FATA
        file then this method returns 0.  Otherwise this method
        returns with a non-zero error code.
    */
    virtual int initialize();

    /** Method to obtain human-readable name for this EST analyzer

        This method provides a human-readable string identifying the
        EST analyzer.  This string is typically used for
        display/debugging purposes (particularly via the PEACE
        Interactive Console).

        \return This method returns the string "CLU" identifiying this
        analyzer.
    */    
    virtual std::string getName() const { return "CLU"; }

    virtual float getValidMetric() const { return 1000; }
    
    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        builds the hash table (used by CLU for searching/comparison)
        for the reference analysis.  The reference EST is also called
        Sq in CLU literature.

        \note This method must be called only after the initialize()
        method is called.  This method overrides the implementation in
        the base class to perform its own custom operation.

        \return If the reference estIdx is invalid then this method
        returns with 1.  Otherwise it pouplates the referenceTable
        array with the necessary information and returns 0.
    */
    virtual int setReferenceEST(const int estIdx);
    
protected:    
    /** Analyze and obtain a similarity metric.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \note This method overrides the default implementation in the
        base class to perform its own custom operations.
        
        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.  This method
        essentially implements the search() method from CLU's
        implementation.

        \return This method must returns a similarity metric by
        comparing the ESTs by calling the analyze() method.
    */
    using FWAnalyzer::getMetric;
    virtual float getMetric(const int otherEST);
    
    /** Dumps a given EST in 3-column format using a ResultLog.

        This method is a helper method that dumps a given EST out to
        the log.  This method overrides the default implementation in
        the base class to perform its own custom operation.
        
        \param[out] log The log to which the EST is to be dumped.

        \param[in] est The EST to be dumped.  This parameter is never
        NULL.

        \param[in] isReference If this flag is true, then this EST is
        the reference EST to be dumped out.
    */
    virtual void dumpEST(ResultLog& log, const EST* est,
                         const bool isReference = false);

    /** Helper method to build reference and complement hash maps.

        This is a helper method that is used to build the reference
        and complement hash maps for a given EST.  If the reference
        (and complement) hash maps already exist for the given EST
        then this metod exits immediately without rebuilding hash
        maps.  Otherwise it builds the reference and complement hash
        maps using the createCLUHashMap() method and pouplates the
        custom ESTCLUData object associated with the EST.

        \param[in,out] est Pointer to the EST whose reference and
        complement hash maps are to be built.
    */
    void buildHashMaps(EST *est);
        
    /** Helper method to create the CLU hash/look up table.

        This method is invoked from the setReferenceEST() method to
        create the referenceHashTable and complementHashTable required
        for comparing and processing other ESTs with the reference
        EST.  This method operates as follows:

        <ol>

        <li>It initializes the table (if it is NULL) to hold
        4^wordSize hash entries.</li>

        <li>It resets all the entries to 0 (zero) in the table.</li>

        <li>For each wordSize base pairs in the sequence, it computes
        the hash value for the sequence and increments the
        corresponding entry (using hash value as the index) in the
        table.</li>

        <li>It finally calls the filterHashMap() method to filter
        out low-complexity and abundant oligos.</li>

        </ol>

        \param[in,out] table The hash table to be populated by this
        method.

        \param[in] sequence The EST sequence to be used for populating
        the hash table.
    */
    void createCLUHashMap(int* &table, const char *sequence);

    /** Helper method to filter out certain entries from the reference
        hash map.

        This helper method is invoked from the setReferenceEST()
        method to filter out certain entries from the hash map as they
        are not significant when comparing ESTs.  Specifically, this
        method filters out non-informative and low-complexity
        sequences in the following manner:

        <ul>

        <li> First, zero out all simple oligos from consideration.
        The simple oligos are sequences of the form:
        "AAAAAACCCCCCGGGGGGTTTTTT".</li>

        <li>Next it removes abundant sequences (sequence that occur
        too frequently) from the table.  Such abundant sequences are
        not informative in EST comparisons.</li>

        <li>Finally, all non-zero entries are normalized to 1.</li>
        
        </ul>

        \param[in,out] table The table of hash values to be filtered
        and normalized by this method.  This table must have been
        populated by a call to createCLUHashTable() method.
        
        \param[in] sequenceLength The length of the EST sequence from
        which the referenceHashMap has been generated.
    */
    void filterHashMap(int *table, const int sequenceLength);

    /** Obtain similarity metric between reference sequence and a
        given sequence.

        This method performs the core task of comparing a given EST
        sequence with the reference sequence given a CLU hash table.
        This method operates as follows:

        <ol>

        <li>It computes the hash for each word the first frame of the
        given sequence and obtains a categorized distribution (cd)
        index by add number of matching words found in the reference
        sequence.</li>

        <li>For sequent frames in the sequence, it uses the values
        computed for the first sequence to compute the categorized
        distribution (cd) index  value. </li>

        <li>It adds one to the corresponding cd entry as indicated by
        the resulting cdIndex value. </li>

        <li>Finally it determines the sum of multiplying the cd values
        with the thresholds (constant values obtained earlier by the
        original authors via Monte-Carlo type simulations) to obtain
        the similarity metric.</li>

        </ol>

        \param[in] hashTable The hash table from the reference
        sequence (either the referenceHashTable or
        complementHashTable) to be used for searching/comparison.

        \param[in] sequence The other sequence with which the
        reference sequence is to be compared.

        \return A similarity score for the given sequence against the
        reference sequence.
    */
    int getSimilarity(const int* const hashTable,
                      const char* const sequence) const;

    /** Container for CLU's hash maps for ESTs

        This class is used to encapsulate the reference and complement
        hash maps needed by CLU to compare ESTs.  This class extends
        ESTCustomData so that the hash maps can be directly associated
        with the EST they refer to.  New instances of this object are
        created in the initialize method.  The place holder objects
        are automatically deleted once analysis is complete.
    */
    class CLUESTData : public ESTCustomData {
    public:
        /** The default constructor.

            The default constructor merely initializes the
            referenceHashMap and complementHashMap to NULL.
        */
        CLUESTData() : referenceHashMap(NULL), complementHashMap(NULL) {}

        /** The destructor.

            The destructor deletes the reference and complement hash
            map in this class if they are not NULL.
        */
        virtual ~CLUESTData();
        
        /** Hash table to accelerate the comparison of all words in a
            sequence.
            
            <p> Quotation from CLU paper: "To accelerate the comparison all
            words found in sequence Sq are presented in a table (H). The
            table H is a linear array, where the offset itself is a hash
            value of certain oligonucleotide. Each element of this array
            contains a number, associated with a corresponding
            oligonucleotide".</p>
            
            This vector has a size 4^wordSize (for example, if wordSize is
            6, this table has 4096 entries).  These entries are populated
            in the setReferenceEST() method.
        */
        int *referenceHashMap;
        
        /** Hash table to accelerate the comparison of all words in a
            complement (reverse) sequence.
            
            <p> Quotation from CLU paper: "To accelerate the comparison all
            words found in sequence Sq are presented in a table (H). The
            table H is a linear array, where the offset itself is a hash
            value of certain oligonucleotide. Each element of this array
            contains a number, associated with a corresponding
            oligonucleotide".</p>
            
            This vector has a size 4^wordSize (for example, if wordSize is
            6, this table has 4096 entries).  These entries are populated
            in the setReferenceEST() method using the reverse of a given
            EST sequence.
        */
        int *complementHashMap;

    protected:
        // Currently the ESTCLUData class has no protected members.

    private:
        // Currently the ESTCLUData class has no private members.
    };
    
    /** Parameter to define fraction value to compute abundance
        metric.

        This instance variable is used to track the fraction of values
        (with respect to sequence length) after which a specific
        oligonucleotide (word) must be considered to be abundant.  The
        default value is 10.  This value can be changed by the user
        via a command line parameter.
    */
    static int abundanceFraction;
    
private:
    /** The set of arguments specific to the CLU program

        This instance variable contains a static list of arguments
        that are specific only to the CLU analyzer class.  This
        argument list is statically defined and shared by all
        instances of this class.

	\note Use of static arguments and parameters makes CLU class
	hierarchy not MT-safe.
    */
    static arg_parser::arg_record argsList[];    

    /** The default constructor.

        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ESTAnalyzerFactory is a
        friend of this class, an object can be instantiated via the
        ESTAnalyzerFactory::create() method.

	
	\param[in] refESTidx The reference EST index value to be used
	when performing EST analysis.  This parameter should be >= 0.
	This value is simply passed onto the base class.  This
	parameter is not really used for this analyzer is used for
	clustering.
        
	\param[in] outputFile The name of the output file to which the
	EST analysis data is to be written.  This parameter is ignored
	if this analyzer is used for clustering.  If this parameter is
	the empty string then output is written to standard output.
	This value is simply passed onto the base class.
    */
    CLU(const int refESTidx, const std::string& outputFile);

    /** A simple array to map characters A, T, C, and G to 0, 1, 2,
        and 3 respectively.

        This is a simple array of 255 entries that are used to convert
        the base pair encoding characters A, T, C, and G to 0, 1, 2,
        and 3 respectively to compute the hash as defined by CLU.
        This array is initialized in the constructor and is never
        changed during the life time of this class.
     */
    static char CharToInt[];
};

#endif
