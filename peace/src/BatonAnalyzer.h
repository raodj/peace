#ifndef BATON_ANALYZER_H
#define BATON_ANALYZER_H

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

#include "ESTAnalyzer.h"
#include "Baton.h"
#include "BatonESTData.h"
#include <string>
#include <vector>

// Forward declaration to keep compiler happy
class EST;
class ResultLog;

// TODO: write documentation/get rid of FWAnalyzer traces
// also, this is in no shape to be inherited from

/** BatonAnalyzer: Frame and Word based Analyzer.

    <p>This analyzer provides a common base class for all EST
    analyzers that use a frame and word concept for analyzing ESTs.
    The total number of base pairs to be compared is called a Frame.
    A frame is broken into a sequence of fixed size (in bp) words.
    The frame size and word size (in terms of number of base pairs) is
    specified as command line arguments.

    <p>This class has been implemented by extending the ESTAnalyzer
    base class.  The ESTAnalyzer base class provides most of the
    standard functionality involved in reading FASTA files and
    generating formatted output.  This class adds functionality to
    compare EST's using the concept of frames and words</p>

    \note This class is never directly instantiated.  Instead, one of
    the derived classes are instantiated (via the
    ESTAnalyzerFactory::create) method and used.
*/
class BatonAnalyzer : public ESTAnalyzer {
  friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~BatonAnalyzer();
    
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
        
        \param[inout] argc The number of command line arguments to be
        processed.

        \param[inout] argc The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.  This method checks to ensure that a valid
        frame size and a valid word size have been specified.
    */
    virtual bool parseArguments(int& argc, char **argv);

    /** Method to begin EST analysis.

        This method is invoked just before commencement of EST
        analysis.  This method loads the list of ESTs from a given
        input multi-FASTA file and pouplates the list of ESTs.

        \return If the ESTs were successfully loaded from the FATA
        file then this method returns 0.  Otherwise this method
        returns with a non-zero error code.
    */
    virtual int initialize();

    virtual std::string getName() const { return "baton"; };

    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        extracts the start frame from the reference EST and sets it in
        the referenceFrame instance variable in this class.

        \note This method must be called only after the initialize()
        method is called.

        \return If the extraction of the reference EST frame was
        successful, then this method returns 0.  Otherwise this method
        returns an error code.
    */
    virtual int setReferenceEST(const int estIdx);
    
    /** Method to begin EST analysis.
        
        This method is used to perform the core tasks of EST analysis
        for all BatonAnalyzer classes.  This method operates in the
        following manner:

        <ol>

        <li>First it loads the necessary EST information from the
        supplied FASTA file using the initialize() method.  If the EST
        data is not successfully loaded then this method returns right
        away with 1.<li>

        <li>Upon successfully loading the EST data, the reference EST
        is set via the setReferenceEST() method.  If the reference EST
        is not correctly determined, then this method immediately
        returns with 2.</li>

        <li>For each EST in the list of ESTs it performs the following
        tasks:

        <ol> <li> It extracts a frame from the reference EST and
        current EST. <li>

        <li> Next, it calls the polymorphic analyze() method to obtain
        similarity metric. </li>

        <li>It logs the similarity metric using suitable methods in
        the ESTAnalyzer base class.<li>

        </ol>

        <li>If all the processing proceeds successfully, this method
        returns 0 (zero).
        
        </ol>

        \return This method returns zero if all the processing
        proceeded successfully. On errors this method returns a
        non-zero value.
    */
    virtual int analyze();
    
protected:    
    /** Analyze and obtain a similarity metric.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.

        \return This method must returns a similarity metric by
        comparing the ESTs by calling the analyze() method.
    */
    virtual float getMetric(const int otherEST);

    virtual float getValidMetric() const { return 100; }

    virtual void resetBatonCollector();

    virtual void buildESTData(BatonESTData** data, const int estIdx);

    int alignBetweenSections(int sect1, int sect2, int totalAllowedErrors,
			     int* begin1, int* begin2);

    bool quickCheckAlignment(int* est1, int* est2, int beg1, int beg2);

    int alignESTs(int* est1, int* est2, int est1size, int est2size,
		  int beg1, int beg2, int nTotalMutationError);

    int rightExtension(int nMer, int* est1, int* est2, int beg1, int beg2,
		       int end1, int end2, int rStepSoFar);

    int leftExtension(int nMer, int* est1, int* est2, int beg1, int beg2,
		      int lStepSoFar);

    int** highestFrequencySection(int** commonBatonFreqs, int criticalV,
				  int* arrLength, const int cbfLength,
				  const int cbfSubLength);

    void commonBatons(int** identicalBatonNumber, int n1, int n2);

    static int nMer;

    static int sectionWidth;

    static int maxErrorAllowed;

    int maxMerValue;

    int* lastPos;

    std::vector<Baton>* batonArr;

    BatonESTData* refESTdata;

    BatonESTData* otherESTdata;
    
    /** The default constructor.

        The default constructor for this class.  The constructor is
        made protected so that this class cannot be directly
        instantiated.  Instead one of the derived analyzer classes
        must be instantiated (via the ESTAnalyzerFactory::create())
        method and used.

        \param[in] name The human readable name for this EST analyzer.
        This name is used when generating errors, warnings, and other
        output messages for this analyzer.  This value is simply
        passed-on to the base class without any checks.

	\param[in] refESTidx The reference EST index value to be used
	when performing EST analysis.  This parameter should be >= 0.
	This value is simply passed onto the base class.
	
	\param[in] outputFile The name of the output file to which the
	EST analysis data is to be written.  This parameter is ignored
	if this analyzer is used for clustering.  If this parameter is
	the empty string then output is written to standard output.
	This value is simply passed onto the base class.	
    */
    BatonAnalyzer(const int refESTidx, const std::string& outputFile);

private:   
    /** The set of common arguments for all BatonAnalyzer instances.

        This instance variable contains a static list of arguments
        that are common all the Frame-Word analyzers.  The common argument
        list is statically defined and shared by all EST instances.

	\note This makes BatonAnalyzer class hierarchy not MT-safe.
    */
    static arg_parser::arg_record commonArgsList[];    
};


#endif
