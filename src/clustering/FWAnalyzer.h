#ifndef FW_ANALYZER_H
#define FW_ANALYZER_H

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
#include <string>
#include <vector>

// Forward declaration to keep compiler happy
class EST;
class ResultLog;
class ESTList;

/** FWAnalyzer: Frame and Word based Analyzer.

    <p>This analyzer provides a common base class for all EST
    analyzers that use a frame-and-word concept for analyzing ESTs.
    Frames are also called <i>windows</i> and these two terms are
    synonymously used throughout PEACE.  A contiguous subsequence of a
    given cDNA fragment is called a Frame, aka window.  A frame is
    broken into a sequence of fixed size (in bp) words.  The frame
    size and word size (in terms of number of base pairs) is typically
    specified as command line arguments.</p>

    <p>This class has been implemented by extending the ESTAnalyzer
    base class.  The ESTAnalyzer base class provides most of the
    standard Application Program Interface (API) to be supported by
    all ESTAnalyzers.  This class adds functionality to compare EST's
    using the concept of frames and words</p>

    \note This class is never directly instantiated.  Instead, one of
    the derived classes are instantiated (via the
    ESTAnalyzerFactory::create()) method and used.
*/
class FWAnalyzer : public ESTAnalyzer {
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~FWAnalyzer();

    /** Obtains a frame number of bp from a given EST sequence.

        This method can used to obtain a frame number of base pairs
        from a given EST sequence.  The frame is extracted either from
        the beginning or the end of an EST depending on the start
        parameter.

        \param[in] est The EST from which a frame number of bp must be
        extracted.

        \param[in] start If this flag is \c true then this method
        extracts base pairs from the beginning of the EST
        sequence. Otherwise a frame is extracted from the end of the
        EST sequence.
    */
    std::string getFrame(const EST* est, bool start = true);
    
    /** Display valid command line arguments for this analyzer.

        This method must be used to add all valid command line options
        that are supported by this analyzer.  Note that derived
        classes may override this method to display additional command
        line options that are applicable to it.  This method is
        typically called when the clustering sub-system is
        initialized.

        \note Derived EST analyzer classes <b>must</b> override this
        method to display help for their custom command line
        arguments.  When this method is overridden don't forget to
        call the corresponding base class implementation to display
        common options.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);

    /** Method to begin EST analysis.
        
        <p>This method is invoked when the SubSystem (that logically
        owns this component) is being initialized.  It is invoked just
        before commencement of EST analysis.  This method validates
        command-line parameters, establishes the ESTAnalyzer::estList
        convenience pointer and initializes the heuristic chain
        established for this analyzer.</p>

        \note When derived classes override this method, they must
        call the base class implementation to ensure that general
        initialization is performed.
        
        \return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    virtual bool initialize();
    
    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        extracts the start frame from the reference EST and sets it in
        the referenceFrame instance variable in this class.

        \note This method must be called only after the initialize()
        method is called.

		\param[in] est A pointer to an immutable EST object to be used
		as the reference.  The reference EST will be subsequently
		analyzed with a number of other ESTs via the analyze() method.
		
        \return If the extraction of the reference EST frame was
        successful, then this method returns 0.  Otherwise this method
        returns an error code.
    */
    virtual int setReferenceEST(const EST *est);
    
    /** Method to begin EST analysis.
        
        This method is used to perform the core tasks of EST analysis
        for all FWAnalyzer classes.  This method operates in the
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

        <ol>

        <li> It extracts a frame from the reference EST and current
        EST. <li>

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

    /** Determine preferred dummy EST lengths to be used with this
        analyzer.

        \note For more detailed description of the motivation for
        dummy ESTs please refer to the documentation for the
        corresponding method in the base class --
        getPreferredDummyESTLength().
        
        \return This method overrides the default implementation in
        the base class to return twice the length of the frame (aka
        window) size.
    */
    virtual int getPreferredDummyESTLength() const
    { return frameSize * 2; }
    
protected:    
    /** Analyze and obtain a similarity metric.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \param[in] otherEST Pointer to an immutable EST object with
        which the reference EST is to be compared.

        \return This method must returns a similarity metric by
        comparing the ESTs by calling the analyze() method.
    */
    virtual float getMetric(const EST* otherEST);
    
    /** Method to compare two frames and compute similarity.

        This method must be overridden by derived Frame-Word analyzers
        (see FMWSCA.h) to compare two frames and report a similarity
        metric.

        \param[in] refFrame The reference frame for comparison
        purposes.  Note that the reference frame is always a constant
        in a given set of caparisons.  Consequently, certain analyzers
        can pre-compute and reuse metrics to make analysis fast.

        \param[in] otherFrame The other frame for comparison.  This
        frame is always guaranteed to be from a different EST than the
        refFrame.

        \param[in] wordSize The size of a word within the given frame.
        This value is always greater than 0 (zero) and less than frame
        size.
        
        \return This method is expected to return a similarity metric
        between the given frame and the refFrame.

        \note The default implementation of this method simply returns
        0.  Derived FWAnalyzer-based classes must override this method
        to perform the necessary operations.
    */
    virtual float analyzeFrame(const std::string& refFrame,
                               const std::string& otherFrame,
                               const int wordSize);
    
    /** Helper method to dump result log header.

        This is a helper method that is invoked from the analyze()
        method to dump a result log header.  This method was
        introduced to keep the code clustter in the analyze method to
        a minimum.

        This method dumps some of the analysis parameters to the
        supplied log.

        \param[out] log The log to which the header is to be dumped.

        \param[in] mean The overall mean similarity for this set of
        ESTs.

        \param[in] variance The overall variance in similarity for the
        given set of ESTs currently analyzed.
    */
    virtual void dumpHeader(ResultLog& log, const double mean,
                            const double variance);
    
    /** Helper method to dump post analysis EST list to a log.
        
        This is a helper method that is invoked from the analyze()
        method to dump the list of analyzed ESTs to a log.  This
        method was introduced to keep the code clustter in the analyze
        method to a minimum.  In addition, it provides the derived
        classes a chance to customize the working of the class.
       
        This method dumps the list of ESTs to the supplied log.

        \param[in] estList The list of ESTs that must be dumped out.

        \param[in] refEST The reference EST.
        
        \param[out] log The log to which the header is to be dumped.
    */    
    virtual void dumpESTList(const ESTList& estList,
                             const EST* refEST,
                             ResultLog& log);

    /** Dumps a given EST in 3-column format.

        This method is a helper method that dumps a given EST out to
        the log.
		
        \param[out] log The log to which the EST is to be dumped.
		
        \param[in] est The EST to be dumped.  This parameter is never
        NULL.
		
        \param[in] isReference If this flag is true, then this EST is
        the reference EST to be dumped out.
    */
    virtual void dumpEST(ResultLog& log, const EST* est,
                         const bool isReference = false);
    
    /** The frame size to be used by this analyzer.

        The frame size (in bp) that must be used for comparisons.  The
        default value is set to 0. However, the value is changed by
        command-line arguments depending on the actual value specified
        by the user.
        
        \note Derived (child) analyzers are not in any way required to
        utilize the user-supplied frame size -- see specific analyzer
        classes for those details.
    */
    int frameSize;

    /** The word size to be used by this analyzer.
        
        The word size (in bp) that must be used for comparisons.  The
        default value is set to 0. However, the value can be changed
        by the user method depending on the actual value specified via
        command-line arguments.

        \note The word size must be smaller than the frame size.
    */    
    int wordSize;
   
    /** The reference frame to be used for EST comparisons.

        This instance variable is set to the reference frame once the
        setReferenceFrame() method is called.  This referenceFrame is
        used in subsequent analysis() methods.
    */
    std::string referenceFrame;
    
    /** The default constructor.

        The default constructor for this class.  The constructor is
        made protected so that this class cannot be directly
        instantiated.  Instead one of the derived analyzer classes
        must be instantiated (via the ESTAnalyzerFactory::create())
        method and used.

        \param[in] analyzerName The human readable name for this EST
        analyzer.  This name can be used when generating errors,
        warnings, and other output messages for this analyzer.  This
        value is simply passed-on to the base class without any
        checks.
    */
    FWAnalyzer(const std::string& analyzerName);

    // These two need to be moved to the output sub-system.
    bool htmlLog;
    std::string outputFileName;
	
private:
	// Currently this class does not have any private members.
};

#endif
