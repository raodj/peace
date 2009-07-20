#ifndef FW_ANALYZER_H
#define FW_ANALYZER_H

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include "ESTAnalyzer.h"
#include <string>
#include <vector>

// Forward declaration to keep compiler happy
class EST;
class ResultLog;

/** FWAnalyzer: Frame and Word based Analyzer.

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
class FWAnalyzer : public ESTAnalyzer {
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~FWAnalyzer();

    /** Obtains a frame number of bp from a given EST sequence.

        This method is used to obtain a frame number of base pairs
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
    virtual void dumpESTList(const std::vector<EST*>& estList,
                             const EST* refEST,
                             ResultLog& log);

    /** Dumps a given EST in 3-column format using R

        This method is a helper method that dumps a given EST out to
        the log.

        \param[out] The log to which the EST is to be dumped.

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
        the parseArguments method depending on the actual value
        specified by the user.
    */
    static int frameSize;

    /** The word size to be used by this analyzer.

        The word size (in bp) that must be used for comparisons.  The
        default value is set to 0. However, the value is changed by
        the parseArguments method depending on the actual value
        specified by the user.

        \note The word size must be smaller than the frame size.
    */    
    static int wordSize;

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
    FWAnalyzer(const std::string& analyzerName, const int refESTidx,
	       const std::string& outputFile);

private:   
    /** The set of common arguments for all FWAnalyzer instances.

        This instance variable contains a static list of arguments
        that are common all the Frame-Word analyzers.  The common argument
        list is statically defined and shared by all EST instances.

	\note This makes FWAnalyzer class hierarchy not MT-safe.
    */
    static arg_parser::arg_record commonArgsList[];
};


#endif
