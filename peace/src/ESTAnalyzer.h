#ifndef EST_ANALYZER_H
#define EST_ANALYZER_H

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

#include "arg_parser.h"

/** The base class of all EST analyzers.

    This class must be the base class of all EST analyzers in the
    system. This class provides some default functionality that can be
    readily used by the EST analyzers.
*/
class ESTAnalyzer {
public:
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
        returns \c false.  This method returns true if all arguments
        are consumed successfully and if a valid estID and estFileName
        have been specified.
    */
    virtual bool parseArguments(int& argc, char **argv);   

    /** Method to begin EST analysis.

        This method is invoked just before commencement of EST
        analysis.  This method typically loads the list of ESTs from a
        given input file.  In addition, it may perform any
        pre-processing as the case may be.

        \note Derived classes must override this method.

        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns with a
        non-zero error code.
    */
    virtual int initialize() = 0;

    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  Setting the
        reference EST provides analyzer's an opportunity to optimize
        certain operations, if possible.

        \note This method must be called only after the initialize()
        method is called.

        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns an error code.
    */
    virtual int setReferenceEST(const int estIdx) = 0;

    /** Analyze and obtain a similarity metric.

        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.

        \return This method returns a similarity metric by comparing
        the ESTs.
    */
    virtual float analyze(const int otherEST) = 0;
        
    /** Method to perform EST analysis.

        This method must be used to perform EST analysis.  This method
        is a pure-virtual method.  Therefore all EST analyzers must
        override this method to perform all the necessary operations.
        Typically, this method performs the following operations:

        <ol>

        <li>This method calls initialize.</li>

        <li>Set's the reference EST via a call to the setReferenceEST()
        method.</li>

        <li>Repeatedly uses the analyze(const int) method to compare
        ESTs.</li>

        <li>Generates analysis reports at the end of analysis.</li>

        </ol>
    */
    virtual int analyze() = 0;

    /** Method to load EST information from a FASTA file.

        This method can be used to load information regarding ESTs
        from a FASTA file.  The file name from where the data is to be
        loaded must be passed in as the parameter.

        \param[in] fileName The file name of the FASTA file from where
        the EST information is to be uploaded.

        \param[in] unpopulate If this parameter is true then the
        header and sequence information in each EST is discarded to
        minimize memory foot print.

        \return This method returns true if all the ESTs were
        successfully loaded from the given file.
    */
    bool loadFASTAFile(const char *fileName, const bool unpopulate = false);

    /** Obtain the input file name.

	This method returns the input file from where the EST data was
	read.

	\return The input file from where the EST data was read.  If
	an input file was not specified, then this method returns
	NULL.
     */
    const char* getInputFileName() const { return estFileName; }

    /** Determine if this EST analyzer provides distance metrics or
        similarity metrics.

        This method can be used to determine if this EST analyzer
        provides distance metrics or similarity metrics.  If this
        method returns \c true, then this EST analyzer returns
        distance metrics (smaller is better).  On the other hand, if
        this method returns \c false, then this EST analyzer returns
        similarity metrics (bigger is better).

        \note Derived classes that operate using distance metrics must
        overload this method to return \c true.
        
        \return This method returns \c false (by default) to indicate
        that this EST analyzer operates using similarity metrics.  If
        it operates using distance metrics then this method returns \c
        true.
    */
    virtual bool isDistanceMetric() const { return false; }

    /** Obtain an invalid (or the worst) metric generated by this
	analyzer.

	This method can be used to obtain an invalid metric value for
	this analyzer.  This value can be used to initialize metric
	values. By default this method returns -1, which should be
	ideal for similarity-based metrics.

	\note Dervied distance-based metric classes must override this
	method to provide a suitable value.

	\return This method returns an invalid (or the worst) metric
	for this EST analyzer.
    */
    virtual float getInvalidMetric() const { return -1; }

    /** Method to compare two metrics generated by this class.

        This method provides the interface for comparing metrics
        generated by this ESTAnalyzer when comparing two different
        ESTs.  This method returns \c true if \c metric1 is
        comparatively better than or equal to \c metric2.

	\note EST analyzers that are based on distance measures \b
	must override this method.
	
        \param[in] metric1 The first metric to be compared against.

        \param[in] metric2 The second metric to be compared against.

        \return This method returns \c true if metric1 is
        comparatively better then or equal to \c metric2.
    */
    virtual bool compareMetrics(const float metric1, const float metric2) const
    { return (metric1 > metric2); }
    
    /** The destructor.

        The destructor frees memory allocated for holding any EST data
        in the base class.
    */
    virtual ~ESTAnalyzer();
    
protected:
    /** The default constructor.

        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived ESTAnalyzer classes must be instantiated via the
        ESTAnalyzerFactor API methods.

        \param[in] name The human readable name for this EST analyzer.
        This name is used when generating errors, warnings, and other
        output messages for this analyzer.

	\param[in] refESTidx The reference EST's index in a given
	multi-FASTA file.  Index values start with 0 (zero).  The
	refESTidx is supplied as a global argument that is processed
	in the main() method.  This value is simply copied to the
	refESTidx member in this class.

	\param[in] outputFileName The file name to which output must
	be written.  If a valid output file is not specified, then
	results are written to standard output.  The outputFileName is
	simply copied to the outputFileName member object.
    */
    ESTAnalyzer(const std::string& analyzerName, const int refESTidx,
		const std::string& outputFileName);
    
    /** Flag to indicate if a read ahead thread must be used.

        This boolean value is by default set to false. However, the
        value is changed by the parseArguments method depending on
        wether the use whishes to use a read-ahead feature.
    */
    static bool readAhead;

    /** The index of the reference EST in a given file. 

        This member object is used to hold the index of a reference
        EST in a given file.  The index values begin from 0 (zero).
        This member is initialized in the constructor and is changed
        by the setReferenceEST() id.
    */
    int refESTidx;

    /** The FASTA file from where EST data is to be read.

        This member object is used to hold the file name from where
        all the EST data is to be loaded.  This member is initialized
        in the constructor and is never changed during the life time
        of this class.
    */
    static char* estFileName;

    /** Flag to indicate if output results must be in HTML format.

	This member is initialized to false.  However, the value is
        changed by the parseArguments method depending on the actual
        value specified by the user.
    */
    static bool htmlLog;

    /** The file to which results must be written.

        This member object is used to hold the file name to which all
        the analysis results are to be written.  This member is
        initialized to NULL.  However, the value is changed by the
        parseArguments method depending on the actual value specified
        by the user.
    */
    const std::string outputFileName;
    
    /** The name of this analyzer.

        This instance variable contains the human recognizable name
        for this analyzer.  This value is set when the analyzer is
        instantiated (in the constructor) and is never changed during
        the life time of this analyzer.  This information is used when
        generating errors, warnings, and other output messages.
    */
    const std::string analyzerName;

private:
    /** The set of common arguments for all EST analyzers.

        This instance variable contains a static list of arguments
        that are common all the EST analyzers.  The common argument
        list is statically defined and shared by all EST instances.

	\note This makes ESTAnalyzer class hierarchy not MT-safe.
    */
    static arg_parser::arg_record commonArgsList[];

    /** A dummy operator=

        The operator=() is supressed for this class as it has constant members
        whose value is set when the object is created.  These values cannot be
        changed during the lifetime of this object.

        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.

        \return Reference to this.
    */
    ESTAnalyzer& operator=(const ESTAnalyzer& src);
};

#endif