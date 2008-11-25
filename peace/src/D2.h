#ifndef D2_H
#define D2_H

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

/** \brief EST Analyzer that uses the d2 algorithm to compute
    distances between two ESTs.

    <p>This analyzer provides the mechanism to use vanilla D2
    algorithm to compute the distance values between a pair of
    ESTs. The D2 implementation has been adapted from the
    implementations of WCD and CLU.<p>

    \note This D2 analyzer uses a word size of 8 base pairs.  This is
    hard coded into the algorithm in the form of various assumptions
    (data types of variables etc.). A similar approach is used by
    other D2 implementations as well.  This is a compromise between
    performance and flexibility of the implementation.
*/
class D2 : public ESTAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~D2();
    
    /** Display valid command line arguments for this analyzer.

        This method must be used to display all valid command line
        options that are supported by this analyzer.  Currently, this
        analyzer does not require any special command line parameters.

        \note The ESTAnalyzer base class requires that derived EST
        analyzer classes <b>must</b> override this method to display
        help for their custom command line arguments.  When this
        method is overridden don't forget to call the corresponding
        base class implementation to display common options.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.

        <p> This method is used to process command line arguments
        specific to this EST analyzer.  This method is typically used
        from the \c main method just after the EST analyzer has been
        instantiated.  This method consumes all valid command line
        arguments.  If the command line arguments were valid and
        successfully processed, then this method returns \c true. </p>

        <p>Currently, this EST analyzer does not require any
        additional command line parameters.  Consequently, it simply
        calls the corresponding method in the base class.</p>

        \note The \c ESTAnalyzer base class requires that derived EST
        analyzer classes <b>must</b> override this method to process
        any command line arguments that are custom to their operation.
        When this method is overridden don't forget to call the
        corresponding base class implementation to display common
        options.
        
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
        analysis.  This method currently does not have any specific
        tasks to perform.  It simply returns 0.

        \return Currently, this method always returns 0 (zero) to
        indicate initialization was successfully completed.
    */
    int initialize();

    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        currently saves the index in the instance variable for further
        look up.  Next it computes the D2 reference hash map for use
        in the analyze() method.  Note that the same D2 reference hash
        map is used for a batch of EST analysis to minimize the
        overhead of recomputing the hash map.  This is an
        implementation optimization and not an algorithmic one.

        \note This method must be called only after the initialize()
        method is called.

        \return This method returns \c true if the estIdx was within
        the given range of values.  Otherwise this method returns a
        non-zero value as the error code.
    */
    virtual int setReferenceEST(const int estIdx);

    /** Analyze and obtain a distance metric.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \param[in] otherEST The index (zero based) of the EST with
        which the reference EST is to be compared.

        \return This method returns the distance value reported by the
        D2 algorithm.
    */
    virtual float analyze(const int otherEST);
    
    /** Method to pefrom a batch of EST analysis.
        
        The \c ESTAnalyzer base class requires this method to be
        overloaded in the dervied class(es).  This method is used to
        perform the core tasks of EST analysis for the D2 EST
        analyzer.  This method operates in the following manner:

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

            <li> It extracts a frame from the reference EST and
            current EST. <li>

            <li> Next, it calls the polymorphic analyze() method to obtain
            similarity metric. </li>

            <li>It logs the similarity metric using suitable methods in
            the ESTAnalyzer base class.<li>

            </ol>
        
        <li>If all the processing proceeds successfully, this method
        returns 0 (zero).
        
        </ol>
    */
    virtual int analyze();
    
protected:
    /** Method to compare two metrics generated by this class.

        This method provides the interface for comparing metrics
        generated by this ESTAnalyzer when comparing two different
        ESTs.  This method returns \c true if \c metric1 is
        comparatively better than or equal to \c metric2.

        \note As per the ESTAnalyzer API requirements, EST analyzers
        that are based on distance measures (such as this D2 analyzer)
        \b must override this method.
        
        \param[in] metric1 The first metric to be compared against.

        \param[in] metric2 The second metric to be compared against.
        
        \return This method returns \c true if metric1 is
        comparatively better than \c metric2.
    */
    inline bool compareMetrics(const float metric1, const float metric2) const
    { return (metric1 < metric2); }
    
    /** Obtain an invalid (or the worst) metric generated by this
        analyzer.
        
        This method can be used to obtain an invalid metric value for
        this analyzer.  This value can be used to initialize metric
        values.
        
        \note Derived distance-based metric classes (such as this D2
        analyzer) \b must override this method to provide a suitable
        value.
        
        \return This method returns an invalid (or the worst) metric
        of 1e7 for this EST analyzer.
    */
    inline float getInvalidMetric() const { return 1e7; }

    /** Determine if this EST analyzer provides distance metrics or
        similarity metrics.
        
        This method can be used to determine if this EST analyzer
        provides distance metrics or similarity metrics.  If this
        method returns \c true, then this EST analyzer returns
        distance metrics (smaller is better).  On the other hand, if
        this method returns \c false, then this EST analyzer returns
        similarity metrics (bigger is better).

        \return This method returns \c true to indicate that this EST
        analyzer operates using distance metrics.
    */
    bool isDistanceMetric() const { return true; }
    
private:
    /* The default constructor for this class.
       
        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ESTAnalyzerFactory is a
        friend of this class; therefore it can instantiate the D2
        analyzer.  Accordingly, the ESTAnalyzerFactory::create()
        method must be used to instantiate this class.

	\param[in] refESTidx The reference EST index value to be used
	when performing EST analysis.  This parameter should be >= 0.
	This value is simply passed onto the base class.
        
	\param[in] outputFile The name of the output file to which the
	EST analysis data is to be written.  This parameter is ignored
	if this analyzer is used for clustering.  If this parameter is
	the empty string then output is written to standard output.
	This value is simply passed onto the base class.
    */
    D2(const int refESTidx, const std::string& outputFileName);
};

#endif
