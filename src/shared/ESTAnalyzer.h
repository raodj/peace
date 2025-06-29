#ifndef EST_ANALYZER_H
#define EST_ANALYZER_H

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
#include "Utilities.h"

// Forward declarations
class HeuristicChain;
class ESTList;
class EST;
class ESTAnalyzer;

/**
   This is a synonym for an ESTAnalyzer. Earlier, PEACE was primarily
   used for EST analysis. However, over time this functionality/scope
   has steadily increased to include other analysis such as
   clustering, assembly, and clustering using completely different
   metrics/measures. Accordingly, the source code is continuing to
   evolve and change with transitions in meanings of the different
   APIs (includes classes and methods).
*/
using Analyzer = ESTAnalyzer;

/** The base class of all EST analyzers.

    <p>This class must be the base class of all EST analyzers in the
    system. This class provides some default functionality that can be
    readily used by the EST analyzers. An ESTAnalyzer (is a Component)
    provides an interface to compare two given cDNA fragments.  The
    result of comparison is a quantitative pseudo-metric that provides
    a measure of either similarity (how close two cDNA fragments are)
    or distance (how dissimilar the given pair of CDNA fragments are).
    This information is used for filtering, clustering, and assembly
    operations in PEACE.  The ESTAnalyzer classes form a central
    component of PEACE.</p>

    <p>Note that ESTAnalyzer objects are used for a variety of
    operations by various sub-systems including:

    <ul>

    <li>In the filtering sub-system, some of the filters (particularly
    the LCFilter) uses analyzers to detect cDNA fragments with
    low-complexity regions.</li>

    <li>The clustering sub-system uses ESTAnalyzer class to identify
    similar cDNA fragments and cluster them together.</li>

    <li>The assembly sub-system uses this API to identify similar cDNA
    fragments to assemble them together to form a contig.</li>

    <ul>

    </p>

    <p>Since a single instnace of this class is shared by several
    independent sub-systems, the various API methods are invoked
    several times from each sub-systems.  For example, the
    initialize() method maybe inovked more than once.  Consequently,
    it maybe important to check for redundant invocations of the
    initialize() method and short circuit initialization
    operation.</p>


    <p>Since the ESTAnalyzer::initialize() method will be invoked
    multiple times from various sub-systems, derived classes can
    optionally short circuit unnecessary intialization operations in
    the following manner:

    \code

    if (isInitialized()) {
        // Already initialized. Nothign further to be done.
        return true;
    }

    \endcode

    </p>

*/
class ESTAnalyzer : public Component {
public:
    /** Add valid command line arguments for this analyzer.

        This method must be used to add all valid command line options
        that are supported by this analyzer.  Note that derived
        classes may override this method to add additional command
        line options that are applicable to it.  This method is
        invoked when the clustering sub-system is initialized.

        \note Derived EST analyzer classes may override this method to
        display help for their custom command line arguments.  When
        this method is overridden don't forget to call the
        corresponding base class implementation to add common options.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);
    
    /** Method to begin EST analysis.

        <p>This method is invoked when the SubSystem (that logically
        owns this component) is being initialized.  It is invoked just
        before commencement of EST analysis.  This method ensures that
        the ESTAnalyzer::estList convenience pointer is set and
        initializes the heuristic chain established for this
        analyzer. As an added convenience, this method also propagates
        the ESTAnalyzer::estList pointer to the heuristic chain, if
        the EST list has not already been set for the heuristic
        chain.</p>

        <p>Since the ESTAnalyzer::initialize() method will be invoked
        multiple times from various sub-systems, derived classes that
        override the initialize method can short circuit unnecessary
        initialization in the the following manner:
        
        \code

		// Do any reinitialization that may be needed here.
		
        if (isInitialized()) {
            // Already initialized. Nothign further to be done.
            return true;
        }

		// Do one time initialization after the above check.
		
        \endcode
        
        \note When derived classes override this method, they must
        call the base class implementation to ensure that general
        initialization is performed.  Furthermore, call the base-class
        implementation only if the isInitialized() method returns \c
        false (to avoid redundant initializations).
        
        \return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    virtual bool initialize();

    /** Set the ESTList to be used by this analyzer.

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
    
    /** Set the reference EST id for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  Setting the
        reference EST provides analyzer's an opportunity to optimize
        certain operations, if possible.

        \note This method must be called only after the initialize()
        method is called.

		\param[in] est A pointer to an immutable EST object to be used
		as the reference.  The reference EST will be subsequently
		analyzed with a number of other ESTs via the analyze() method.
		
        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns an error code.
    */
    virtual int setReferenceEST(const EST* est) = 0;

    /** Wind-up the operations of this analyzer.

        This method is invoked when the SubSystem (that logically owns
        this analyzer) is being finalized.  This method is expected to
        perform any cleanup operations (converse of operations
        performed in the initialize() method).  This method resets the
        internal ESTAnalyzer::estList pointer to NULL.
    */
    virtual void finalize();
    
    /** Analyze and obtain a similarity metric using the attached
        heuristic chain (if one exists) followed by the appropriate
        heavy weight distance/similarity measure associated with
        this ESTAnalyzer.
		
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.
	
        \note This method may return -1, if the otherEST is
        significantly different from the reference EST (possibly
        warranting no further analysis) that a meanigful metric cannot
        be generated.
		
        \param[in] otherEST A pointer to an immutable EST object with
        which the reference EST is to be compared.
		
        \param[in] useHeuristics A directive instructing the ESTAnalyzer
        on whether or not to use its heuristis chain.  Defaults to true.
		
        \param[in] useHeavyWeight A directive instructing the ESTAnalyzer
        on whether or not to use the heavy weight metric.  Defaults to true.
		
        \return This method returns a similarity/distance metric by
        comparing the ESTs. This method may return -1, if the otherEST
        is significantly different from the reference EST (possibly
        warranting no further analysis) that a meanigful metric cannot
        be generated.
    */
    float analyze(const EST* otherEST, const bool useHeuristics = true,
				  const bool useHeavyWeight = true);

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

    /** Get alignment data for the previous call to analyze method.

        This method can be used to obtain alignment data (if any) that
        was obtained typically as an byproduct of the previous call
        tothe analyze() method.

        \param[out] alignmentData The parameter is updated to the
        alignment information generated as a part of the the
        immediately preceding analyze(const int) method call is
        returned in the parameter.

        \note Not all ESTAnalyzer classes may compute additional
        alignment data.  In this case, this method will return \c
        false.  Furthermore, if a previous analyze() method call was
        not made, then the value returned in alignmentData parameter
        is not defined.
        
        \return This method returns \c true if the alignment data is
        actually computed by this ESTAnalyzer. The default
        implementation of this method always returns \c false.
    */
    virtual bool getAlignmentData(int& UNREFERENCED_PARAMETER(alignmentData))
    { return false; }

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

    /** Obtain a valid (or the best) metric generated by this
        analyzer.

        This method can be used to obtain a valid metric value for
        this analyzer.  This value can be used to initialize metric
        values. By default this method returns 0, which should be
        ideal for distance-based metrics.

        \note Dervied similarity-based metric classes must override this
        method to provide a suitable value.

        \return This method returns a valid (or the best) metric
        for this EST analyzer.
    */
    virtual float getValidMetric() const { return 0; }

    /** Determine preferred dummy EST lengths to be used with this
        analyzer.

        <p>This method can be used to determine the preferred dummy
        EST lengths to be used with this EST analyzer.  This method
        may be overridden in derived classes to provide a more
        appropriate dummy EST length.</p>

        <p>Dummy ESTs are used for the following purpose: When
        clustering FASTA data that contains low complexity reads, the
        low complexity reads provide false relationships between ESTs
        giving raise to very large clusters.  These large clusters are
        created because transitive relationships are established
        between ESTs due to low complexity reads.</p>

        <p>In order to avoid super-clusters that get formed due to low
        complexity reads, PEACE adds two dummy ESTs, one with all \c
        "AAAAA...." and another with all \c "CCCCCC...". The length of
        the ESTs must be appropriately chosen based on the type of
        analyzer used. This method helps ClusterMaker hierarchy to
        determine the appropriate dummy EST length.</p>
        
        \return The default implementation of this method always
        returns 128.
    */
    virtual int getPreferredDummyESTLength() const { return 128; }
    
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

    /** Method to attach a heuristic chain to this EST analyzer.

        \param[in] chain The heuristic chain to be attached.

        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns an error code.
    */
    virtual int setHeuristicChain(HeuristicChain* chain);

    /** Method to obtain the heuristic chain set for this EST
        analyzer.

        This method may be used to obtain a pointer to the heuristic
        chain set for use by this analyzer. If a heuristic chain has
        not been set, then this method returns NULL.

        \note The caller must \c not modify or delete the returned
        heuristic pointer.
        
        \return A pointer to the heuristic chain associated set for
        this analyzer.  If a heuristic has not been set, then this
        method returns NULL.
    */
    virtual HeuristicChain* getHeuristicChain() const { return chain; }
    
    /** Method to display performance statistics.

        This method can be used to display any statistics collated by
        this class (and its descendants) regarding their operation and
        performance.  This method was primarily introduced to enable
        derived classes a mechanism to override statistics display and
        print additional information.

        \note The default implementation in the base class prints
        statistics from the heuristic chain associated with this
        analyzer.
		
        \param[out] os The output stream to which the statistics must
        be written.
    */
    virtual void displayStats(std::ostream& os);

	/** Obtain the initial reference EST index.

		This method can be used to obtain the reference EST index that
		is specified by the user.  This value can be changed by the
		user via the \c --estIdx command-line parameter.
		
		\return Obtain the reference EST index specified by the user.
	*/
	inline int getInitialRefESTidx() const { return initialRefESTidx; }
	
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

        \param[in] analyzerName The human readable name for this EST
        analyzer.  This name is used when generating errors, warnings,
        and other output messages for this analyzer.
    */
    ESTAnalyzer(const std::string& analyzerName);

    /** Analyze and compute a similarity or distance metric between
        a given EST and the reference EST using the heavy weight metric
        associated with this ESTAnalyzer.

        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \note This method may return -1, if the otherEST is
        significantly different from the reference EST (possibly
        warranting no further analysis) that a meanigful metric cannot
        be generated.
	
        \param[in] otherEST Pointer to an immutable EST object with
        which the reference EST is to be compared.

        \return This method returns a similarity/distance metric by
        comparing the ESTs. This method may return -1, if the otherEST
        is significantly different from the reference EST (possibly
        warranting no further analysis) that a meanigful metric cannot
        be generated.
    */
    virtual float getMetric(const EST* otherEST) = 0;

    /** Pointer to the reference EST to compare with.

        This member object is used to hold a pointer to a given
        reference EST.  This member is initialized in the constructor
        and is changed by the setReferenceEST() id.
    */
    const EST* refEST;

    /** The heuristic chain associated with this EST analyzer.

        The heuristic chain contains a sequence of heuristics that
        must be used to minimize the number of pairs of ESTs that must
        be actually analyzed (using heavy weight algorithms such as
        D2).  The chain is created in the \c main method via a call to
        HeuristicChain::setupChain method and is set by \c main method
        via a call to setHeuristicChain method.
    */
    HeuristicChain* chain;

    /** A shortcut reference to the shared list of cDNA fragments
        being analyzed.

        This list holds a pointer to the shared list of cDNA fragments
        currently being analyzed.  This pointer is initialized to NULL
        in the constructor.  A valid pointer is filled in when the
        initialize() method is invoked.  The pointer to the ESTList is
        set via a call to setESTList() method in this class.
        Typically, this pointer is set by the
        ClusteringSubSystem::initializeSubComponents() method.
    */
    ESTList* estList;

	/** Command-line argument containing the initial reference EST
		index set by the user.

		This instance variable is used to track the initial reference
		EST index value specified by the user via the command-line
		argument \c --estIdx.  The default value for this instance
		variable is \c -1 (an invalid index).

		\note This value must not be changed by the derived classes.
	*/
	int initialRefESTidx;
	
private:
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
