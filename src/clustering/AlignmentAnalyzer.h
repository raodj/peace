#ifndef ALIGNMENT_ANALYZER_H
#define ALIGNMENT_ANALYZER_H

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

/** An analyzer to report similarity between 2 strains based on global
   alignment between a pair of strains.

   This is a more heavy-weight (that is, uses more CPU-time) analyzer
   that uses alignment information between a pair of reads to estimate
   similarity between strains.  This analyzer uses the NW score and
   normalizes the score from 0 (identical) to 1.0 (completely
   different) based on the aligned length.
*/
class AlignmentAnalyzer : public ESTAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~AlignmentAnalyzer();

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

    /** Method to initialize analyses.

        This method is currently just a placeholder for performing any
        initialization operations. Currently this method does not have
        any specific initialization operations to perform.

        \note As per the API requirements, this method calls the base
        class's initialize method to allow the base class to perform
        its own initialization operations as needed.
     */
    virtual bool initialize();

    /** Set the reference EST id for comparative analysis.
        
        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        just makes a note of the reference est for future calls to the
        analyze method.

        \param[in] est A pointer to an immutable EST object to be used
        as the reference.  The reference EST will be subsequently
        analyzed with a number of other ESTs via the analyze() method.
	
        \return Currently, this method always returns zero.
    */
    virtual int setReferenceEST(const EST* est) {
        this->refEST = est;
        return 0;
    }

    /** Method to perform exhaustive analysis to generate a detailed
        report. Currently, this method is not implemented for this
        class.

        \return This method currently simply returns 1 to indicate an
        error.
    */
    virtual int analyze() {
        return 1;   // error
    }
    
    /** Analyze and compute a similarity metric between a given EST
        and the reference EST using global alignment between the 2
        reads.

        This method is used to compare a given EST with the reference
        EST (set via the call to the setReferenceEST()) method.

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
    virtual float getMetric(const EST* otherEST);

    /** Determine if this EST analyzer provides distance metrics or
        similarity metrics.
        
        \return This method returns \c true to indicate
        that this EST analyzer operates using distance metrics.
    */
    virtual bool isDistanceMetric() const { return true; }

    /** Method to compare two metrics generated by this class.

        This method overrides the default implementation in the base
        class (that performs similarity-based comparison) to use
        distance based comparisons.

        \note As per the ESTAnalyzer API requirements, EST analyzers
        that are based on distance measures (such as this
        Alignment-based analyzer) \b must override this method.
        
        \param[in] metric1 The first metric to be compared against.

        \param[in] metric2 The second metric to be compared against.
        
        \return This method returns \c true if \c metric1 is
        comparatively better than \c metric2.
    */
    inline bool compareMetrics(const float metric1, const float metric2)
        const override
    { return (metric1 < metric2); }

    /** Obtain an invalid (or the worst) metric generated by this
        analyzer.

        This method can be used to obtain an invalid metric value for
        this analyzer.  This value can be used to initialize metric
        values.

        \note Dervied distance-based metric classes must override this
        method to provide a suitable value.

        \return This method returns an invalid (or the worst) metric
        for this EST analyzer.
    */
    virtual float getInvalidMetric() const override { return 1.f; }
    
protected:
    /** The default constructor.

       The default constructor for this class.  The constructor is
       made private so that this class cannot be directly
       instantiated.  However, since the ESTAnalyzerFactory is a
       friend of this class; therefore it can instantiate this2
       analyzer.  Accordingly, the ESTAnalyzerFactory::create() method
       must be used to instantiate this class.
    */
    AlignmentAnalyzer();

    inline int encodeBase(char c) const {
        return baseCodes.at(c);
    }
    
    void setScoringMatrix(int match, int mismatch);

    /**
       A convenient structure to wrap the results from the alignments
       produced by the alignment algorithm(s).
     */
    struct AlignResult {
	int score;
	std::string str1;
	std::string str2;
    };
    
    int getNWScore(const std::string& s1, const std::string& s2) const;
    AlignResult getNWAlignment(const std::string& s1, const std::string& s2) const;
    AlignResult getSWAlignment(const std::string& s1, const std::string& s2) const;    
    
private:
    /** Parameter to set the score associated with a match between two
        bases. This value can be modified via the --matchScore
        command-line argument.
    */
    int matchScore = 1;

    /** Parameter to set the score associated with a mismatch between two
        bases. This value can be modified via the --mismatchScore
        command-line argument.
    */
    int mismatchScore = -1;

    /** Parameter to set the score associated with a gap between two
        bases. This value can be modified via the --gapPenalty
        command-line argument.
    */
    int gapPenalty = -2;

    /**
       The scoring matrix used during alignment.
     */
    std::vector<std::vector<int>> scoreMatrix;

    /**
       This vector is initialized with positional codes for the
       various nucleotides and amino acids to make look-up faster.
    */
    std::vector<int> baseCodes;
};

#endif
