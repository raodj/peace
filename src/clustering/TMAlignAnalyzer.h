#ifndef TMALIGN_ANALYZER_H
#define TMALIGN_ANALYZER_H

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

/** An alias to make this source contextually meaningful */
using PDB = EST;

/** An analyzer to report similarity between 2 strains based on
    structural similarity.

   This is a more heavy-weight (that is, uses more CPU-time) analyzer
   that uses structural similarity score reported by TM-Align as the
   distance between a pair of reads.  TM-Align is a 3rd-party program
   that uses structural information from CIF files to report
   structural similarity. TM-Align scores are in the range 0.0
   (dissimilar) to 1.0 (identical).

   \see TMAlign.cpp
*/
class TMAlignAnalyzer : public ESTAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~TMAlignAnalyzer();

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

    /** Set the reference EST containing the CIF file data for
        comparative analysis.
        
        This method is invoked just before a round of comparisons are
        performed via a call to the analyze(EST *) method.  This
        method just makes a note of the reference structure for future
        calls to the analyze method.

        \param[in] est A pointer to an immutable EST object to be used
        as the reference.  The reference will be subsequently analyzed
        with a number of other ESTs via the analyze() method.
	
        \return Currently, this method always returns zero.
    */
    virtual int setReferenceEST(const PDB* pdb) {
        this->refEST = pdb;
        return 0;
    }

    /** Method to perform exhaustive analysis to generate a detailed
        report. Currently, this method is not implemented for this
        class.

        \return This method currently simply returns -1 to indicate an
        error.
    */
    virtual int analyze() {
        return 1;   // error
    }
    
    /** Analyze and compute the structural distance metric between two
        PDB/CIF files using TM-align.
        
        This method is used to compare a given PDB/CIF files (specifed
        by the information of the EST) against the PDB/CIF file set
        via the call to the setReferenceEST() method.
	
        \param[in] otherPDB Pointer to an immutable object that
        contains the path to the PDB/CIF file to be used for
        comparison.

        \return This method returns a distance metric computed as 1 -
         TMalign_score. The value returned by this method is 0
         (identical) to 1.0 (dissimilar).
    */
    virtual float getMetric(const PDB* otherPDB);

    /** Determine if this analyzer provides distance metrics or
        similarity metrics.
        
        \return This method returns \c true to indicate that this PDB
        analyzer operates using distance metrics.
    */
    virtual bool isDistanceMetric() const { return true; }

    /** Method to compare two metrics generated by this class.

        This method overrides the default implementation in the base
        class (that performs similarity-based comparison) to use
        distance based comparisons.

        \note As per the API requirements, analyzers that are based on
        distance measures (such as this one) \b must override this
        method.
        
        \param[in] metric1 The first metric to be compared against.

        \param[in] metric2 The second metric to be compared against.
        
        \return This method returns \c true if \c metric1 is
        comparatively better (i.e., lower distance) than \c metric2.
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
        for this analyzer.
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
    TMAlignAnalyzer();
    
private:
    // Possibly expose the pertinent TM-algin parameters here instead
    // of just using the default settings in TM-align.
};

#endif
