#ifndef MATRIX_FILE_ANALYZER_H
#define MATRIX_FILE_ANALYZER_H

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

/** MatrixFileAnalyzer: EST Analyzer that simply obtains distances
    from a matrix file for processing.

    <p>This analyzer provides a simple interface for using precomputed
    distance/similarity values from a given data file.  A matrix data
    file must have the following format:

    \code

    # Lines starting with '#' character are assumed to be comments and
    # they are ignored. The first non-comment line must be a number
    # indicating the number of ESTs for which data is present in the
    # file. For example here is a file for 3 ESTs
    3
    
    # After the number of EST's there must be nxn matrix of values
    # where n is number of EST's.
    0.0 10   20
    10  0    15
    20  15.2 0

    \endcode

*/
class MatrixFileAnalyzer : public ESTAnalyzer {
    friend class ESTAnalyzerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~MatrixFileAnalyzer();
    
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

        This method is invoked just before commencement of EST
        analysis.  This method loads the list of distance values from
        the given input file and pouplates the matrix
        MatrixFileAnalyzer::distanceValues for futher use in the
        analyze() method.

        \return If the matrix data was successfully loaded from the
        data file then this method returns \c true.  Otherwise this
        method reports an error to std::cerr and returns \c false.
    */
    virtual bool initialize();

    /** Set the reference EST  for analysis.

        This method is invoked just before a batch of ESTs are
        analyzed via a call to the analyze(EST *) method.  This method
        currently saves the index in the instance variable for further
        look up.

        \note This method must be called only after the initialize()
        method is called.

		\param[in] est A pointer to an immutable EST object to be used
		as the reference.  The reference EST will be subsequently
		analyzed with a number of other ESTs via the analyze() method.
		
        \return This method returns \c true if the estIdx was within
        the range of values that were loaded from the data file.
        Otherwise this method returns 1 as the error code.
    */
    virtual int setReferenceEST(const EST* est);
    
    /** Method to pefrom a batch of EST analysis.
        
        The \c ESTAnalyzer base class requires this method to be
        overloaded in the dervied class(es).  This method is used to
        perform the core tasks of EST analysis for teh
        MatrixFileAnalyzer.  This method operates in the following
        manner:

        <ol>

        <li>First it loads the necessary distance information from the
        supplied data file file using the initialize() method.  If the
        data is not successfully loaded then this method returns right
        away with 1.<li>

        <li>Upon successfully loading the data, the reference EST is
        set via the setReferenceEST() method.  If the reference EST is
        not correctly determined, then this method immediately returns
        with 2.</li>

        <li>For each EST in the list of ESTs it performs the following
        tasks:

        <li>It logs the similarity metric using suitable methods in
        the ESTAnalyzer base class.<li>

        <li>If all the processing proceeds successfully, this method
        returns 0 (zero).
        
        </ol>
    */
    virtual int analyze();
    
protected:
    /** Analyze and obtain a distance (or similarity) metric.
        
        This method can be used to compare a given EST with the
        reference EST (set via the call to the setReferenceEST())
        method.

        \param[in] otherEST Pointer to an immutable EST object with
        which the reference EST is to be compared.

        \return This method returns the distance (or similarity as the
        case may be) value loaded by the data file (loading is
        performed by the initialize() method). 
    */
    virtual float getMetric(const EST* otherEST);

    /** The number of EST's for which we have data in \c distanceValues
        matrix.

        This instance variable's value is set by the initialize method
        in this class.
    */
    int estCount;
    
    /** The array of distance values.

        This matrix contains the distance values between a given pair
        of ESTs.  The zero-based index of the EST is used to look up
        values in this matrix.  For example, given a pair of ESTs
        <est1, est2> \c distanceValues[est1][est2] provides the
        distance from est1 to est2 while \c distanceValues[est2][est1]
        provides the distance from est2 to est1.  Note that distances
        do not have to be symmetric.
    */
    float **distanceValues;

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
        comparatively better than \c metric2.
    */
    bool compareMetrics(const float metric1, const float metric2) const
    { return (metric1 < metric2); }
    
    /** Obtain an invalid (or the worst) metric generated by this
        analyzer.
        
        This method can be used to obtain an invalid metric value for
        this analyzer.  This value can be used to initialize metric
        values.
        
        \note Dervied distance-based metric classes must override this
        method to provide a suitable value.
        
        \return This method returns an invalid (or the worst) metric
        of 1e7 for this EST analyzer.
    */
    float getInvalidMetric() const { return 1e7; }

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
    /*
        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ESTAnalyzerFactory is a
        friend of this class, an object can be instantiated via teh
        ESTAnalyzerFactory::create() method.
    */
    MatrixFileAnalyzer();

    /** Method to read a line from a given EST file.

        This method is a helper method that is used to load a given
        line from the file.

        \param[in,out] fp The file pointer from where the data is to be
        read.

        \return The line read from the file.  If EOF was reached then
        this method returns an empty line.
    */
    std::string readLine(FILE *fp);

    /** Utility method to read metrics from a given line into a given
        array from a given starting position.

        This method is a utility method that is used to read parse in
        a line containing space separated set of values into the array
        values. Data is stored into the values array starting with
        startPos.

        \param[in] line The line whose contents is to be updated.

        \param[out] values The array into which the values must be
        stored.

        \param[in] startPos The starting position in the array from
        where values must be stored. 

        \param[in] maxValues The maximum number of items that must be
        processed from line.
        
        \return This method returns the number of values
        actually processed and stored into the values array.
    */
    int parseMetrics(const char* line, float *values,
                     const int startPos, const int maxValues);

    /** Helper method to process EST count information from a line.

        This method is a helper method that is used to read the number
        of EST's from the file and create the distanceValues matrix.

        \param[in] line The line from where the EST count information
        is to be read.

        \return This method returns \c true if the EST count was
        processed successfully.  Otherwise this method returns \c
        false.
    */
    bool parseESTCount(const char *line);
    
    /** The matrix data file from where distance metrics are to be
        read.

        This member object is used to hold the file name from where
        all the distance matrix data is to be loaded.  The value is
        set to the value supplied by the user via a suitable command
        line argument by the parseArguments() method.
    */
    std::string dataFileName;
};

#endif
