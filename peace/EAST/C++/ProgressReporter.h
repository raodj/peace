#ifndef PROGRESS_REPORTER_H
#define PROGRESS_REPORTER_H

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

#include <fstream>
#include <string>

/** A process-wide unique singleton class used to report progress
   statistics.

   <p>This class provides a convenient mechanism for various
   components in EAST to peridocially report some progress information
   that can be utilized by the PEACE-GUI to provide feedback to the
   user.  The progress data is written to a file. The file name (where
   progress is to be written) is specified by the user as a
   command-line argument ("--progress prog.dat").  The file name is
   set in this class by the main() method.</p>

   <p>Currently, the progress information reported by EAST is a bit
   coarse. It reports progress as major steps of six-tuple generation,
   contig-formation, and reconstruction-are done. Later on, hopefully
   these steps can be somehow be further refined (at a later date) to
   provide a more meanigful progress information to the user.</p>

   \note This class is implemented as a singleton pattern. It is not
   meant to be instantiated. Instead, use the static get() method in
   this class to obtain a reference to the lgobal progress reporter
   object.
*/
class ProgressReporter {
public:
    /** Obtain a reference to the progress reporter object.

        This method must be used to obtain an instance of the
        process-wide unique instance of the progress reporter object.

        \return A reference to the progress reporter object that must
        be used to report progress.
    */
    static ProgressReporter& get() {
        return uniquePR;
    }

    /** Setup the destination file and maximum number of steps.

        This method must be used to setup the destination file to
        where the progress information is to be written.  Typically,
        this method is invoked from the main() method just before
        assembly commences after command-line arguments have been
        processed.  This method performs the following tasks:

        <ol>

        <li>It opens the progress file for writing. If errors occur it
        reports the error and aborts the process.</li>

        <li>On successfully opening the file, it reports an initial
        zero progress (via the reportProgress() method) </li>

        </ol>

        \param[in] outFile The name of the output file to which
        progress information must be written.

        \param[in] maxSteps The maximum number of steps that would be
        present in the progress information.
    */
    void initialize(const std::string& outFile, const int maxSteps);
    
    /** Method to report updated progress information.

        This method must be used to report updates to the progress
        information.  This method is frequently (but not so frequently
        that most of the time is spent in reporting progress) invoked.

        \param[in] progress The progress that has been made thus
        far. This value must be in the range 0 &lte; maxSteps.
    */
    void reportProgress(const int progress);
    
protected:
    /** The maximum number of steps of progress that will be reported.

        This value is an initial maximum progress value that is set
        just before EAST starts assembly.  This value is set via the
        initialize() method and is typically never changed during the
        lifetime of the progress.
     */
    int maxSteps;

    /** The output file to which progress information is written.

        This file stream is created and setup in the initialize()
        method. This stream is used periodically report progress
        updates (via the reportProgress() method) to a given progress
        file.
    */
    std::ofstream progFile;
    
private:
    /**
       The default and only constructor for this class.

       The constructor has been made private to ensure that this class
       is never directly instantiated.  Instead the get() method in
       this class must be used to obtain a reference to the globally
       unique instance.
    */
    ProgressReporter() {}

    /** The destructor.

        The destructor is called just before the EAST job
        completes. The destructor closes the handle to the progress
        file (if one has been opened).  The destructor has been made
        private to ensure that the globally-unique object is not
        deleted (even by accident).
    */
    ~ProgressReporter();

    /**
       The globally (process-wide) unique singelton instance of the
       progress reporter object that must be used to report progress
       information.
    */
    static ProgressReporter uniquePR;
};

#endif
