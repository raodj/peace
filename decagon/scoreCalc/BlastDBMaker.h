#ifndef BLAST_DB_MAKER_H
#define BLAST_DB_MAKER_H

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

/** \file BlastDBMaker.h

	\brief A simple class to check and build BLAST database for source
	transcripts.

	This file contains the declaration for the BlastDBMaker class.  It
	is used to run BLAST to build a database (used by MegaBlast).  The
	database is built only if needed.
*/

#include "BlastCommon.h"

/** Helper class to build a BLAST data base for a given source
    transcript file.
    
	<p>This is a helper class that is used to run \c formatdb command
	from BLAST to build a BLAST data base for a given source
	transcript file. The source transcript file is typically a file
	from which cDNA reads have been synthetically created (using a
	Java program called MetaSim). The BLAST database is required to
	run Mega Blast.</p>

    <p>Note that the BLAST database for a given source transcript file
    can be rused and does not need to be created for each
    run. Consequently, this class checks to see if the necessary BLAST
    database files (files with extensions \c _blast_db.nhr, \c
    _blast_db.nin, and \c _blast_db.nsq) exist for a given transcript
    file. If the files exist then the database creation step is
    short-circuited.</p>
*/
class BlastDBMaker : public BlastCommon {
public:
    /** The default constructor.

        The constructor does not perform any special operation other
        than to initialize members to default initial value.
    */
    BlastDBMaker();

    /** The destructor.

        The destructor has no specific tasks. It is present merely to
        adhere to coding conventions.
    */
    ~BlastDBMaker();

    /** Primary method to check and run \c formatdb.

        This method uses the supplied source transcript file name
        (parameter srcTransFile) without file type extension (such as:
        \c .fa or \c .fasta) as the base path to check if the files
        with the following names exist and are newer than the source
        transcript file:

        <ul>
        <li>\i srcTransFileNoExt \c _blast_db.nht</li>
        <li>\i srcTransFileNoExt \c _blast_db.nin</li>
        <li>\i srcTransFileNoExt \c _blast_db.nsq</li>        
        </ul>

        If the files exist and are newer than the source file then
        this method assumes they are up to date and returns
        successfully.

        <p>Otherwise this method uses BlastCommon::getBlastDBCmdLine()
        method to obtain the command to be run to generate the BLAST
        data base. The BLAST data base is written in the same
        directory as the source file.</p>

        \param[in] srcTransFile The full path to the source transcript
        file for which a BLAST data base is to be created (as
        needed).

		\return This method returns the path to the BLAST DB that was
		created and/or verified by this method.
		
        \exception DecagonException This method throws a decagon
        exception if errors occur when performing the various tasks.
    */
	std::string checkMakeDB(const std::string& srcTransFile)
        throw(DecagonException);
                 
protected:
    /** Helper method to obtain last modification time for a given
        file.
        
        This method is called from the checkMakeDB() to verify a file
        exists and obtain its last modification time.

        \param[in] fileName The path to the file whose last
        modification time is to be retrieved.

        \return This method returns the last modification time (in
        milliseconds since epoch) for the file, assuming the file
        exists.

        \exception DecagonException This method throws a decagon
        exception if the file does not exist.
    */
    time_t getTimestamp(const std::string& fileName) const
        throw(DecagonException);

    /** Helper method to maintain the base file name to be used for
        the BLAST data base.

        This method is a helper method that is invoked once from the
        checkMakeDB() method.  It is used to obtain the base path to
        be used for the BLAST DB.  This method essentially removes any
        trailing extensions by looking for a period (\c '.') from the
        end of the file and removing all characters (including the
        period) after it.

        \param[in] srcTransFile The supplied path to the source
        transcript file to be used to generate the BLAST data base
        name.

        \return The base file name to be used for the BLAST data base.
    */
    std::string getDBName(const std::string& srcTransFile) const;
    
private:
    // Currently this class does not have any private members.
};

#endif
