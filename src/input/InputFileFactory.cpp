#ifndef INPUT_FILE_FACTORY_CPP
#define INPUT_FILE_FACTORY_CPP

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

#include "InputFileFactory.h"
#include "FASTAFile.h"
#include "SFFReader.h"
#include "PDBListFile.h"
#include "SubSystem.h"
#include "OnDemandESTList.h"
#include "RuntimeContext.h"

InputFileFactory::InputFileFactory() : Component("InputFileFactory") {
    // Initialize command line parameters
    noMaskBases     = false;
    randomizeNbases = false;
    noOnDemand      = false;
    noNormalizeNTs  = false;
}

InputFileFactory::~InputFileFactory() {
    // Empty constructor begets an empty destructor
}

InputFile*
InputFileFactory::create(const std::string& fileName) {
    if (fileName.size() < 1) {
        return NULL;
    }
    if (SFFReader::isSFF(fileName)) {
        // This is a valid SFF file. So process it as one.
        return create(fileName, "sff");
    }
    if (FASTAFile::isFASTA(fileName)) {
        // Appears to be a valid FASTA file.
        return create(fileName, "fasta");
    }
    // Could not auto detect.
    std::cerr << "Unable to auto detect file type for: " << fileName
              << std::endl;
    return NULL;
}

InputFile*
InputFileFactory::create(const std::string& fileName,
                         const std::string& format) {
    if (fileName.size() < 1) {
        return NULL;
    }
    InputFile* inputFile = NULL;
    
    if (format == "fasta") {
        inputFile = new FASTAFile(fileName);
    } else if (format == "sff") {
        inputFile = new SFFReader(fileName);
    } else if (format == "pdblist") {
        inputFile = new PDBListFile(fileName);
    }

    if ((inputFile == NULL) || (!inputFile->good())) {
        // Could not open file successfully.
        std::cerr << "Unable to open " << fileName << " as a "
                  << format << " file for reading."
                  << std::endl;
        if (inputFile != NULL) {
            delete inputFile;
            inputFile = NULL;
        }
    } else {
        // File was opened successfully. Setup options
        inputFile->setOptions(noMaskBases, randomizeNbases, !noNormalizeNTs);
    }
    // Return file pointer (if any) back to the caller
    return inputFile;
}

void
InputFileFactory::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"--no-mask-bases", "Don't mask out lower case nucleotides in reads",
         &noMaskBases, ArgParser::BOOLEAN},    
        {"--rnd-n-bases", "Change 'N' bases to a random ACGT base",
         &randomizeNbases, ArgParser::BOOLEAN},  
        {"--fastaFiles", "FASTA files (space separated) to be processed",
         &fastaFileNames, ArgParser::STRING_LIST},
        {"--sffFiles", "SFF files (space separated) to be processed",
         &fastaFileNames, ArgParser::STRING_LIST},
        {"--fastaFile", "Deprecated. Now it is an alias for --fastaFiles",
         &fastaFileNames, ArgParser::STRING_LIST},
        {"--sffFile", "Deprecated. Now it is an alias for --sffFiles",
         &sffFileNames, ArgParser::STRING_LIST},
        {"--estFile", "Deprecated. Now it is an alias for --estFiles",
         &estFileNames, ArgParser::STRING_LIST},
        {"--estFiles", "Space separated data files (format auto-detected) to be processed",
         &estFileNames, ArgParser::STRING_LIST},
        {"--pdbListFiles", "Set of PDB-list file path(s)",
         &pdbListNames, ArgParser::STRING_LIST},        
        {"--no-ondemand", "Hold all reads in memory (faster; uses more RAM)",
         &noOnDemand, ArgParser::BOOLEAN},
        {"--no-normalize", "Use sequences as it and don't normalize them to just nucleotide values (use for protiens or other inputs)",
         &noNormalizeNTs, ArgParser::BOOLEAN},        
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(Arguments);
}

bool
InputFileFactory::initialize() {
    // Load all fasta files.
    if (!loadFiles(fastaFileNames, "fasta")) {
        // Error occured during file loading.
        return false;
    }
    // Load all SFF files
    if (!loadFiles(sffFileNames, "sff")) {
        // Error occured during file loading.
        return false;
    }
    if (!loadFiles(pdbListNames, "pdblist")) {
        // Error occured during file loading.
        return false;
    }
    // Load all files via auto-detection
    if (!loadFiles(estFileNames, "")) {
        // Error occured during file loading.
        return false;
    }
    // All of the data was successfully loaded
    return true;
}

bool
InputFileFactory::loadFiles(ArgParser::StringList fileNames,
                            const std::string& format,
                            const long startIdx, const long endIndex) {
    ASSERT ( subSystem != NULL ); 
    RuntimeContext *context = subSystem->getContext();   
    ASSERT ( context != NULL );
    ASSERT ( context->getESTList() != NULL );
    OnDemandESTList* estList =
        dynamic_cast<OnDemandESTList*>(context->getESTList());
    if (estList == NULL) {
        std::cerr << "InputFileFactory::loadFiles() was expected to find a "
                  << "OnDemandESTList object in the context but did not "
                  << "find one. Exiting prematurely!" << std::endl;
        return false;
    }

    // Override for start index if noOnDemand flag is set and if the
    // startIndex is the default MAX_READS value.
    const int startIndex = (noOnDemand && (startIdx == MAX_READS) ? 0 :
                            startIdx);
    for(size_t i = 0; (i < fileNames.size()); i++) {
        // Open the data file specified by the user using the user
        // supplied format or through auto-detection.
        InputFile* inputFile = (format == "") ?
            create(fileNames[i]) : create(fileNames[i], format);
        if (inputFile == NULL) {
            // Error opening this file with specified format
            return false;
        }
        // Try loading fragments from this file into the est list
        if (!estList->add(inputFile, startIndex, endIndex)) {
            // Error loading fragments in to the EST list.
            return false;
        }
        // Data from this file was loaded successfully. Add meta data
        // about the input file to the runtime context.
        char timeStamp[128];
        getTimeStamp(fileNames[i].c_str(), timeStamp);
        std::string fileInfo = fileNames[i] + " (" + std::string(timeStamp)
            + ")";
        context->addConfig("Input file", fileInfo);
    }
    // All the files were successfully loaded
    return true;
}

void
InputFileFactory::finalize() {
    // Nothing to be done for now.
}

#endif
