#ifndef CLUSTER_WRITER_CPP
#define CLUSTER_WRITER_CPP

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

#include "ClusterWriter.h"
#include "ClusterMaker.h"
#include "ArgParser.h"
#include "Utilities.h"
#include "MSTCluster.h"
#include "RuntimeContext.h"
#include <fstream>

ClusterWriter::ClusterWriter() : Component("ClusterWriter") {
    // Initialize command-line arguments to default value
    prettyPrint = false;
    guiPrint    = false;
    flatPrint   = false;
}


ClusterWriter::~ClusterWriter() {
    // Nothing to be done for now.
}

void
ClusterWriter::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"--pretty-print", "Print a pretty cluster tree.",
         &prettyPrint, ArgParser::BOOLEAN},
        {"--gui-print", "Print the cluster tree for GUI processing.",
         &guiPrint, ArgParser::BOOLEAN},
        {"--flat-print", "Print the cluster tree in verbose flat format.",
         &flatPrint, ArgParser::BOOLEAN},        
        {"--output-cls-file", "File to which clustering output must be written",
         &clusterFileName, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add arguments for this component to the parser
    argParser.addValidArguments(Arguments);
}

bool
ClusterWriter::write(RuntimeContext* context) {
    const ClusterMaker *clusterMaker = context->getClusterMaker();
    if ((clusterMaker == NULL) || (clusterMaker->getClusters() == NULL) ||
        (clusterMaker->getClusters()->isEmpty())) {
        // Nothing further to be done.
        return true;
    }
    ASSERT( context->getESTList() != NULL );
    const ESTList* estList = context->getESTList();
    // Redirect cluster output to outputFile as needed.
    std::ofstream outFile;
    if (!clusterFileName.empty()) {
        outFile.open(clusterFileName.c_str());
        if (!outFile.good()) {
            std::cerr << "Error opening output file " << clusterFileName
                      << " for writing cluster data.\n"
                      << "Cluster data will be dumped on stdout.\n";
        }
    }
    // Choose output stream depending on condition of outFile.
    std::ostream& dest = (outFile.is_open() && outFile.good()) ?
        outFile : std::cout;
    // Obtain the root of the cluster tree
    const MSTCluster* root = clusterMaker->getClusters();    
    if (prettyPrint) {
        root->printClusterTree(*estList, dest);
    } else if (guiPrint) {
        dest << "# Cluster Data\n";
        context->printConfig(dest);
        dest << "# Data format: <E,estIdx,clstrId> | "
             << "<C,clstrID,ParntClstrID>\n";
        root->guiPrintClusterTree(dest);
    } else if (flatPrint) {
        // Flat print the cluster information.
        root->printFlatClusterTree(*estList, dest);
    } else {
        // No pretty printing. Just dump the info out.
        dest << *root << std::endl;
    }
    // Return status depending on status of out file
    return outFile.good();
}

#endif
