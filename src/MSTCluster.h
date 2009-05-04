#ifndef MST_CLUSTER_H
#define MST_CLUSTER_H

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

#include "MST.h"

class MSTCluster;

/** \def ClusterList

    \brief Typedef for std::vector<MSTCluster>.

    This typedef is a convenience definition to handle a vector of
    Cluster objects associated with a given MSTCluster.
*/
typedef std::vector<MSTCluster*> ClusterList;


/** Class to encapsulate cluster information and aid in building
    clusters using Minimum Spanning Tree (MST) data.

    This class was introduced to encapsulate and manage the data
    associated with various clusters built using MST information.
*/
class MSTCluster {
    // Insertion operator to make dumping MST for debugging easier.
    friend std::ostream& operator<<(std::ostream&, const MSTCluster&);  
public:
    MSTCluster(MSTCluster* owner = NULL);
    ~MSTCluster();

    double makeClusters(NodeList& nodeList, const double percentile);
    void add(const MSTNode& node);

    void printClusterTree(std::ostream& os = std::cout,
                          const std::string& prefix = "") const;
               
protected:
    double calculateThreshold(const NodeList& nodeList,
                              const double percentile) const;
    
private:
    ClusterList clusterList;
    const MSTCluster* parent;
    NodeList members;
};

/** \func operator<<

    Insertion operator to stream MSTCluster information to a given
    output stream.  This method provides a convenient mechanism to
    dump the complete MSTCluster information for debugging purposes.
    The clusters in the MSTCluster are displayed in the order in which
    they were created.
*/
extern std::ostream& operator<<(std::ostream&, const MSTCluster&);

#endif
