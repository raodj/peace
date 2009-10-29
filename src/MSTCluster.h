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
#include "ESTAnalyzer.h"

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

    double makeClusters(NodeList& nodeList, const double percentile,
			const int analysisCount, const ESTAnalyzer* analyzer);
    void add(const MSTNode& node);

    void printClusterTree(std::ostream& os = std::cout,
                          const std::string& prefix = "") const;

    void guiPrintClusterTree(std::ostream& os = std::cout,
							 const char *srcFile = NULL) const;

    void makeMergedClusters(const int size, int* parent, bool* root);
               
    //protected:
    double calculateThreshold(const int nodeCount,
                              const double percentile,
			      const int analysisCount,
			      const ESTAnalyzer* analyzer) const;

    inline int getClusterID() const { return clusterID; }
    
protected:
	void guiPrintTree(std::ostream& os) const;
	
private:
    ClusterList clusterList;
    const MSTCluster* parent;
    NodeList members;
    const int clusterID;
    
	/** Instance variable to track the next available cluster ID.

		This instance variable is used to generate unique cluster ID
		values for each newly created cluster.  It is intialized to
		zero. Each time a cluster is instantiated, the constructor
		uses this value to set the clusterID and then increments this
		value.
	*/
	static int clusterIDSequence;
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
