#ifndef MST_H
#define MST_H

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

#include "MSTNode.h"
#include <vector>

/** \def NodeList

    \brief Typedef for std::vector<MSTNode>.

    This typedef is a convenience definition to handle a vector of
    MSTNode objects associated with a given MST.
*/
typedef std::vector<MSTNode> NodeList;

/** Class to represent a Minimum Spanning Tree (MST).

    This class is used to represent a MST that is built by the
    MSTClusterMaker.  It essentially maintains a vector of MSTNode
    objects.  The order of nodes in the MST reflects the order in
    which nodes were added to create the tree.  Refer to the
    documentation on the MSTNode class for details on the data
    encapsulated by it.
*/
class MST {
    // Insertion operator to make dumping MST for debugging easier.
    friend std::ostream& operator<<(std::ostream&, const MST&);    
public:
    /** The default constructor.

        The default constructor for the MST.  The constructor performs
        no special operations and is present merely to adhere to
        coding conventions.
    */
    MST(const int maxNodes);

    /** The destructor.

        The destructor performs no special tasks as the base class
        (std::vector) performs all the necessary cleanups.
    */
    ~MST() {}

    /** Convenience method to add a MSTNode to the EST.

        This method may be used to append a node to the EST.  The
        first node to be added must be root node (with parentIndex set
        to -1).  This method performs no special checks to verify the
        integrity of the data.
        
        \param[in] parentIndex The index of the parent EST for
        this Node.  If a parent does not exist, then this value
        should be -1.
        
        \param[in] estIndex The index of the EST.  This value must
        be the index of the corresponding EST in the list of ESTs
        returned by EST::getESTList() method.
        
        \param[in] similarityMetric The similarity metric between
        this EST and its parent EST.
    */
    void addNode(const int parentIdx, const int estIdx,
                 const float similarity);

    /** Obtain a mutable list of nodes in this MST.

	This method can be used to obtain a mutable list of nodes in
	contained by this MST.

	\return A reference to the list of nodes encapsulated by this
	class.
    */
    NodeList& getNodes() { return nodeList; }

    /** Helper method to write MST data to a given file.

        This method is a helper method that is used to write MST data
        to a given output file.  This method is typically called from
        the makeClusters() method after MST construction has been
        completed.
        
        \param[in] fileName The file to which the MST data is to be
        written.
        
	\param[in] srcFile The source file from where data would be
	written.
    */
    void serialize(const char *fileName, const char* srcFile) const;
    
    /** Method to create and load data for MST from a text file.
        
        This method creates an MST and populates it with node data
        from a given text file.  The text file must have been created
        through an earlier call to the serialize() method.  If the MST
        was successfully populated then this method returns a pointer
        to the newly created MST.  On errors this method returns NULL.

        \param[in] fileName The text file from which the data for the
        MST must be read.

        \return On success this method returns a valid pointer to the
        newly created (and populated) MST.  On errors it returns NULL.
    */
    static MST* deSerialize(const char *fileName);
    
private:
    /** The list of nodes in this MST.

        This vector is used to track the set of nodes added to this
        minimum spanning tree.  The maximum nodes in the tree is
        typically set in the constructor.  Nodes are added via the
        addNode() method in this class.
    */
    NodeList nodeList;
};

/** \func operator<<

    Insertion operator to stream MST information to a given output
    stream.  This method provides a convenient mechanism to dump the
    complete MST information for debugging purposes.  The nodes in the
    MST are displayed in the order in which they were inserted into
    the MST.
*/
extern std::ostream& operator<<(std::ostream&, const MST&);

#endif
