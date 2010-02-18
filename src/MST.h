#ifndef MST_H
#define MST_H

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
// Authors:   Dhananjai M. Rao              raodm@muohio.edu
//
//---------------------------------------------------------------------

#include "MSTNode.h"
#include <vector>

/** \def NodeList

    \brief Typedef for std::vector<MSTNode>.

    This typedef is a convenience definition to handle a vector of
    MSTNode objects associated with a given MST.
*/
typedef std::vector<MSTNode> NodeList;

/** \brief Class to represent a Minimum Spanning Tree (MST).

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
    /** The constructor.

        The default constructor for the MST.  The constructor performs
        no special operations and is present merely to adhere to
        coding conventions.

        \param[in] maxNodes The number of nodes for which space must
        be initially reserved.  Reserving space reduces the number of
        times the MST has to reallocate memory as nodes are added to
        it.

		\param[in] haveAlignmentMetric If this value is \c true, then
		that indicates that this MST also includes additional
		alignment information with each node in the MST.
    */
    MST(const int maxNodes, const bool haveAlignmentMetric);

    /** The destructor.

        The destructor performs no special tasks as the base class
        (std::vector) performs all the necessary cleanups.
    */
    ~MST() {}

    /** Determine if nodes in this MST have alignment metric values.

        This method can be used to determine if the nodes in this MST
        have alignment information saved for each node.
        
        \return This method returns \c true if the nodes in this MST
        have alignment metric.
    */
    inline bool hasAlignmentMetric() const { return haveAlignmentMetric; }
    
    /** Convenience method to add a MSTNode to the EST.

        This method may be used to append a node to the EST.  The
        first node to be added must be root node (with parentIndex set
        to -1).  This method performs no special checks to verify the
        integrity of the data.
        
        \param[in] parentIdx The index of the parent EST for this
        Node.  If a parent does not exist, then this value should be
        -1.
        
        \param[in] estIdx The index of the EST.  This value must be
        the index of the corresponding EST in the list of ESTs
        returned by EST::getESTList() method.
        
        \param[in] similarity The similarity metric between this EST
        and its parent EST.
        
        \param[in] alignmentInfo Additional alignment information
        between the parentIdx and the estIdx to be stored in this
        node.
    */
    void addNode(const int parentIdx, const int estIdx,
                 const float similarity, const int alignmentInfo);

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
        
	\param[in] srcFile The source file from where the ESTs were
	read. This file name is simply used for cross reference in the
	data generated for this MST.

	\param[in] threshold The threshold value to be used in clustering
	of the MST data.  Important for downstream analysis such as that
	done by assembly tools.
    */
    void serialize(const char *fileName, const char* srcFile,
		   const float threshold) const;
    
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

    /** Obtain the total distance of this MST.

        This method can be used to determine the total distance of the
        MST. This method iterates over all the nodes present in the
        MST and adds the distance/similarity metric together.  It then
        returns the total.

        \return The total MST distance for this MST.  If the MST has
        one (or zero) nodes then this method returns 0 (zero).
    */
    float getMSTDistance() const;
    
private:
    /** The list of nodes in this MST.

        This vector is used to track the set of nodes added to this
        minimum spanning tree.  The maximum nodes in the tree is
        typically set in the constructor.  Nodes are added via the
        addNode() method in this class.
    */
    NodeList nodeList;

    /** Instance variable to indicate if alignment metric in nodes is
	valid.

	This instance variable is initialized when this MST is
	created.  If this instance variable is \c true then each Node
	in this MST has a valid alignment information.
    */
    bool haveAlignmentMetric;
};

/** \fn std::ostream& operator<<(std::ostream&, const MST&)

    Insertion operator to stream MST information to a given output
    stream.  This method provides a convenient mechanism to dump the
    complete MST information for debugging purposes.  The nodes in the
    MST are displayed in the order in which they were inserted into
    the MST.
*/
extern std::ostream& operator<<(std::ostream&, const MST&);

#endif
