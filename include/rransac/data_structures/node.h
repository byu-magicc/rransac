#ifndef RRANSAC_DATA_STRUCTURES_NODE_H_
#define RRANSAC_DATA_STRUCTURES_NODE_H_

#include <vector>

#include "rransac/common/measurement/measurement_factory.h"

namespace rransac
{
/**
 * \class Node
 * The Node class contains pointers to children nodes or measurement objects.
 * This class is only used by the DataTree class data structure, where each node
 * of the tree is a Node Object.
 */ 
class Node
{
public:
    /**
     * Node Constructor.
     */
    Node();

    /**
     * Node Destructor. Delete all children node pointers to avoid memory leak.
     */
    ~Node();

    /**
     * Adds Node to vector of children nodes
     * @param[in] Pointer to new node
     */
    void AddChildNode(std::shared_ptr<Node> new_node);

    /**
     * Removes all children nodes, but doesn't delete them from memory.
     * This is handled in the destructor of the RRStarTree.
     */
    void RemoveAllChildNodes();
    
    /**
     * Adds a measurement to this node.
     * @param[in] Meas struct object to be added to node
     */
    void AddMeasurement(Meas measurement);
    
    /**
     * Removes all measurements from node.
     */
    void RemoveAllMeasurements();

    /**
     * Compresses or expands node_bounding_box_ depending on
     * whether or not a child node/meas object was added to
     * this node.
     */
    void UpdateNodeBoundingBox();

    /**
     * Sets the OverflowTreatment boolean variable.
     * @param[in/out] Sets overflow_treatment_was_called_
     */
    void SetOverflowTreatment(bool overflow_treatment_setter);

    /**
     * Returns vector of children nodes of this node.
     */
    std::vector<std::shared_ptr<Node>> GetNodeChildren();

    /**
     * Returns number of children nodes.
     */
    int GetNumberOfChildren();

    /**
     * Returns node_bounding_box_ around node.
     */
    Eigen::MatrixXd GetNodeBounds();

    /**
     * Returns all measurements contained in node.
     */
    std::vector<Meas> GetMeasurements();

    /**
     * Returns number of measurement objects in node.
     */
    int GetNumberOfMeasurements();

    /**
     * Returns whether or not overflow treatment has 
     * been called on this node from RRStarTree class.
     */
    bool WasOverflowTreatmentCalled();

    /**
     * Returns whether or not this node is a leaf node.
     */
    bool IsNodeALeafNode();

    /**
     * Returns center coordinates of node_bounding_box_.
     */
    Eigen::MatrixXd GetCenterOfBoundingBox();
    
private:

    std::vector<std::shared_ptr<Node>> children_nodes_; /**< Vector filled with pointers to children Nodes of this Node. */
    
    Eigen::MatrixXd node_bounding_box_; /** < Matrix containing outer bounds of each dimension.
                                            Formatted such that [first_dim_min_val, first_dim_max_val;
                                                                second_dim_min_val, second_dim_max_val;
                                                                        ...                 ...        ],
                                            are the meaning of the values at different positions within matrix.*/
    
    std::vector<Meas> measurements_; /**< Vector filled with Meas children objects. */
    
    bool overflow_treatment_was_called_; /**<   A flag that indicates if Node has been used in the OverflowTreatment
                                                function within the RRStarTree class. */
    
    bool is_leaf_node_; /**< A flag that indicates if the Node is a leaf Node. */
};

} // namespace rransac

#endif // RRANSAC_DATA_STRUCTURES_NODE_H_