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
 * of the tree is a Node object.
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
     * @param[in] new_node pointer to new Node
     */
    void AddChildNode(std::shared_ptr<Node> new_node);

    /**
     * Removes all children nodes, but doesn't delete them from memory.
     * This is handled in the destructor of the RRStarTree.
     */
    void RemoveAllChildNodes();
    
    /**
     * Adds a measurement to this Node.
     * @param[in] Meas measurement struct object to be added to Node
     */
    void AddMeasurement(Meas measurement);
    
    /**
     * Removes all measurements from Node.
     */
    void RemoveAllMeasurements();

    /**
     * Compresses or expands node_bounding_box_ depending on
     * whether or not a child Node/Meas object was added to
     * this Node.
     */
    void UpdateNodeBoundingBox();

    /**
     * Sets the OverflowTreatment boolean variable.
     * @param[in/out] sets overflow_treatment_was_called_ variable to true or false
     */
    void SetOverflowTreatment(bool overflow_treatment_setter);

    /**
     * \return Returns vector of children nodes of this Node
     */
    std::vector<std::shared_ptr<Node>> GetNodeChildren();

    /**
     * \return Returns number of children nodes.
     */
    int GetNumberOfChildren();

    /**
     * \return Returns node_bounding_box_ a round Node.
     */
    Eigen::MatrixXd GetNodeBounds();

    /**
     * \return Returns all measurements contained in Node.
     */
    std::vector<Meas> GetMeasurements();

    /**
     * \return Returns number of measurement objects in Node.
     */
    int GetNumberOfMeasurements();

    /**
     * \return Returns whether or not the OverflowTreatment function from the RStarTree class was called on this Node
     */
    bool WasOverflowTreatmentCalled();

    /**
     * \return Returns whether or not this Node is a leaf Node.
     */
    bool IsNodeALeafNode();

    /**
     * \return Returns center coordinates of node_bounding_box_.
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