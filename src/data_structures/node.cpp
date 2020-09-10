#include "rransac/data_structures/node.h"

namespace rransac
{

//-----------------------------------------------------------------------------

Node::Node() = default;

//-----------------------------------------------------------------------------

Node::~Node()
{
    for (auto i : children_nodes_) {
        delete i;
    }
    children_nodes_.clear();
}

//-----------------------------------------------------------------------------

void Node::AddChildNode(Node* new_node)
{
    children_nodes_.push_back();
}

//-----------------------------------------------------------------------------

void Node::RemoveAllChildNodes() 
{
    children_nodes_.clear();
}

//-----------------------------------------------------------------------------

void Node::AddMeasurement(Meas measurement)
{
    measurements_.push_back(measurement);
}

//-----------------------------------------------------------------------------

void Node::RemoveAllMeasurements()
{
    measurements_.clear();
}

//-----------------------------------------------------------------------------

void Node::UpdateNodeBoundingBox()
{
    std::vector<double> mins, maxs;

    if (IsNodeALeafNode() == true) {
        //For each data point
        for (int i = 0; i < measurements_.size(); ++i) {
            Eigen::MatrixXd current_measurement_data = measurements_.at(i).data;

            //For each dimension
            for (int dimension = 0; dimension < current_measurement_data.rows(); ++dimension) {
                if (i == 0) {
                    mins.push_back(current_measurement_data(dimension,0));
                    maxs.push_back(current_measurement_data(dimension,0));
                } else {
                    if (current_measurement_data(dimension,0) < mins.at(dimension)) {
                    mins.at(dimension) = current_measurement_data(dimension,0);
                    }
                    if (current_measurement_data(dimension,1) > maxs.at(dimension)) {
                    maxs.at(dimension) = current_measurement_data(dimension,0);
                    }
                }
            }
        }
    } else {
        //For each child Node
        for (int i = 0; i < children_nodes_.size(); ++i) {
          Eigen::MatrixXd current_node_bounding_box = children_nodes_.at(i)->GetNodeBounds();

          if (!currNodeBoundingBox.empty()) {
            //For each dimension
            for (int dimension = 0; dimension < current_node_bounding_box.rows(); ++dimension) {
                if (i == 0) {
                    mins.push_back(current_node_bounding_box(dimension,0));
                    maxs.push_back(current_node_bounding_box(dimension,1));
                } else {
                    if (current_node_bounding_box(dimension,0) < mins.at(dimension)) {
                    mins.at(dimension) = current_node_bounding_box(dimension,0);
                    }
                    if (current_node_bounding_box(dimension,1) > maxs.at(dimension)) {
                    maxs.at(dimension) = current_node_bounding_box(dimension,1);
                    }
                }
            }
          }
        }
    }

    Eigen::MatrixXd new_node_bounding_box;
    for (int dimension = 0; dimension < mins.size(); ++dimension) {
        new_node_bounding_box << mins.at(dimension), maxs.at(dimension);
    }
    node_bounding_box_.clear();
    node_bounding_box_ = new_node_bounding_box;
}

//-----------------------------------------------------------------------------

void Node::SetOverflowTreatment(bool overflow_treatment_setter)
{
    overflow_treatment_was_called_ = overflow_treatment_setter;
}

//-----------------------------------------------------------------------------

std::vector<Node*> Node::GetNodeChildren()
{
    return children_nodes_;
}

//-----------------------------------------------------------------------------

int Node::GetNumberOfChildren()
{
    return children_nodes_.size();
}

//-----------------------------------------------------------------------------

Eigen::MartixXd Node::GetNodeBounds()
{
    return node_bounding_box_;
}

//-----------------------------------------------------------------------------

std::vector<Meas> Node::GetMeasurements()
{
    return measurements_;
}

//-----------------------------------------------------------------------------

int Node::GetNumberOfMeasurements()
{
    return measurements_.size();
}

//-----------------------------------------------------------------------------

bool Node::WasOverflowTreatmentCalled()
{
    return overflow_treatment_was_called_;
}

//-----------------------------------------------------------------------------

bool Node::IsNodeALeafNode()
{
    if (measurements_ > 0) {
        return true;
    } else {
        return false;
    }
}

//-----------------------------------------------------------------------------

Eigen::MartrixXd Node::GetCenterOfBoundingBox()
{
    Eigen::MatrixXd center_of_bounding_box = Eigen::MatrixXd(node_bounding_box_.rows(),1);
    for (int dimension = 0; dimension < node_bounding_box_.rows(); ++dimension) {
        double center_value = (node_bounding_box_(dimension,0) + node_bounding_box_(dimension,1)) / 2;
        center_of_bounding_box << center_value;
    }
    return center_of_bounding_box;
}

};

} // namespace rransac