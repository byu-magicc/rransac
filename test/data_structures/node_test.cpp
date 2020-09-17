#include <gtest/gtest.h>
#include <math.h>

#include "common/measurement/measurement_base.h"
#include "common/parameters.h"
#include "data_structures/node.h"

namespace rransac 
{

//---------------------------------------------------------------------------------------

TEST(NodeTest, addandremovechildnodes) 
{
    std::shared_ptr<Node> n1(new Node());
    std::shared_ptr<Node> firstChildOfn1(new Node());
    std::shared_ptr<Node> secondChildOfn1(new Node());

    n1->AddChildNode(firstChildOfn1);
    EXPECT_EQ(n1->GetNumberOfChildren(),1) << "Node n1 should only have 1 child";
    n1->AddChildNode(secondChildOfn1);
    EXPECT_EQ(n1->GetNumberOfChildren(),2) << "Node n1 should now have 2 children";
    n1->RemoveAllChildNodes();
    EXPECT_EQ(n1->GetNumberOfChildren(),0) << "Node n1 should now have 0 children after removing them";
}

//---------------------------------------------------------------------------------------

TEST(NodeTest, addandremovechildmeasurements) 
{
    std::shared_ptr<Node> n1(new Node());
    Meas m1;
    Meas m2;

    n1->AddMeasurement(m1);
    EXPECT_EQ(n1->GetNumberOfMeasurements(),1) << "Node n1 should only have 1 child Meas object";
    n1->AddMeasurement(m2);
    EXPECT_EQ(n1->GetNumberOfMeasurements(),2) << "Node n1 should now only have 2 Meas object children";
    n1->RemoveAllMeasurements();
    EXPECT_EQ(n1->GetNumberOfMeasurements(),0) << "Node n1 should now have 0 Meas object children";
}

//---------------------------------------------------------------------------------------

TEST(NodeTest, updatenodebb) 
{
    std::shared_ptr<Node> n1(new Node());
    std::shared_ptr<Node> n2(new Node());
    std::shared_ptr<Node> root(new Node());

    // TESTING: Node Bounding Box size matches one added measurement data
    for (int i = 0; i < 10; ++i) {
        Meas m1;
        m1.data = Eigen::MatrixXd::Random(3,1);
        n1->AddMeasurement(m1);
        EXPECT_DOUBLE_EQ((n1->GetNodeBounds().col(0) - m1.data).norm(),0) << "Node bounds should equal data in Meas object";
        n1->RemoveAllMeasurements();
    }
    
    // TESTING: Node Bounding Box size matches multiple added measurement data
    for (int i = 0; i < 10; ++i) {
        Meas m1;
        Meas m2;
        Meas m3;
        int number_of_added_measurements = 3;
        m1.data = Eigen::MatrixXd::Random(3,1);
        m2.data = Eigen::MatrixXd::Random(3,1);
        m3.data = Eigen::MatrixXd::Random(3,1);
        n1->AddMeasurement(m1);
        n1->AddMeasurement(m2);
        n1->AddMeasurement(m3);

        // Find minimum and maximum values in each dimension (test is 3D) to correctly test Node Bounding Box
        Eigen::MatrixXd correct_bounding_box = Eigen::MatrixXd::Zero(3,2); // 3 dimensions, min & max columns
        double min;
        double max;
        for (int dim = 0; dim < 3; ++dim) {
            min = m1.data(dim,0);
            max = m1.data(dim,0);
            if (m2.data(dim,0) < min) {
                min = m2.data(dim,0);
            }
            if (m2.data(dim,0) > max) {
                max = m2.data(dim,0);
            }
            if (m3.data(dim,0) < min) {
                min = m3.data(dim,0);
            }               
            if (m3.data(dim,0) > max) {
                max = m3.data(dim,0);
            }
            correct_bounding_box(dim,0) = min;
            correct_bounding_box(dim,1) = max;
        }
        //TEST
        Eigen::MatrixXd nodeBB(n1->GetNodeBounds());
        Eigen::MatrixXd correctBB(correct_bounding_box);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                EXPECT_DOUBLE_EQ(nodeBB(i,j),correctBB(i,j)) << "Node bounds should have updated properly for multiple Meas objects";
            }
        }
        n1->RemoveAllMeasurements();
    }

    // TESTING: Node Bounding Box size matches one added Node with measurement child
    for (int i = 0; i < 10; ++i) {
        Meas m1;
        m1.data = Eigen::MatrixXd::Random(3,1);
        n1->AddMeasurement(m1);
        root->AddChildNode(n1);
        //TEST
        Eigen::MatrixXd nodeBB(root->GetNodeBounds());
        Eigen::MatrixXd correctBB(n1->GetNodeBounds());
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                EXPECT_DOUBLE_EQ(nodeBB(i,j),correctBB(i,j)) << "Node bounds should have same bounds as 1 child node";
            }
        }
        n1->RemoveAllMeasurements();
        root->RemoveAllChildNodes();
    }
    
    // TESTING: Node Bounding Box size matches multiple added Nodes with measurement children
    for (int i = 0; i < 10; ++i) {
        Meas m1;
        Meas m2;
        Meas m3;
        int number_of_added_measurements = 3;
        m1.data = Eigen::MatrixXd::Random(3,1);
        m2.data = Eigen::MatrixXd::Random(3,1);
        m3.data = Eigen::MatrixXd::Random(3,1);
        n1->AddMeasurement(m1);
        n1->AddMeasurement(m2);
        n2->AddMeasurement(m3);
        root->AddChildNode(n1);
        root->AddChildNode(n2);

        // Find minimum and maximum values in each dimension (test is 3D) to correctly test Node Bounding Box
        Eigen::MatrixXd correct_bounding_box = Eigen::MatrixXd::Zero(3,2); // 3 dimensions, min & max columns
        double min;
        double max;
        for (int dim = 0; dim < 3; ++dim) {
            min = m1.data(dim,0);
            max = m1.data(dim,0);
            if (m2.data(dim,0) < min) {
                min = m2.data(dim,0);
            }
            if (m2.data(dim,0) > max) {
                max = m2.data(dim,0);
            }
            if (m3.data(dim,0) < min) {
                min = m3.data(dim,0);
            }               
            if (m3.data(dim,0) > max) {
                max = m3.data(dim,0);
            }
            correct_bounding_box(dim,0) = min;
            correct_bounding_box(dim,1) = max;
        }
        //TEST
        Eigen::MatrixXd nodeBB(root->GetNodeBounds());
        Eigen::MatrixXd correctBB(correct_bounding_box);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                EXPECT_DOUBLE_EQ(nodeBB(i,j),correctBB(i,j)) << "Node bounds should have updated properly for multiple Meas & Node objects";
            }
        }
        n1->RemoveAllMeasurements();
        n2->RemoveAllMeasurements();
        root->RemoveAllChildNodes();
    }

}

//---------------------------------------------------------------------------------------

TEST(NodeTest, isnodealeafnode) 
{
    Meas m1;
    m1.data = Eigen::MatrixXd::Random(3,1);
    std::shared_ptr<Node> leafNode(new Node());
    std::shared_ptr<Node> parentNode(new Node());
    parentNode->AddChildNode(leafNode);

    EXPECT_TRUE(leafNode->IsNodeALeafNode()) << "Node has no children Nodes, thus it should be a leaf node";
    EXPECT_FALSE(parentNode->IsNodeALeafNode()) << "Node has a child Nodes, thus it is not a leaf node";
    leafNode->AddMeasurement(m1);
    EXPECT_TRUE(leafNode->IsNodeALeafNode()) << "Node now has one child measurements, but it should remain a leaf node";    
}

//---------------------------------------------------------------------------------------

TEST(NodeTest, getcenterofbb) 
{
    std::shared_ptr<Node> n1(new Node());

    // TESTING: Center of Node bounding box with only 1 measurement
    for (int i = 0; i < 10; ++i) {
        Meas m1;
        m1.data = Eigen::MatrixXd::Random(3,1);
        n1->AddMeasurement(m1);

        Eigen::MatrixXd nodeBB(n1->GetCenterOfBoundingBox());
        Eigen::MatrixXd correctBB(n1->GetNodeBounds());
        for (int i = 0; i < 3; ++i) {
            EXPECT_DOUBLE_EQ(nodeBB(i,0),correctBB(i,0)) << "With 1 measurement, center of Node's bounding box should be the same as the Node's bounding box";
        }
        n1->RemoveAllMeasurements();
    }

    // TESTING: Center of Node bounding box with multiple measurements
    for (int i = 0; i < 10; ++i) {
        Meas m1;
        Meas m2;
        Meas m3;
        int number_of_added_measurements = 3;
        m1.data = Eigen::MatrixXd::Random(3,1);
        m2.data = Eigen::MatrixXd::Random(3,1);
        m3.data = Eigen::MatrixXd::Random(3,1);
        n1->AddMeasurement(m1);
        n1->AddMeasurement(m2);
        n1->AddMeasurement(m3);

        // Find minimum and maximum values in each dimension (test is 3D) to correctly test Node Bounding Box
        Eigen::MatrixXd correct_center_bounding_box = Eigen::MatrixXd::Zero(3,1); // 3 dimensions
        double min;
        double max;
        for (int dim = 0; dim < 3; ++dim) {
            min = m1.data(dim,0);
            max = m1.data(dim,0);
            if (m2.data(dim,0) < min) {
                min = m2.data(dim,0);
            }
            if (m2.data(dim,0) > max) {
                max = m2.data(dim,0);
            }
            if (m3.data(dim,0) < min) {
                min = m3.data(dim,0);
            }               
            if (m3.data(dim,0) > max) {
                max = m3.data(dim,0);
            }
            correct_center_bounding_box(dim,0) = ( min + max ) / 2;
        }
        //TEST
        Eigen::MatrixXd nodeBB(n1->GetCenterOfBoundingBox());
        Eigen::MatrixXd correctBB(correct_center_bounding_box);
        for (int i = 0; i < 3; ++i) {
            EXPECT_DOUBLE_EQ(nodeBB(i,0),correctBB(i,0)) << "Center of node bounding box should equal average of max and min values of bounding box";
        }
        n1->RemoveAllMeasurements();
    }
}

} // namespace rransac