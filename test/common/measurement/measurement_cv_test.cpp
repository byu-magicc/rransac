#include <gtest/gtest.h>
#include <math.h>

#include "common/measurement/measurement_cv.h"
#include "common/parameters.h"
#include "common/transformation.h"

namespace rransac 
{

//---------------------------------------------------------------------------------------

TEST(Measurement_CVTest, transformmeasurement) 
{
    Measurement_CV test_obj;
    DistanceType distType = kSpatial;
    Parameters P;

    Meas m1, m2;
    Transformation transformation;
    double dt(1.0);
    m1.data = Eigen::MatrixXd(3,1);
    m1.data << 1,1,1;
    m2.data = Eigen::MatrixXd(3,1);
    m2.data << 1,1,1;
    transformation.T = Eigen::MatrixXd(3,3);
    transformation.T << 1,0,0,
                        0,1,0,
                        0,0,1;

    // Transformation is applied to m1 and should do nothing
    test_obj.TransformMeasurement(m1,transformation,dt);
    EXPECT_FLOAT_EQ((float)(m1.data-m2.data).norm(),0.0) << "Should be 0.";

    // 90 degree X-direction rotation Transformation is applied to m1
    transformation.T << 1,0,0,
                        0,cos(M_PI),-sin(M_PI),
                        0,sin(M_PI),cos(M_PI);
    test_obj.TransformMeasurement(m1,transformation,dt);
    EXPECT_FLOAT_EQ((float)(m1.data-m2.data).norm(),2*sqrt(2)) << "Should be 2*sqrt(2).";
}

//---------------------------------------------------------------------------------------

TEST(Measurement_CVTest, spatialdistance)
{
    Measurement_CV test_obj;
    DistanceType distType = kSpatial;
    Parameters P;

    Meas m1;
    Meas m2;
    m1.data = Eigen::MatrixXd(3,1);
    m1.data << 1,1,1;
    m2.data= Eigen::MatrixXd(3,1);
    m2.data << 1,1,1;

    // Testing same point location == 0.0 distance
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),0.0) << "Spatial distance should equal 0.";
    
    // Testing different points including negative data entry values
    m2.data(0,0) = 2;
    m2.data(1,0) = -3;
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),(m1.data-m2.data).norm()) << "Spatial distances should be equal.";

    // nD Matrix tests:
    Meas m3;
    Meas m4;
    m3.data = Eigen::MatrixXd(3,3);
    m3.data << 1,4,1,
                -2,3,1,
                3,1,5;
    m4.data= Eigen::MatrixXd(3,3);
    m4.data << 3,-1,10,
                1,1,1,
                2,1,-12;
    
    // Testing different points nD including negative data entry values
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m3,m4,distType,P),(m3.data-m4.data).norm()) << "Spatial distances should be equal.";
}

//---------------------------------------------------------------------------------------

TEST(Measurement_CVTest, temporaldistance)
{
    Measurement_CV test_obj;
    DistanceType distType = kTemporal;
    Parameters P;

    Meas m1;
    Meas m2;
    m1.time_stamp = 0.0;
    m2.time_stamp = 0.0;

    // Testing same temporal values == 0.0 distance
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),0.0) << "Temporal Distance should equal 0.";
    
    // Testing different integer time steps
    m2.time_stamp = 4.0;
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),4) << "Temporal distance should equal 4.";

    // Testing different decimal time steps
    m1.time_stamp = 0.5;
    m2.time_stamp = 4.25;
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),3.75) << "Temporal distance should equal 3.75.";
}

//---------------------------------------------------------------------------------------

TEST(Measurement_CVTest, totaldistance)
{
    Measurement_CV test_obj;
    DistanceType distType = kTotal;
    Parameters P;

    Meas m1;
    Meas m2;
    m1.data = Eigen::MatrixXd(3,1);
    m1.data << 1,1,1;
    m1.time_stamp = 0.0;
    m2.data = Eigen::MatrixXd(3,1);
    m2.data << 1,1,1;
    m2.time_stamp = 0.0;
    
    // Testing same temporal and spatial values == 0.0 distance
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),0.0) << "Total Distance should equal 0.";
    
    // Testing different integer temporal and spatial values
    m2.time_stamp = 4.0;
    m2.data(1,0) = 3.0;
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),sqrt(pow(4,2)+pow(2,2))) << "Total distance should equal sqrt(20).";

    // Testing different decimal temporal and spatial values
    m1.time_stamp = 0.5;
    m1.data(0,0) = 1.5;
    m2.time_stamp = 4.0;
    m2.data(1,0) = 3.5;
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),sqrt(pow(3.5,2)+pow(sqrt(pow(2.5,2)+pow(.5,2)),2))) << "Total distance should be equal.";
    
    // nD Matrix tests:
    Meas m3;
    Meas m4;
    m3.time_stamp = 1.5;
    m4.time_stamp = 1;
    m3.data = Eigen::MatrixXd(3,3);
    m3.data << 1,4,1,
                -2,3,1,
                3,1,5;
    m4.data= Eigen::MatrixXd(3,3);
    m4.data << 3,-1,10,
                1,1,1,
                2,1,-12;
    
    // Testing different points nD including negative data entry values
    EXPECT_FLOAT_EQ(test_obj.GetDistance(m3,m4,distType,P),sqrt(pow((m3.data-m4.data).norm(),2)+pow(.5,2))) << "Total distances should be equal.";

}

} // namespace rransac