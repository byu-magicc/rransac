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

    for (int i = 0; i < 100; ++i) 
    {
        Meas m1;
        Meas m2;
        double dt(1.0);
        Transformation transformation;

        m1.data = Eigen::MatrixXd::Random(3,1);
        m2.data = Eigen::MatrixXd(3,1);
        m2.data = m1.data;
        transformation.T = Eigen::MatrixXd::Random(3,3);

        test_obj.TransformMeasurement(m1,transformation,dt);
        EXPECT_FLOAT_EQ(m1.data.norm(),(transformation.T * m2.data).norm()) << "Both transformations should be equal";
    }
}

//---------------------------------------------------------------------------------------

TEST(Measurement_CVTest, spatialdistance)
{
    Measurement_CV test_obj;
    DistanceType distType = kSpatial;
    Parameters P;

    for (int i = 0; i < 100; ++i) 
    {
        Meas m1;
        Meas m2;

        m1.data = Eigen::MatrixXd::Random(3,1);
        m2.data= Eigen::MatrixXd::Random(3,1);

        EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),(m1.data - m2.data).norm()) << "Spatial distance should be equal";
    }
}

//---------------------------------------------------------------------------------------

TEST(Measurement_CVTest, temporaldistance)
{
    Measurement_CV test_obj;
    DistanceType distType = kTemporal;
    Parameters P;

    for (int i = 0; i < 100; ++i) 
    {
        Meas m1;
        Meas m2;
        double X(500.0);

        m1.time_stamp = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / X));
        m2.time_stamp = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / X));

        EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),abs(m1.time_stamp - m2.time_stamp)) << "Temporal distance should be equal";
    }
}

//---------------------------------------------------------------------------------------

TEST(Measurement_CVTest, totaldistance)
{
    Measurement_CV test_obj;
    DistanceType distType = kTotal;
    Parameters P;

    for (int i = 0; i < 100; ++i) 
    {
        Meas m1;
        Meas m2;
        double X(500.0);
        
        m1.data = Eigen::MatrixXd::Random(3,1);
        m1.time_stamp = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / X));
        m2.data= Eigen::MatrixXd::Random(3,1);
        m2.time_stamp = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / X));

        EXPECT_FLOAT_EQ(test_obj.GetDistance(m1,m2,distType,P),
                        sqrt(pow((m1.data - m2.data).norm(),2) + pow(abs(m1.time_stamp - m2.time_stamp),2))) << "Total distance should be equal";
    }
}

} // namespace rransac