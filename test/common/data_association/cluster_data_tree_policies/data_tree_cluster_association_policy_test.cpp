#include <gtest/gtest.h>
#include "system.h"
#include "common/models/model_RN.h"
#include "common/sources/source_RN.h"
#include "common/transformations/transformation_null.h"
#include "common/data_association/data_association_host.h"
#include "common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h"



using namespace lie_groups;
using namespace rransac;

// This is a dummy policy that does nothing. It is just need.
template<typename tModel>
class ModelDummyPolicy {

public:

static void  PolicyDataAssociationModel(System<tModel>& sys){}

};

TEST(DataTreeClusterAssociationPolicyTest, PolicyDataAssociationClusterDataTree) {  

    typedef ModelRN<R2_r2, TransformNULL> Model;
    typedef Meas<double> Measurement;

    SourceParameters source_params;
    source_params.expected_num_false_meas_ = 0.1;
    source_params.meas_cov_ = Eigen::Matrix2d::Identity();
    source_params.gate_probability_ = 0.8;
    source_params.source_index_ = 0;
    source_params.probability_of_detection_ = 0.8;
    source_params.type_ = MeasurementTypes::RN_POS;

    typename Model::Source source;
    source.Init(source_params);

    System<Model> sys;
    sys.params_.cluster_position_threshold_ = 0.1;
    sys.params_.cluster_time_threshold_ = 1;
    sys.params_.cluster_velocity_threshold_ = 0.2;
    sys.sources_.push_back(source);
    DataAssociationHost<Model, ModelDummyPolicy, DataTreeClusterAssociationPolicy> data_association;
    
    Measurement m;
    m.source_index = 0;
    unsigned int num_meas= 1000;
    std::list<Measurement> measurements(num_meas);

    for(auto iter = measurements.begin(); iter != measurements.end(); ++iter) {
        iter->time_stamp = 0;
        iter->pose = Eigen::Matrix<double,2,1>::Random();
        iter->source_index = 0;
        
    }



    sys.new_meas_ = measurements;
    data_association.AssociateNewMeasurements(sys);

    ASSERT_EQ(sys.data_tree_.Size(), num_meas);
    ASSERT_EQ(sys.new_meas_.size(), 0);




    // Generate a bunch of measurements;


}