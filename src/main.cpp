
#include <rransac/rransac.h>
#include <lie_groups/state.h>
#include <rransac/common/transformations/trans_homography.h>
#include <rransac/common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h>
#include <rransac/common/data_association/model_policies/model_pdf_policy.h>
#include <rransac/common/measurement/measurement_base.h>
#include <rransac/common/sources/source_base.h>

#if TRACKING_SE2
  #include <rransac/common/sources/source_SEN_pos_vel.h>
  #include <rransac/common/models/model_SEN_pos_vel.h>
  #include <rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h>
  #include <rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h> 
#else // R2
  #include <rransac/common/sources/source_RN.h>
  #include <rransac/common/models/model_RN.h>
  #include <rransac/track_initialization/seed_policies/null_seed_policy.h>
  #include <rransac/track_initialization/lmle_policies/linear_lmle_policy.h>
#endif


#if TRACKING_SE2
  typedef rransac::RRANSACTemplateParameters<lie_groups::SE2_se2, rransac::SourceSENPosVel,rransac::TransformHomography, rransac::ModelSENPosVel,rransac::SE2PosSeedPolicy,rransac::NonLinearLMLEPolicy,rransac::ModelPDFPolicy,rransac::DataTreeClusterAssociationPolicy> RRTemplateParams;

#else
  typedef rransac::RRANSACTemplateParameters<lie_groups::R2_r2, rransac::SourceRN,rransac::TransformHomography, rransac::ModelRN,rransac::NULLSeedPolicy,rransac::LinearLMLEPolicy,rransac::ModelPDFPolicy,rransac::DataTreeClusterAssociationPolicy> RRTemplateParams;
#endif

typedef typename RRTemplateParams::Model RR_Model;
typedef typename RRTemplateParams::State RR_State;
typedef typename rransac::System<RR_Model> RR_System;
typedef typename rransac::RRANSAC<RRTemplateParams> RR_RRANSAC;

#include "blah.h"

int main() {

    RR_RRANSAC rransac_; /**< RRANSAC is a multiple target tracking algorithm. */
    rransac::Parameters rransac_params_;              /**< The parameters needed for RRANSAC. */
    const RR_System* rransac_sys_; /**< A constant pointer to all of the data of RRANSAC. */


    rransac_params_.process_noise_covariance_ = Eigen::Matrix<double,4,4>::Identity();
    rransac_.SetSystemParameters(rransac_params_);

    Blah();

 return 0;   
}