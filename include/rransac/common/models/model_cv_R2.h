#ifndef RRANSAC_COMMON_MODELS_CV_R2_H_
#define RRANSAC_COMMON_MODELS_CV_R2_H_

#include "common/models/model_base.h"
#include "common/measurement/measurement_factory.h"

namespace rransac {

/**
 * \class ModelCvNipR2
 * This is a constant velocity model on R2
 * 
 * In order to work with any Lie group, the model base is a template class that requires 
 * The lie group G, lie algebra U and the measurement M
 * 
 * Note: You can use any model. 
 */ 

class ModelCvR2 : public ModelBase<lie_groups::Rn<2>, lie_groups::rn<2>, Meas> {

    typedef lie_groups::Rn<2> G;
    typedef lie_groups::rn<2> U;
    typedef Meas M;
    typedef lie_groups::State<G,U> State;
    typedef Eigen::Matrix<double, G::size1_ + U::size1_, G::size1_ + U::size1_> Mat;

public:

    void Init(const Parameters& params) override;    

    State PropagateState(const State& state, const double dt) override;

    Mat GetLinTransFuncMatState(const State& state, const double dt) override;

    Mat GetLinTransFuncMatNoise(const State& state, const double dt) override;

    void PropagateModel(const double dt) override;

    void UpdateModel(const Parameters& param) override;

    Eigen::MatrixXd GetLinObsMatState(const unsigned int source_ID) override;

    Eigen::MatrixXd GetLinObsMatSensorNoise(const double dt, const unsigned int source_ID) override;

    Eigen::MatrixXd GetEstMeas(const State& state, const unsigned int source_ID) override;

    void TransformModel(const Transformation& T, const double dt) override;

private: 

    Mat F_; /** < The jacobian of the system function w.r.t. the state. */
    Mat G_; /** < The jacobian of the system function w.r.t. the noise. */
    Eigen::Matrix2d identity_2d_; /** < 2x2 identity matrix */


};


}


#endif // RRANSAC_COMMON_MODELS_CV_R2_H_