#ifndef RRANSAC_COMMON_SOURCES_SOURCE_R2_POS_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_R2_POS_H_

#include "common/sources/source_base.h"

namespace rransac
{

/**
 * \class SourceR2Pos
 * This measurement source is used for a target whose manifold is R2 and only the position is observed.
 */ 

class SourceR2Pos : public SourceBase<lie_groups::Rn<2>,lie_groups::rn<2>> {

    typedef lie_groups::Rn<2> G;
    typedef lie_groups::rn<2> U;
    typedef lie_groups::State<G,U> State;

public:

/** Initializes the measurement source */
void Init(const SourceParameters& params, unsigned int source_id) override; 

/** Returns the jacobian of the observation function w.r.t. the state. */
Eigen::MatrixXd GetLinObsMatState() override {return H_;}  

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd GetLinObsMatSensorNoise() override {return V_;} 

/** Returns an estimated measurement according to the state. */
Eigen::MatrixXd GetEstMeas(const State& state) override {return state.g_.data_;}
    

private:

Eigen::Matrix<double,2,4> H_; /** < Constant jacobian observation matrix w.r.t. the states */
Eigen::Matrix<double,2,2> V_; /** < Constant jacobian observation matrix w.r.t. the noise */


};
    
} // namespace rransac


#endif // RRANSAC_COMMON_SOURCES_SOURCE_R2_POS_H_
