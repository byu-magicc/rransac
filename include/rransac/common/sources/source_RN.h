#include "common/sources/source_base.h"
#include <typeinfo>
#include "common/sources/source_base.h"

namespace rransac {

template<class S>
class SourceRN : public SourceBase<S,SourceRN<S>> {

public:

// typedef S type_;

/** Initializes the measurement source. This function must set the parameters.  */
void Init(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states */
Eigen::MatrixXd GetLinObsMatState(S const& state) {return this->H_;}                        

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state) {return this->V_;}                         

/** Computes the estimated measurement given a state */
Eigen::MatrixXd GetEstMeas(const S& state) {return state.g_.data_;} /** Returns an estimated measurement according to the state. */



};


template<class S>
void SourceRN<S>::Init(const SourceParameters& params) {

    const unsigned int sizeg = S::g_type_::size1_;
    const unsigned int sizeu = S::u_type_::size1_;

    // Make sure that the state is a valid type
    if (typeid(S).name() != typeid(lie_groups::State<lie_groups::Rn<sizeg>,lie_groups::Rn<sizeu>>).name())
    {
        std::runtime_error("SourceRNPos::Init State is not supported by this source type");
    }

    // Construct the Jacobians
    switch (params.type_)
    {
    case MeasurementTypes::RN_POS:
        this->H_ = Eigen::Matrix<double,sizeg,sizeg+sizeu>::Zero();
        this->H_.block(0,0,sizeg,sizeg).setIdentity();
        this->V_ = Eigen::Matrix<double,sizeg,sizeg>::Identity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        this->H_ = Eigen::Matrix<double,sizeg+sizeu,sizeg+sizeu>::Zero();
        this->H_.block(0,0,sizeg,sizeg).setIdentity();
        this->V_ = Eigen::Matrix<double,sizeg+sizeu,sizeg+sizeu>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN:Init Measurement type not supported.");
        break;
    }

    this->params_ = params;

}

// template<class S>
// Eigen::MatrixXd SourceRN<S>::GetLinObsMatState(S const& state) {return this->H_;}


} // namesapce rransac