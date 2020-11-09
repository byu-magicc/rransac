#include "common/sources/source_RN.h"

namespace rransac
{
    

template<class S>
void SourceRN<S>::Init(const SourceParameters& params) {

    const unsigned int sizeg = S::g_type_::size1;
    const unsigned int sizeu = S::u_type_::size1;

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

    params_ = params;

}


} // namespace rransac
