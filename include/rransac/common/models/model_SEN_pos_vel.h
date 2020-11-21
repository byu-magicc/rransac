// state = Model::Model::MB_State::Random();

// // Construct a system that is observable
// if (Model::MeasType1 == MeasurementTypes::SEN_POS || Model::MeasType1 == MeasurementTypes::SEN_POS_VEL) {

//     state.g_.R_.block(0,0,state.u_.p_.rows(),1) = state.u_.p_.normalized(); 

//     if(state.g_.dim_ == 3) {
//         state.g_.R_.block(0,1,state.u_.p_.rows(),1) << - state.g_.data_(1,0), state.g_.data_(0,0);
//     } else {

//         state.g_.R_.block(0,1,state.u_.p_.rows(),1) << -state.g_.data_(1,0), state.g_.data_(0,0), state.g_.data_(2,0);
//         state.g_.R_.block(0,2,state.u_.p_.rows(),1) = Model::Model::MB_State::g_type_::rot_algebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1))*state.g_.data_.block(0,1,state.u_.p_.rows(),1);

//     }

//     double px = state.u_.p_.norm();
//     state.u_.p_.setZero();
//     state.u_.p_(0,0) = px;
// }