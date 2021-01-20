
// #include <typeinfo>
// #include <Eigen/Core>
// #include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <cmath>
// #include <chrono>
// #include "system.h"
#include "state.h"
// #include "common/models/model_RN.h"
// #include "common/sources/source_RN.h"
// #include "common/transformations/transformation_null.h"
// #include "common/measurement/measurement_base.h"


// #include "common/models/model_base.h"
// #include "common/models/model_RN.h"
#include "common/models/model_SEN_pos_vel.h"
// #include "common/models/model_SEN_pose_twist.h"
// #include "common/sources/source_base.h"
// #include "common/sources/source_RN.h"
// #include "common/sources/source_SEN_pos_vel.h"
// #include "common/sources/source_SEN_pose_twist.h"
// #include "parameters.h"
// #include "common/transformations/transformation_base.h"
// #include "common/transformations/trans_homography.h"
#include "common/transformations/transformation_null.h"



// template<typename T>
// struct Meh<Blah<T>>{
//     void print() {Blah<T>::print(); Blah<T>::print();}
// };

template<typename tModel>
struct Blah {

    typedef typename tModel::State State;
    typedef typename tModel::template ModelTemplate<float, State::template StateTemplate>::State StateF;
    

    void Something() {
        StateF state;
        double dt = 0.1;
        tModel::template ModelTemplate<float, State::template StateTemplate>::PropagateState(state,dt);
    }

};


// template <typename T, typename S>
// struct Blah {};

// template <typename T, template<typename> typename S>
// struct Hah{};

using namespace rransac;
using namespace lie_groups;
int main(){

// Hah<double,Blah<S>>;

// typedef State<Rn,double,3> HMM;
// HMM::StateTemplate<float> state_f;
// HMM::StateTemplate<double> state_d;
// typedef ModelRN<HMM, TransformNULL> Model;
// Model::ModelTemplate<float,HMM::StateTemplate> model_f;

// typedef State<SE2,double,3> SMM;
// typedef ModelSENPosVel<SMM,TransformNULL> ModelS;
// ModelS::ModelTemplate<float,SMM::StateTemplate> models_f;

// Blah<ModelS> blah;
// blah.Something();

std::cout << std::abs(-1.178965) << std::endl;

return 0;

} 

