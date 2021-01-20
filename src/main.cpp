
// #include <typeinfo>
// #include <Eigen/Core>
// #include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
// #include <chrono>
// #include "system.h"
#include "state.h"
// #include "common/models/model_RN.h"
// #include "common/sources/source_RN.h"
// #include "common/transformations/transformation_null.h"
// #include "common/measurement/measurement_base.h"


// #include "common/models/model_base.h"
#include "common/models/model_RN.h"
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
    typedef typename tModel::template ModelTemplate<float, State::StateTemplate> Model;
    

    void Something() {
        
    }

};


using namespace rransac;
using namespace lie_groups;
int main(){

typedef State<Rn,double,3> HMM;
HMM::StateTemplate<float> state_f;
HMM::StateTemplate<double> state_d;
typedef ModelRN<HMM, TransformNULL> Model;
Model::ModelTemplate<float,HMM::StateTemplate> model_f;

typedef State<SE2,double,3> SMM;
typedef ModelSENPosVel<SMM,TransformNULL> ModelS;
ModelS::ModelTemplate<float,SMM::StateTemplate> models_f;

Blah<ModelS> blah;
blah.Something();

return 0;

} 

