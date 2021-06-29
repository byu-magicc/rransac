
#include "lie_groups/state.h"
#include <Eigen/Dense>
// #include "rransac/common/models/model_SEN_pos_vel.h"
// #include "rransac/common/sources/source_SEN_pos_vel.h"
// #include "rransac/common/transformations/transformation_null.h"
#include <iostream>
#include <tuple>
#include <variant>

using namespace lie_groups;


template<typename State>
struct Test {

typedef Test<State> TMP;
typedef typename State::template StateTemplate<int> S;


template<typename T>
using TestTemplate = Test<typename State::template StateTemplate<T>>;

};





int main() {


lie_groups::R2_r2 state;

lie_groups::R2_r2::StateTemplate<float> state1;
typedef lie_groups::R2_r2 State;

Test<State::StateTemplate<float>> test;





 return 0;   
}