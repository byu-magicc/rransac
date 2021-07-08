
#include "lie_groups/state.h"
#include <Eigen/Dense>
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/common/measurement/measurement_base.h"
#include <iostream>
#include <tuple>
#include <variant>

using namespace lie_groups;

using namespace rransac;


// typedef R2_r2 State;
typedef SourceRN<R2_r2,MeasurementTypes::RN_POS,TransformNULL> Source;
typedef SourceContainer<Source> SC;
typedef ModelRN<SC> Model;



int main() {


Model track1;
typedef typename Model::template ModelTemplate<float> ModelT;
ModelT track2;

typedef typename R2_r2::template StateTemplate<float> StateT;

StateT state2;

typedef typename Source::template SourceTemplate<float> SourceT;

SourceT source2;

std::cout << "hi " << std::endl;




 return 0;   
}