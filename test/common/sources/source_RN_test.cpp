#include <gtest/gtest.h>
#include "common/sources/source_RN.h"

namespace rransac
{
    

TEST(Source_RN, INIT){

// Make sure that not every type can be initialized
SourceParameters params;
SourceRN<lie_groups::SE2_se2> source_se2;
ASSERT_ANY_THROW(source_se2.Init(params));






}






} // namespace rransac
