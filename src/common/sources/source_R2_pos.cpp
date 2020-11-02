#include <common/sources/source_R2_pos.h>

namespace rransac
{

void SourceR2Pos::Init(const SourceParameters& params, unsigned int source_id) {
    params_ = params;
    source_id_ = source_id;
    H_.block(0,0,2,2).setIdentity();
    H_.block(0,2,2,2).setZero();
    V_.setIdentity();
}


} // namespace rransac
