#ifndef RRANSAC_SYSTEM_H_
#define RRANSAC_SYSTEM_H_

#include "parameters.h"
#include "common/sources/source_base.h"

namespace rransac {


/**
 * \class System
 * This class holds all of the system information including parameters, sources, and models. 
 * 
 */ 
template <typename S> 
class System {



public:

    Parameters params_;               /** < System parameters */
    std::vector<SourceBase> sources_; /** < Contains all of the instantiated sources. */
    std::vector<Meas> new_meas_;

private:

};

template <typename Source>
class node {

void computeBoundingBox<const std::vector<Source>& sources>

    Source::Mat_A = sources[0].ToAlgebra(meas)

}

}

#endif // RRANSAC_SYSTEM_H_
