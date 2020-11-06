#ifndef RRANSAC_RRANSAC_H_
#define RRANSAC_RRANSAC_H_

#include "system.h"

namespace rransac {


/**
 * \class RRANSAC
 * This class is the interface for the user. 
 * 
 */ 

class RRANSAC {

public:

    System sys_;

    /**
     * This method adds a source to the vector Parameters::sources_. Before adding a source,
     * it verifies that the source has a unique ID and that the member variables of the source
     * are set properly.
     * @return Returns true if the source was added.
     */
    bool AddSource(const SourceParameters& params);

    /**
     * \detail Sets all of the parameters
     * \param[in] new_params The new parameters.
     * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
     */
    bool SetParameters(const Parameters &new_params) { return sys_.params_.SetParameters(new_params);  }



private:


};

}

#endif // RRANSAC_RRANSAC_H_