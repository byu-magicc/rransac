#include <gtest/gtest.h>
#include "system.h"
#include "common/models/model_RN.h"
#include "common/sources/source_RN.h"
#include "common/transformations/transformation_null.h"
#include "common/data_association/data_association_host.h"
#include "common/data_association/model_policies/model_data_assoc_pdf_policy.h

template<typename tModel>
class DummyClusterDataTreeAssocPolicy {
public:

    static void  PolicyDataAssociationClusterDataTree(System<tModel>& sys) {
        // Store the size as the current time so we can verify it.
        sys.current_time_ = sys.new_meas_.size();
        sys.new_meas_.clear();
    }

};


TEST(ModelPDFPolicyTest, PolicyDataAssociationModel) {









}