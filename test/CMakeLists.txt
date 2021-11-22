cmake_minimum_required(VERSION 3.14)

########################################
# Install googletest in local directory
########################################

# set(GTEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)

# # Download and install GoogleTest
# configure_file(
# 	"${PROJECT_SOURCE_DIR}/cmake/googletest.cmake.in"
# 	"${GTEST_DIR}/CMakeLists.txt"
# 	IMMEDIATE @ONLY)

# execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${GTEST_DIR} )
# execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${GTEST_DIR})

# # Prevent overriding the parent project's compiler/linker settings on Windows
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# # Add googletest directly to the build. This defines the `gtest` and `gtest_main` targets. Don't add these
# # targets to the ALL target so that `make install` will ignore googletest installation targets
# add_subdirectory(${GTEST_DIR}/googletest-src ${GTEST_DIR}/googletest-build EXCLUDE_FROM_ALL)

include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


include(CTest)
########################################
# Define Tests
########################################

include_directories(../include
					../liegroups/include
					../test
					${EIGEN3_INCLUDE_DIR}
					${Boost_INCLUDE_DIRS}
					${OpenCV_INCLUDE_DIRS})

# ---------------------------------------
# Parameter Test
# ---------------------------------------

add_executable(parameterTest 
parameters_test.cpp
../src/parameters.cpp)
target_link_libraries(parameterTest gtest_main)
add_test(NAME AllTestsInParameterTest COMMAND parameterTest)

# ---------------------------------------
# Cluster 
# ---------------------------------------

add_executable(cluster_Test 
data_containers/cluster_test.cpp
../src/parameters.cpp)
target_link_libraries(cluster_Test gtest_main )
add_test(NAME AllTestsInClusterTest COMMAND cluster_Test)


# ---------------------------------------
# Source tests
# ---------------------------------------
add_executable(source_baseTest
common/sources/source_base_test.cpp
../src/parameters.cpp)
target_link_libraries(source_baseTest gtest_main )
add_test(NAME AllTestsInsSourceBaseTest COMMAND source_baseTest)


add_executable(source_RNTest
common/sources/source_RN_test.cpp
../src/parameters.cpp)
target_link_libraries(source_RNTest gtest_main )
add_test(NAME AllTestsInSourceRNTest COMMAND source_RNTest)

add_executable(source_R2_R3_radarTest
common/sources/source_R2_R3_radar_test.cpp
../src/parameters.cpp)
target_link_libraries(source_R2_R3_radarTest gtest_main )
add_test(NAME AllTestsInSourceRNRadarTest COMMAND source_R2_R3_radarTest)


add_executable(source_SEN_Pos_VelTest
common/sources/source_SEN_pos_vel_test.cpp
../src/parameters.cpp)
target_link_libraries(source_SEN_Pos_VelTest gtest_main )
add_test(NAME AllTestsInSourceSENPosVelTest COMMAND source_SEN_Pos_VelTest)


add_executable(source_SEN_Pose_TwistTest
common/sources/source_SEN_pose_twist_test.cpp
../src/parameters.cpp)
target_link_libraries(source_SEN_Pose_TwistTest gtest_main )
add_test(NAME AllTestsInSourceSENPoseTwistTest COMMAND source_SEN_Pose_TwistTest)

add_executable(source_SE3_CamDepthTest
common/sources/source_SE3_cam_depth_test.cpp
../src/parameters.cpp)
target_link_libraries(source_SE3_CamDepthTest gtest_main )
add_test(NAME AllTestsInSourceSE3CamDepthTest COMMAND source_SE3_CamDepthTest)

add_executable(source_SE2_SE3_RadarTest
common/sources/source_SE2_SE3_radar_test.cpp
../src/parameters.cpp)
target_link_libraries(source_SE2_SE3_RadarTest gtest_main )
add_test(NAME AllTestsInSourceSE2SE3RAdarTest COMMAND source_SE2_SE3_RadarTest)

add_executable(source_containerTest
common/sources/source_container_test.cpp
../src/parameters.cpp)
target_link_libraries(source_containerTest gtest_main )
add_test(NAME AllTestsInSourceContainerTest COMMAND source_containerTest)

# ---------------------------------------
# Consensus Set Test
# ---------------------------------------
add_executable(consensus_setTest 
data_containers/consensus_set_test.cpp)
target_link_libraries(consensus_setTest gtest_main)
add_test(NAME AllTestsInConsensusSetTest COMMAND consensus_setTest)

# ---------------------------------------
# Transformations
# ---------------------------------------

add_executable(trans_homographyTest
common/transformations/transformation_homography_test.cpp
../src/parameters.cpp)
target_link_libraries(trans_homographyTest gtest_main )

add_executable(trans_radar_R2_R3_with_SE2_SE3Test
common/transformations/trans_radar_R2_R3_with_SE2_SE3_test.cpp
../src/parameters.cpp)
target_link_libraries(trans_radar_R2_R3_with_SE2_SE3Test gtest_main )

add_executable(trans_radar_SE2_SE3_with_SE2_SE3Test
common/transformations/trans_radar_SE2_SE3_with_SE2_SE3_test.cpp
../src/parameters.cpp)
target_link_libraries(trans_radar_SE2_SE3_with_SE2_SE3Test gtest_main )

add_executable(trans_SE3_cam_depthTest
common/transformations/trans_SE3_cam_depth_test.cpp
../src/parameters.cpp)
target_link_libraries(trans_SE3_cam_depthTest gtest_main )

# ---------------------------------------
# Model
# ---------------------------------------

add_executable(centralized_measurement_fusionTest 
common/models/centralized_measurement_fusion_test.cpp
../src/parameters.cpp)
target_link_libraries(centralized_measurement_fusionTest gtest_main )

add_executable(model_baseTest
common/models/model_base_test.cpp
../src/parameters.cpp)
target_link_libraries(model_baseTest gtest_main )

add_executable(model_baseTest2
common/models/model_base_test2.cpp
../src/parameters.cpp)
target_link_libraries(model_baseTest2 gtest_main )

add_executable(model_rn_tangentTest
common/models/model_rn_tangent_test.cpp
../src/parameters.cpp)
target_link_libraries(model_rn_tangentTest gtest_main )

add_executable(model_managerTest 
common/models/model_manager_test.cpp
../src/parameters.cpp)
target_link_libraries(model_managerTest gtest_main )




# ---------------------------------------
# Data Trees
# ---------------------------------------

add_executable(data_tree_clustersTest 
data_containers/data_tree/data_tree_clusters_test.cpp
../src/parameters.cpp)
target_link_libraries(data_tree_clustersTest gtest_main )


# ---------------------------------------
# Data Association
# ---------------------------------------

add_executable(data_association_host_test 
common/data_association/data_association_host_test.cpp
../src/parameters.cpp)
target_link_libraries(data_association_host_test gtest_main )

add_executable(validation_region_policy_tests 
common/data_association/validation_region_policies/validation_region_policy_tests.cpp
../src/parameters.cpp)
target_link_libraries(validation_region_policy_tests gtest_main )

add_executable(tli_ipdaf_policy_test 
common/data_association/track_likelihood_info_policies/tli_ipdaf_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(tli_ipdaf_policy_test gtest_main )

add_executable(mw_ipdaf_policy_test 
common/data_association/measurement_weight_policy/mw_ipdaf_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(mw_ipdaf_policy_test gtest_main )

# ---------------------------------------
# Track Initialization
# ---------------------------------------

add_executable(linear_lmle_policyTest 
track_initialization/lmle_policies/linear_lmle_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(linear_lmle_policyTest gtest_main )

add_executable(nonlinear_lmle_policyTest 
track_initialization/lmle_policies/nonlinear_lmle_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(nonlinear_lmle_policyTest gtest_main Ceres::ceres)

add_executable(SE2_pos_policyTest 
track_initialization/seed_policies/SE2_pos_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(SE2_pos_policyTest gtest_main
)

add_executable(SE3_cam_depth_seed_policyTest 
track_initialization/seed_policies/SE3_cam_depth_seed_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(SE3_cam_depth_seed_policyTest gtest_main
)

add_executable(radar_R2_R3_seed_policyTest 
track_initialization/seed_policies/radar_R2_R3_seed_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(radar_R2_R3_seed_policyTest gtest_main
)

add_executable(radar_SE2_SE3_seed_policyTest 
track_initialization/seed_policies/radar_SE2_SE3_seed_policy_test.cpp
../src/parameters.cpp)
target_link_libraries(radar_SE2_SE3_seed_policyTest gtest_main
)

add_executable(ransacTest 
track_initialization/ransac_test.cpp
../src/parameters.cpp)
target_link_libraries(ransacTest gtest_main )


add_executable(ransacTest2 
track_initialization/ransac_test2.cpp
../src/parameters.cpp)
target_link_libraries(ransacTest2 gtest_main Ceres::ceres)


add_executable(ransacTest3 
track_initialization/ransac_test3.cpp
../src/parameters.cpp)
target_link_libraries(ransacTest3 gtest_main Ceres::ceres)

# ---------------------------------------
# RRANSAC
# ---------------------------------------

add_executable(rransacTest
rransac_test.cpp
../src/parameters.cpp)
target_link_libraries(rransacTest gtest_main)

add_executable(rransacFullTest
rransac_full_test.cpp
../src/parameters.cpp)
target_link_libraries(rransacFullTest gtest_main Ceres::ceres)











