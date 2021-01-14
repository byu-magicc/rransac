cmake_minimum_required(VERSION 3.8)

########################################
# Install googletest in local directory
########################################

set(GTEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)

# Download and install GoogleTest
configure_file(
	"${PROJECT_SOURCE_DIR}/cmake/googletest.cmake.in"
	"${GTEST_DIR}/CMakeLists.txt"
	IMMEDIATE @ONLY)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${GTEST_DIR} )
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${GTEST_DIR})

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to the build. This defines the `gtest` and `gtest_main` targets. Don't add these
# targets to the ALL target so that `make install` will ignore googletest installation targets
add_subdirectory(${GTEST_DIR}/googletest-src ${GTEST_DIR}/googletest-build EXCLUDE_FROM_ALL)
  
########################################
# Define Tests
########################################

#---------------------------------------
# Parameter Test
#---------------------------------------

# add_executable(parameterTest 
# parameters_test.cpp
# ../src/parameters.cpp)
# target_link_libraries(parameterTest gtest_main)

#---------------------------------------
# Consensus Set Test
#---------------------------------------
# add_executable(consensus_setTest 
# data_containers/consensus_set_test.cpp)
# target_link_libraries(consensus_setTest gtest_main)

#---------------------------------------
# Source tests
#---------------------------------------
# add_executable(source_baseTest
# common/sources/source_base_test.cpp
# ../src/parameters.cpp)
# add_dependencies(source_baseTest lie_groups)
# target_link_libraries(source_baseTest gtest_main lie_groups)

# add_executable(source_RNTest
# common/sources/source_RN_test.cpp
# ../src/parameters.cpp)
# add_dependencies(source_RNTest lie_groups)
# target_link_libraries(source_RNTest gtest_main lie_groups)

# add_executable(source_SEN_Pos_VelTest
# common/sources/source_SEN_pos_vel_test.cpp
# ../src/parameters.cpp)
# add_dependencies(source_SEN_Pos_VelTest lie_groups)
# target_link_libraries(source_SEN_Pos_VelTest gtest_main lie_groups)

# add_executable(source_SEN_Pose_TwistTest
# common/sources/source_SEN_pose_twist_test.cpp
# ../src/parameters.cpp)
# add_dependencies(source_SEN_Pose_TwistTest lie_groups)
# target_link_libraries(source_SEN_Pose_TwistTest gtest_main lie_groups)

#---------------------------------------
# Transformations
#---------------------------------------

# add_executable(trans_homographyTest
# common/transformations/transformation_homography_test.cpp
# ../src/parameters.cpp)
# add_dependencies(trans_homographyTest lie_groups)
# target_link_libraries(trans_homographyTest gtest_main lie_groups)

#---------------------------------------
# Model
#---------------------------------------

# add_executable(model_baseTest
# common/models/model_base_test.cpp
# ../src/parameters.cpp)
# add_dependencies(model_baseTest lie_groups)
# target_link_libraries(model_baseTest gtest_main lie_groups)

# add_executable(model_baseTest2
# common/models/model_base_test2.cpp
# ../src/parameters.cpp)
# add_dependencies(model_baseTest2 lie_groups)
# target_link_libraries(model_baseTest2 gtest_main lie_groups)

# add_executable(model_managerTest 
# common/models/model_manager_test.cpp
# ../src/parameters.cpp)
# add_dependencies(model_managerTest lie_groups)
# target_link_libraries(model_managerTest gtest_main lie_groups)

#---------------------------------------
# Cluster 
#---------------------------------------

# add_executable(cluster_listTest 
# data_containers/cluster_test.cpp
# ../src/parameters.cpp)
# add_dependencies(cluster_listTest lie_groups)
# target_link_libraries(cluster_listTest gtest_main lie_groups)

#---------------------------------------
# Data Trees
#---------------------------------------

# add_executable(data_tree_clustersTest 
# data_containers/data_tree/data_tree_clusters_test.cpp
# ../src/parameters.cpp)
# add_dependencies(data_tree_clustersTest lie_groups)
# target_link_libraries(data_tree_clustersTest gtest_main lie_groups)


#---------------------------------------
# Data Association
#---------------------------------------
# add_executable(cluster_data_tree_data_associationTest 
# common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy_test.cpp
# ../src/parameters.cpp)
# add_dependencies(cluster_data_tree_data_associationTest lie_groups)
# target_link_libraries(cluster_data_tree_data_associationTest gtest_main lie_groups)

# add_executable(model_pdf_policyTest 
# common/data_association/model_policies/model_pdf_policy_test.cpp
# ../src/parameters.cpp)
# add_dependencies(model_pdf_policyTest lie_groups)
# target_link_libraries(model_pdf_policyTest gtest_main lie_groups)

#---------------------------------------
# Track Initialization
#---------------------------------------

# add_executable(linear_lmle_policyTest 
# track_initialization/lmle_policies/linear_lmle_policy_test.cpp
# ../src/parameters.cpp)
# add_dependencies(linear_lmle_policyTest lie_groups)
# target_link_libraries(linear_lmle_policyTest gtest_main lie_groups)

# add_executable(ransacTest 
# track_initialization/ransac_test.cpp
# ../src/parameters.cpp)
# add_dependencies(ransacTest lie_groups)
# target_link_libraries(ransacTest gtest_main lie_groups)







add_executable(rransacTest
rransac_test.cpp
../src/parameters.cpp
../src/rransac.cpp)
target_link_libraries(rransacTest gtest_main lie_groups)













# add_executable(nodeTest
# data_structures/node_test.cpp
# ../src/data_structures/node.cpp)
# target_link_libraries(nodeTest gtest_main)

# add_custom_target(gtests_param parameterTest)
# add_custom_target(gtests_meas measurement_cvTest)
# add_custom_target(gtests_factory measurement_factoryTest)
# add_custom_target(gtests_node nodeTest)
# add_custom_target(gtests_source_R2_posTest source_R2_posTest)