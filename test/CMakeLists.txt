cmake_minimum_required(VERSION 2.8.12)

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
# Source tests
#---------------------------------------
add_executable(source_baseTest
common/sources/source_base_test.cpp
../src/parameters.cpp)
add_dependencies(source_baseTest lie_groups)
target_link_libraries(source_baseTest gtest_main lie_groups)

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

# add_executable(rransacTest
# rransac_test.cpp
# ../src/parameters.cpp
# ../src/rransac.cpp)
# target_link_libraries(rransacTest gtest_main lie_groups)







# add_executable(consensus_setTest 
# data_structures/consensus_set_test.cpp)
# target_link_libraries(consensus_setTest gtest_main)



# add_executable(measurement_cvTest
# common/measurement/measurement_cv_test.cpp
# ../src/common/measurement/measurement_cv.cpp
# ../src/common/parameters.cpp)
# target_link_libraries(measurement_cvTest gtest_main)

# add_executable(measurement_factoryTest
# common/measurement/measurement_factory_test.cpp
# ../src/common/measurement/measurement_factory.cpp
# ../src/common/measurement/measurement_cv.cpp)
# target_link_libraries(measurement_factoryTest gtest_main)

# add_executable(nodeTest
# data_structures/node_test.cpp
# ../src/data_structures/node.cpp)
# target_link_libraries(nodeTest gtest_main)

# add_custom_target(gtests_param parameterTest)
# add_custom_target(gtests_meas measurement_cvTest)
# add_custom_target(gtests_factory measurement_factoryTest)
# add_custom_target(gtests_node nodeTest)
# add_custom_target(gtests_source_R2_posTest source_R2_posTest)