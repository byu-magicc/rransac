cmake_minimum_required(VERSION 3.8)
project(rransac VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_BUILD_TYPE RELEASE)
# set(CMAKE_BUILD_TYPE DEBUG)

IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(Ceres REQUIRED)
find_package(OpenCV 4 REQUIRED )

IF(CMAKE_BUILD_TYPE MATCHES Debug)
message(STATUS "Build type is Debug")
message(STATUS ${CMAKE_BUILD_TYPE})
add_compile_definitions(DEBUG_BUILD)
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)



# message(STATUS "${EIGEN3_INCLUDE_DIR}")
# message(STATUS "here")

# include_directories(
#     ${EIGEN3_INCLUDE_DIR}
#     ${Boost_INCLUDE_DIRS}
#     ${OpenCV_INCLUDE_DIRS}
#     include/rransac
#     lie_groups/liegroups/include
#     )







# ########################
# ## Install ##
# ########################

set(lie_groups_header_path "${PROJECT_SOURCE_DIR}/lie_groups/liegroups/include")

set(lie_group_header_files 
${lie_groups_header_path}/lie_algebras/algebra_base.h
${lie_groups_header_path}/lie_algebras/rn.h 
${lie_groups_header_path}/lie_algebras/se2.h 
${lie_groups_header_path}/lie_algebras/se3.h 
${lie_groups_header_path}/lie_algebras/so2.h 
${lie_groups_header_path}/lie_algebras/so3.h 
${lie_groups_header_path}/lie_groups/group_base.h 
${lie_groups_header_path}/lie_groups/Rn.h 
${lie_groups_header_path}/lie_groups/SE2.h 
${lie_groups_header_path}/lie_groups/SE3.h 
${lie_groups_header_path}/lie_groups/SO2.h 
${lie_groups_header_path}/lie_groups/SO3.h 
${lie_groups_header_path}/state.h 
${lie_groups_header_path}/utilities.h
)


set(rransac_header_path "${PROJECT_SOURCE_DIR}/include")

set(rransac_header_files 
${rransac_header_path}/rransac/common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h
${rransac_header_path}/rransac/common/data_association/model_policies/model_pdf_policy.h
${rransac_header_path}/rransac/common/data_association/data_association_host.h 
${rransac_header_path}/rransac/common/measurement/measurement_base.h 
${rransac_header_path}/rransac/common/models/model_base.h 
${rransac_header_path}/rransac/common/models/model_RN.h 
${rransac_header_path}/rransac/common/models/model_SEN_pos_vel.h 
${rransac_header_path}/rransac/common/models/model_SEN_pose_twist.h 
${rransac_header_path}/rransac/common/models/model_manager.h 
${rransac_header_path}/rransac/common/sources/source_base.h 
${rransac_header_path}/rransac/common/sources/source_RN.h 
${rransac_header_path}/rransac/common/sources/source_SEN_pos_vel.h 
${rransac_header_path}/rransac/common/sources/source_SEN_pose_twist.h 
${rransac_header_path}/rransac/common/transformations/transformation_base.h
${rransac_header_path}/rransac/common/transformations/trans_homography.h 
${rransac_header_path}/rransac/common/transformations/transformation_null.h 
${rransac_header_path}/rransac/common/utilities.h 
${rransac_header_path}/rransac/data_containers/data_tree/data_tree_base.h 
${rransac_header_path}/rransac/data_containers/data_tree/data_tree_cluster.h 
${rransac_header_path}/rransac/data_containers/cluster.h 
${rransac_header_path}/rransac/data_containers/consensus_set.h 
${rransac_header_path}/rransac/track_initialization/lmle_policies/linear_lmle_policy.h 
${rransac_header_path}/rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h 
${rransac_header_path}/rransac/track_initialization/seed_policies/null_seed_policy.h 
${rransac_header_path}/rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h 
${rransac_header_path}/rransac/track_initialization/ransac.h 
${rransac_header_path}/rransac/parameters.h 
${rransac_header_path}/rransac/rransac.h 
${rransac_header_path}/rransac/system.h 
)

set(headers 
${lie_group_header_files}
${rransac_header_files})

set(sources 
src/parameters.cpp)

add_library(${PROJECT_NAME}_lib ${sources} ${headers})

set(include_dest "include/${PROJECT_NAME}-${PROJECT_VERSION}")
set(main_lib_dest "lib/${PROJECT_NAME}-${PROJECT_VERSION}")

install(TARGETS ${PROJECT_NAME}_lib EXPORT ${PROJECT_NAME}_lib DESTINATION ${main_lib_dest})
install(FILES ${headers} DESTINATION ${include_dest})

install(EXPORT ${PROJECT_NAME}_lib DESTINATION ${main_lib_dest})

target_include_directories(${PROJECT_NAME}_lib PUBLIC
$<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
$<INSTALL_INTERFACE:${EIGEN3_INCLUDE_DIR}> 
$<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
$<INSTALL_INTERFACE:${Boost_INCLUDE_DIRS}> 
$<BUILD_INTERFACE:${OpenCV_INCLUDE_DIRS}>
$<INSTALL_INTERFACE:${OpenCV_INCLUDE_DIRS}> 
$<BUILD_INTERFACE:${rransac_header_path}>
$<BUILD_INTERFACE:${lie_groups_header_path}>
$<INSTALL_INTERFACE:${include_dest}> 
$<INSTALL_INTERFACE:${main_lib_dest}> 
)

get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/rransac_lib.cmake)

#############
## Testing ##
#############

# add_executable(main src/main.cpp
#  src/parameters.cpp)
# target_link_libraries(main ${OpenCV_LIBS})
# ${Boost_LIBRARIES})
# src/common/sources/source_base.cpp)


# option(BUILD_TESTS "Build the Params GTests" OFF)
# if(BUILD_TESTS)
  # add_subdirectory(test)
# endif()


# option(BUILD_SIMS "Build the simulations" OFF)
# if(BUILD_SIMS)
  # add_subdirectory(simulations)
# endif()