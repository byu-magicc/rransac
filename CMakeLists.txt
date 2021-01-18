cmake_minimum_required(VERSION 3.8)
project(rransac)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_BUILD_TYPE RELEASE)
# set(CMAKE_BUILD_TYPE DEBUG)

IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

find_package(Eigen3)
find_package(Boost)
find_package(Ceres REQUIRED)

# message(STATUS "${EIGEN3_INCLUDE_DIR}")
# message(STATUS "here")

include_directories(
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    include/rransac
    lie_groups/liegroups/include
    include
    src)


add_library(lie_groups
lie_groups/liegroups/src/lie_algebras/so2.cpp
lie_groups/liegroups/src/lie_algebras/so3.cpp
lie_groups/liegroups/src/lie_algebras/se2.cpp
lie_groups/liegroups/src/lie_algebras/se3.cpp
lie_groups/liegroups/src/lie_groups/SO2.cpp
lie_groups/liegroups/src/lie_groups/SO3.cpp
lie_groups/liegroups/src/lie_groups/SE2.cpp
lie_groups/liegroups/src/lie_groups/SE3.cpp
)


# target_link_libraries(helloworld Ceres::ceres)

# ########################
# ## Lie Groups Library ##
# ########################
# include(ExternalProject)
# ExternalProject_Add (
#   lie_groups
#   GIT_REPOSITORY git@magiccvs.byu.edu:mark_petersen/liegroups.git
#   GIT_TAG master
# )

#############
## Testing ##
#############

# add_executable(main src/main.cpp
# src/parameters.cpp)
# ${Boost_LIBRARIES})
# src/common/sources/source_base.cpp)


# option(BUILD_TESTS "Build the Params GTests" OFF)
# if(BUILD_TESTS)
  add_subdirectory(test)
# endif()
