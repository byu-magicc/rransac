cmake_minimum_required(VERSION 3.1)
project(rransac)

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3)

include_directories(
    ${EIGEN3_INCLUDE_DIR}
    include/rransac)



#############
## Testing ##
#############


# option(BUILD_TESTS "Build the Params GTests" OFF)
# if(BUILD_TESTS)
  add_subdirectory(test)
# endif()
