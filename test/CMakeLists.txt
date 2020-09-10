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

add_executable(parameterTest 
common/parameters_test.cpp
../src/common/parameters.cpp)
target_link_libraries(parameterTest gtest_main)

add_executable(measurement_cvTest
common/measurement/measurement_cv_test.cpp
../src/common/measurement/measurement_cv.cpp
../src/common/parameters.cpp)
target_link_libraries(measurement_cvTest gtest_main)

add_executable(measurement_factoryTest
common/measurement/measurement_factory_test.cpp
../src/common/measurement/measurement_factory.cpp
../src/common/measurement/measurement_cv.cpp)
target_link_libraries(measurement_factoryTest gtest_main)

add_custom_target(gtests_param parameterTest)
add_custom_target(gtests_meas measurement_cvTest)
add_custom_target(gtests_factory measurement_factoryTest)