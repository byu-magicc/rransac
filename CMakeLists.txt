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











# ########################
# ## Install ##
# ########################


set(lie_groups_header_path "${PROJECT_SOURCE_DIR}/include")

set(rransac_header_path "${PROJECT_SOURCE_DIR}/include")

set(sources 
src/parameters.cpp)


set(INSTALL_INCLUDE_DIR "include/${PROJECT_NAME}-${PROJECT_VERSION}" CACHE PATH "Installation directory for header files.")
set(INSTALL_LIB_DIR "lib/${PROJECT_NAME}-${PROJECT_VERSION}" CACHE PATH "Installation directory for libraries. ")
set(INSTALL_BIN_DIR "bin/${PROJECT_NAME}-${PROJECT_VERSION}" CACHE PATH "Installation directory for executable files.") 
set(INSTALL_CMAKE_DIR "lib/cmake/${PROJECT_NAME}-${PROJECT_VERSION}" CACHE PATH "Installation directory for cmake files.")


# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()


add_library(${PROJECT_NAME} SHARED ${sources} )

target_include_directories(${PROJECT_NAME} PUBLIC
$<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
$<INSTALL_INTERFACE:${EIGEN3_INCLUDE_DIR}> 
$<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
$<INSTALL_INTERFACE:${Boost_INCLUDE_DIRS}> 
$<BUILD_INTERFACE:${OpenCV_INCLUDE_DIRS}>
$<INSTALL_INTERFACE:${OpenCV_INCLUDE_DIRS}> 
$<BUILD_INTERFACE:${rransac_header_path}>
$<BUILD_INTERFACE:${lie_groups_header_path}>
$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}> 
$<INSTALL_INTERFACE:${INSTALL_LIB_DIR}> 
)


install(TARGETS ${PROJECT_NAME}
EXPORT ${PROJECT_NAME}
LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
RUNTIME DESTINATION "${INSTALL_BIN_DIR}")

install(EXPORT ${PROJECT_NAME} DESTINATION ${INSTALL_CMAKE_DIR})


install(DIRECTORY include/ DESTINATION "${INSTALL_INCLUDE_DIR}")
install(DIRECTORY liegroups/include/ DESTINATION "${INSTALL_INCLUDE_DIR}")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)



##################################
#   CMake Package Installation   #
##################################

# Create the rransac-config.cmake and rransac-config-version.cmake files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
configure_file("${CMAKE_MODULE_PATH}/rransac-config.cmake.in" "${PROJECT_BINARY_DIR}/rransac-config.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${RRANSAC_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file("${CMAKE_MODULE_PATH}/rransac-config.cmake.in" "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/rransac-config.cmake" @ONLY)
# ... for both
configure_file("${CMAKE_MODULE_PATH}/rransac-config-version.cmake.in" "${PROJECT_BINARY_DIR}/rransac-config-version.cmake" @ONLY)

# Install the rransac-config.cmake and rransac-config-version.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/rransac-config.cmake"
  "${PROJECT_BINARY_DIR}/rransac-config-version.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")


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