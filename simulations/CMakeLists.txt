cmake_minimum_required(VERSION 3.8)

include_directories(../include
					../liegroups/include
					${EIGEN3_INCLUDE_DIR}
					${Boost_INCLUDE_DIRS}
					${OpenCV_INCLUDE_DIRS})

add_executable(sim_se2_pos_vel
se2_pos_vel_sim.cpp
../src/parameters.cpp)
target_link_libraries(sim_se2_pos_vel Ceres::ceres ${OpenCV_LIBS} )


add_executable(r2_tracking_se2_target_sim
r2_tracking_se2_target_sim.cpp
../src/parameters.cpp)
target_link_libraries(r2_tracking_se2_target_sim Ceres::ceres ${OpenCV_LIBS} )