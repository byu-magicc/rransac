cmake_minimum_required(VERSION 3.8)



add_executable(sim_se2_pos_vel
se2_pos_vel_sim.cpp
../src/parameters.cpp)
target_link_libraries(sim_se2_pos_vel Ceres::ceres ${OpenCV_LIBS} )


add_executable(r2_tracking_se2_target_sim
r2_tracking_se2_target_sim.cpp
../src/parameters.cpp)
target_link_libraries(r2_tracking_se2_target_sim Ceres::ceres ${OpenCV_LIBS} )