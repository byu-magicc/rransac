cmake_minimum_required(VERSION 3.8)



add_executable(sim_se2_pos_vel
se2_pos_vel_sim.cpp
../src/parameters.cpp)

target_link_libraries(sim_se2_pos_vel Ceres::ceres ${OpenCV_LIBS} )