cmake_minimum_required(VERSION 2.8.3)
project(ardrone_pack)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  sensor_msgs
  cv_bridge
  image_transport
  message_generation
)
find_package(OpenCV REQUIRED)

add_message_files(
  FILES
)

generate_messages(
  DEPENDENCIES
)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES ardrone_pack
#  CATKIN_DEPENDS other_catkin_pkg
#  DEPENDS system_lib
  CATKIN_DEPENDS message_runtime
)

include_directories(
  ${catkin_INCLUDE_DIRS}
  include
)

add_executable(yaw_est src/yaw_estimation.cpp)
target_link_libraries(yaw_est
  ${catkin_LIBRARIES}
)


