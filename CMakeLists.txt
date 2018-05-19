project(Path_Planning)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp src/pathplanner.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows") 

find_path(UWS_INCLUDE_DIR uWS/uWS.h)
find_library(UWS_LIBRARY uWS)
include_directories(${UWS_INCLUDE_DIR})
link_libraries(${UWS_LIBRARY})

find_path(ZLIB_INCLUDE_DIR zlib.h)
find_library(ZLIB_LIBRARY zlib)
include_directories(${ZLIB_INCLUDE_DIR})
link_libraries(${ZLIB_LIBRARY}) 

endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

add_executable(path_planning ${sources})

target_link_libraries(path_planning z ssl uv uWS)
