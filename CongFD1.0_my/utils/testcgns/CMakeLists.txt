
# set the project name and version
project(cgnstest VERSION 0.0.0)
add_executable(cgnstest write_grid_str.c)
#target_link_libraries(cgnstest cgns)
target_link_libraries(cgnstest libcongFD)
#target_include_directories(cgnstest PUBLIC ../../src/include)

