
file(GLOB_RECURSE SRCS CONFIGURE_DEPENDS src/*.cpp include/*.hpp)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_Release "${CMAKE_C_FLAGS_Release} -O3")

find_package(OpenMP)
add_library(libcongFD STATIC ${SRCS})
target_link_libraries(libcongFD PUBLIC cgns OpenMP::OpenMP_CXX)
target_include_directories(libcongFD PUBLIC include/)





#target_link_libraries(libcongFD OpenMP::OpenMP_CXX)

