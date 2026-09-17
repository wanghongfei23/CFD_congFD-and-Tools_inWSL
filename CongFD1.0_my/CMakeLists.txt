# 指定构建此项目所需的最低CMake版本为3.15，必须放在文件第一行[1,6](@ref)
# 如果当前CMake版本低于此要求，配置过程将报错停止[7](@ref)
cmake_minimum_required(VERSION 3.15)

# 设置项目名称为congFD，版本号为0.0.0[1](@ref)
# project命令会定义变量PROJECT_NAME(值为congFD)和PROJECT_VERSION(值为0.0.0)[2](@ref)
project(congFD VERSION 0.0.0)

# 配置头文件：将模板文件helloConfig.h.in中的CMake变量替换为实际值，生成helloConfig.h[9](@ref)
# 这样可以在C++代码中使用项目配置信息（如版本号）[1](@ref)
configure_file(helloConfig.h.in helloConfig.h)

# 设置C++编译标准为C++23[1](@ref)
set(CMAKE_CXX_STANDARD 23)

# 这行被注释掉的代码原本用于生成时间戳字符串并存储在COMPILE_TIME变量中
#string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)

# 以下两行被注释的代码涉及添加子目录和库链接
#add_subdirectory(libMeshb)   # 添加libMeshb子目录，CMake会处理其中的CMakeLists.txt[1](@ref)
#list(APPEND EXTRA_LIBS libMeshb)  # 将libMeshb库添加到EXTRA_LIBS变量中

# 设置编译模式为Release（发布模式）[3](@ref)  # Release模式通常会进行更多优化，减小生成二进制文件的大小[2](@ref)
# set(CMAKE_BUILD_TYPE Release)

# 设置C++编译器标志：添加-g（调试信息）、-w（抑制所有警告）、-O3（三级优化，较-O2更激进）[2](@ref)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -w -O3")

# 设置Debug模式的C编译器标志：添加-O3优化级别[3](@ref)
# set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# 设置Release模式的C编译器标志：添加-O3优化级别[3](@ref)
set(CMAKE_C_FLAGS_Release "${CMAKE_CXX_FLAGS_Release} -O3")

# 查找OpenMP包，OpenMP用于支持并行编程[1](@ref)
# 如果找到OpenMP，会定义OpenMP::OpenMP_CXX等目标供后续链接使用
find_package(OpenMP)

# 创建名为${PROJECT_NAME}（即congFD）的可执行文件，源文件为zoneMain.cpp[1](@ref)
add_executable(${PROJECT_NAME} zoneMain.cpp)

# 为congFD目标链接库：libcongFD和OpenMP::OpenMP_CXX[1](@ref)
# PUBLIC关键字表示这些库的依赖关系会传递给链接congFD的其他目标[2](@ref)
target_link_libraries(${PROJECT_NAME} PUBLIC libcongFD OpenMP::OpenMP_CXX)

# 创建名为AccuracyTest的可执行文件，源文件为AccuracyTest.cpp
add_executable(AccuracyTest AccuracyTest.cpp)

# 为AccuracyTest目标链接相同的库
target_link_libraries(AccuracyTest PUBLIC libcongFD OpenMP::OpenMP_CXX)

# 创建名为symmetrytest的可执行文件，源文件为symmetrytest.cpp[1](@ref)
add_executable(symmetrytest symmetrytest.cpp)

# 为symmetrytest目标链接相同的库
target_link_libraries(symmetrytest PUBLIC libcongFD OpenMP::OpenMP_CXX)

# 创建名为dim3KernelTest的可执行文件，源文件为dim3KernelTest.cpp（三维物理层内核测试，任务 3 起）
add_executable(dim3KernelTest dim3KernelTest.cpp)

# 为dim3KernelTest目标链接相同的库
target_link_libraries(dim3KernelTest PUBLIC libcongFD OpenMP::OpenMP_CXX)

# 添加src和utils子目录，CMake会处理这些子目录中的CMakeLists.txt文件[1](@ref)
# 通常这些子目录中定义了libcongFD库的构建规则
add_subdirectory(src)
add_subdirectory(utils)

# 包含CMakePrintHelpers模块，该模块提供了cmake_print_properties等调试辅助命令[1](@ref)
include(CMakePrintHelpers)

# 打印OpenMP::OpenMP_CXX目标的属性，用于调试和验证OpenMP配置[1](@ref)
# 会显示编译选项、头文件目录和链接库信息
cmake_print_properties(
TARGETS
OpenMP::OpenMP_CXX
PROPERTIES
INTERFACE_COMPILE_OPTIONS      # 接口编译选项
INTERFACE_INCLUDE_DIRECTORIES  # 接口头文件目录
INTERFACE_LINK_LIBRARIES       # 接口链接库
)
