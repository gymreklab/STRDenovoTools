cmake_minimum_required(VERSION 3.7)

include(ExternalProject)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_compile_options(-Werror)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


# set the project name
project(MonSTR VERSION 2.0)
configure_file(MonSTRConfig.h.in MonSTRConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

ExternalProject_Add(zlib
        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/zlib
        GIT_REPOSITORY "https://github.com/madler/zlib.git"
        GIT_TAG "v1.2.8"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/thirdparty/zlib/src/zlib/configure --prefix=${CMAKE_BINARY_DIR}/thirdparty/zlib --static
        INSTALL_DIR ${CMAKE_BINARY_DIR}/thirdparty/zlib
        LOG_DOWNLOAD 1
        LOG_INSTALL 1
        )
ExternalProject_Add(htslib
        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/htslib
        GIT_REPOSITORY "https://github.com/samtools/htslib.git"
        GIT_TAG "1.3.1"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make
        INSTALL_COMMAND make install prefix=${CMAKE_BINARY_DIR}/thirdparty/htslib
        LOG_DOWNLOAD 1
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR}/thirdparty/zlib/include)
set(zlib_static ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib/libz.a)
set(htslib_static ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib/libhts.a)
include_directories(${CMAKE_BINARY_DIR}/thirdparty/htslib/include)
add_dependencies(htslib zlib)

add_subdirectory(monstr)

# add the executable
add_executable(MonSTR src/main_sdt.cpp)
# add_executable(MonSTRExamine src/main_examine.cpp)

target_compile_features(MonSTR PRIVATE cxx_range_for)
target_include_directories(MonSTR PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(MonSTR monstr pthread)
install(TARGETS MonSTR DESTINATION bin)

# target_compile_features(MonSTRExamine PRIVATE cxx_range_for)
# target_include_directories(MonSTRExamine PUBLIC "${PROJECT_BINARY_DIR}")
# target_link_libraries(MonSTRExamine monstr pthread)
# install(TARGETS MonSTRExamine DESTINATION bin)

add_dependencies(htslib zlib)
add_dependencies(monstr htslib)
