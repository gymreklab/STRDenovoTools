file(GLOB SOURCES *.cpp *.h)
add_library(lib ${SOURCES})
target_link_libraries(lib ${htslib_static} ${zlib_static})
target_include_directories(lib PUBLIC "${PROJECT_BINARY_DIR}")

add_subdirectory(tests)
