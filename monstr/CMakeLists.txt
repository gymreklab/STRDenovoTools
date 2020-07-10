file(GLOB SOURCES *.cpp *.h)
add_library(monstr ${SOURCES})
target_link_libraries(monstr ${htslib_static} ${zlib_static})
target_include_directories(monstr PUBLIC "${PROJECT_BINARY_DIR}")
