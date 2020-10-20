add_executable(CommonTest CommonTest.cpp)
target_link_libraries(CommonTest lib gtest gmock_main)
add_test(NAME CommonTest COMMAND CommonTest)
