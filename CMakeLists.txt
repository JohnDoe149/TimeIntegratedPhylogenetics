cmake_minimum_required(VERSION 3.10)
project(JohnPhylo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(${PROJECT_SOURCE_DIR}/lib/ncl)
add_subdirectory(${PROJECT_SOURCE_DIR}/lib/taskflow)

include_directories(${PROJECT_SOURCE_DIR}/src)

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(main ${SOURCES})

target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/lib/taskflow)

target_link_libraries(main PRIVATE ncl_static)

enable_testing()

add_executable(
  hello_test
  ${PROJECT_SOURCE_DIR}/src/hello_test.cc
)
target_link_libraries(
  hello_test
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(hello_test)

# comment out for debugging
#set(CMAKE_CXX_FLAGS, "${CMAKE_CXX_FLAGS} -march=native -O3")