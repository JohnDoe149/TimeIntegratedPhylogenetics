include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(
  test_runner
  ${PROJECT_SOURCE_DIR}/testing/matrix_math_test.cc
  ${PROJECT_SOURCE_DIR}/testing/random_tree_test.cc
  ${PROJECT_SOURCE_DIR}/src/modeling/model/RandomTree.cpp
  ${PROJECT_SOURCE_DIR}/src/modeling/model/TransitionProbability.cpp
  ${PROJECT_SOURCE_DIR}/src/modeling/parameters/trees/TreeObject.cpp
  ${PROJECT_SOURCE_DIR}/src/modeling/parameters/trees/Node.cpp
  ${PROJECT_SOURCE_DIR}/src/modeling/parameters/RateMatrix.cpp 
  ${PROJECT_SOURCE_DIR}/src/core/RandomVariable.cpp
  ${PROJECT_SOURCE_DIR}/src/core/Probability.cpp
  ${PROJECT_SOURCE_DIR}/src/core/Msg.cpp
  ${PROJECT_SOURCE_DIR}/src/core/EigenSystem.cpp
)
target_link_libraries(
  test_runner
  GTest::gtest_main
  ncl_static
)
include_directories(${PROJECT_SOURCE_DIR}/src)

include(GoogleTest)
gtest_discover_tests(test_runner)