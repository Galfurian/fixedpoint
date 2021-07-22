# -----------------------------------------------------------------------------
# TEST 1
# -----------------------------------------------------------------------------
# Add the test.
add_executable(${PROJECT_NAME}_test_fixedpoint test_fixedpoint.cpp)
# Liking for the test.
target_link_libraries(${PROJECT_NAME}_test_fixedpoint fplib)
# Add the test.
add_test(${PROJECT_NAME}_test_1 ${PROJECT_NAME}_test_fixedpoint)
