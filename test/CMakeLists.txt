# -----------------------------------------------------------------------------
# TEST 1
# -----------------------------------------------------------------------------
# Add the test.
add_executable(test_fixedpoint test_fixedpoint.cpp)
# Liking for the test.
target_link_libraries(test_fixedpoint fplib)
# Add the test.
add_test(test_1 test_fixedpoint)
