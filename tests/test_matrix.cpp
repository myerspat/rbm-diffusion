#include "unit_test_framework.hpp"

#include "../src/matrix.hpp"


TEST(build_matrix_1) {
  // Initialize matrix
  linalg::Matrix mat(5, 8);

  mat.set(1, 7, 5.6);
  mat.set(4, 7, 5.5);
  mat.set(4, 3, 5.7);

  ASSERT_EQUAL(mat.at(1, 7), 5.6);
  ASSERT_EQUAL(mat.at(4, 7), 5.5);
  ASSERT_EQUAL(mat.at(4, 3), 5.7);
  ASSERT_EQUAL(mat.size().first, 5);
  ASSERT_EQUAL(mat.size().second, 8);
}

TEST_MAIN();
