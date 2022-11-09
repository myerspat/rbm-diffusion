#include "../src/linalg.hpp"
#include "../src/matrix.hpp"
#include "unit_test_framework.hpp"

TEST(power_iteration) {
  linalg::Matrix M(5, 5);
  linalg::Matrix F(5, 5);

  // Create test migration matrix
  M.set(0, 0, 0.93333);
  M.set(0, 1, -0.33333);
  M.set(1, 0, -0.33333);
  M.set(1, 1, 1.26667);
  M.set(1, 2, -0.33333);
  M.set(2, 1, -0.33333);
  M.set(2, 2, 1.26667);
  M.set(2, 3, -0.33333);
  M.set(3, 2, -0.33333);
  M.set(3, 3, 1.26667);
  M.set(3, 4, -0.33333);
  M.set(4, 3, -0.33333);
  M.set(4, 4, 1.21905);

  // Create test fission matrix
  for (size_t i = 0; i < 5; i++) {
    F.set(i, i, 0.66);
  }

  // Run power iteration
  std::pair<std::vector<double>, double> solution =
      linalg::powerIteration(M, F);

  // Expected results
  double correct_flux[5] = {1.37616, 1.27097, 1.06864, 0.78462, 0.44062};
  double correct_k = 1.05519;

  // Assertions
  for (size_t i = 0; i < 5; i++) {
    ASSERT_ALMOST_EQUAL(solution.first[i], correct_flux[i], 0.00001);
  }
  ASSERT_ALMOST_EQUAL(solution.second, correct_k, 0.00001);
}



TEST_MAIN();
