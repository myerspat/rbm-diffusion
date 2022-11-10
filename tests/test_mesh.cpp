#include "rbm/mesh.hpp"
#include "unit_test_framework.hpp"

TEST(build_mesh) {
  // Setting bounds
  std::pair<double, double> left_bound;
  left_bound.first = 0.0;
  left_bound.second = 1.0;
  std::pair<double, double> right_bound;
  right_bound.first = 1.0;
  right_bound.second = -2.0;
  // creating mesh and running
  Mesh my_mesh(80, 30, left_bound, right_bound);
  my_mesh.run(0.11, 0.10, 2.0);
  linalg::Matrix fission_mat = my_mesh.getF();
  linalg::Matrix migration_mat = my_mesh.getM();
  // Checking Results of middle diag values
  for (int i = 1; i < 80 - 1; i++) {
    ASSERT_ALMOST_EQUAL(fission_mat.at(i, i), 0.04125, .00001);
    ASSERT_ALMOST_EQUAL(migration_mat.at(i, i - 1), -5.33333333, 0.00001)
    ASSERT_ALMOST_EQUAL(migration_mat.at(i, i + 1), -5.33333333, 0.00001)
    ASSERT_ALMOST_EQUAL(migration_mat.at(i, i), 10.70416667, 0.00001)
  }
  // Testing for [79][78] and [79][79] positions
  ASSERT_ALMOST_EQUAL(migration_mat.at(79, 79), 5.84844527, 0.8);
  ASSERT_ALMOST_EQUAL(migration_mat.at(79, 78), -5.3333333, 0.8);
  // Testing for [0][0] and [0][1]
  ASSERT_ALMOST_EQUAL(migration_mat.at(0, 0), 5.3708333, 0.8);
  ASSERT_ALMOST_EQUAL(migration_mat.at(0, 1), -5.3333333, 0.8);
}

TEST_MAIN();
