#include "../src/mesh.cpp"
#include "../src/mesh.hpp"
#include "unit_test_framework.hpp"

TEST(build_mesh) {
  std::pair<double, double> left_bound;
  left_bound.first = 0;
  left_bound.second = 0;
  std::pair<double, double> right_bound;
  right_bound.first = 1.0;
  right_bound.second = 2.0;
  Mesh my_mesh(80, 30, left_bound, right_bound);
  my_mesh.run(0.11, 0.10, 2.0);
  linalg::Matrix fission_mat = my_mesh.getF();
  linalg::Matrix migration_mat = my_mesh.getM();

  for (int i = 0; i < 80; i++) {
    ASSERT_EQUAL(fission_mat.at(i, i), 0.04125);
    ASSERT_EQUAL(migration_mat.at(i, i - 1), -5.33333333)
    ASSERT_EQUAL(migration_mat.at(i, i + 1), -5.33333333)
    ASSERT_EQUAL(migration_mat.at(i, i), 10.70416667)
  }
}

TEST_MAIN();
