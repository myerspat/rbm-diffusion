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

TEST(matrix_vector_1) {
  // Initialize matrix
  linalg::Matrix mat(2, 3);
  std::vector<double> vet{{2}, {1}, {0}};

  mat.set(0, 0, 1);
  mat.set(0, 1, -1);
  mat.set(0, 2, 2);
  mat.set(1, 0, 0);
  mat.set(1, 1, -3);
  mat.set(1, 2, 1);

  std::vector<double> vet_correct{{1}, {-3}};
  int correct_size = 2;

  std::vector<double> vect_solution = mat.matvec(vet);

  ASSERT_EQUAL(vect_solution.size(), correct_size);

  // Assertions
  for (size_t i = 0; i < 2; i++) {
    ASSERT_ALMOST_EQUAL(vect_solution[i], vet_correct[i], 0.00001);
  }
}

TEST(determinant_test) {
  // Initialize matrix
  linalg::Matrix mat(3, 3);

  mat.set(0, 0, 2);
  mat.set(0, 1, 1);
  mat.set(0, 2, 3);
  mat.set(1, 0, 6);
  mat.set(1, 1, 5);
  mat.set(1, 2, 7);
  mat.set(2, 0, 4);
  mat.set(2, 1, 9);
  mat.set(2, 2, 8);

  int correct_val = 36;

  int val_solution = mat.determinantOfMatrix();

  ASSERT_ALMOST_EQUAL(val_solution, correct_val, 0.00001);
}

TEST(adjoint_test) {
  // Initialize matrix
  linalg::Matrix mat(4, 4);

  mat.set(0, 0, 5);
  mat.set(0, 1, -2);
  mat.set(0, 2, 2);
  mat.set(0, 3, 7);
  mat.set(1, 0, 1);
  mat.set(1, 1, 0);
  mat.set(1, 2, 0);
  mat.set(1, 3, 3);
  mat.set(2, 0, -3);
  mat.set(2, 1, 1);
  mat.set(2, 2, 5);
  mat.set(2, 3, 0);
  mat.set(3, 0, 3);
  mat.set(3, 1, -1);
  mat.set(3, 2, -9);
  mat.set(3, 3, 4);

  linalg::Matrix mat_correct(4, 4);
  mat_correct.set(0, 0, -12);
  mat_correct.set(0, 1, 76);
  mat_correct.set(0, 2, -60);
  mat_correct.set(0, 3, -36);
  mat_correct.set(1, 0, -56);
  mat_correct.set(1, 1, 208);
  mat_correct.set(1, 2, -82);
  mat_correct.set(1, 3, -58);
  mat_correct.set(2, 0, 4);
  mat_correct.set(2, 1, 4);
  mat_correct.set(2, 2, -2);
  mat_correct.set(2, 3, -10);
  mat_correct.set(3, 0, 4);
  mat_correct.set(3, 1, 4);
  mat_correct.set(3, 2, 20);
  mat_correct.set(3, 3, 12);

  linalg::Matrix mat_solution = mat.adjoint();

  ASSERT_EQUAL(mat_solution.at(0, 0), mat_correct.at(0, 0));
  ASSERT_EQUAL(mat_solution.at(0, 1), mat_correct.at(0, 1));
  ASSERT_EQUAL(mat_solution.at(0, 2), mat_correct.at(0, 2));
  ASSERT_EQUAL(mat_solution.at(0, 3), mat_correct.at(0, 3));
  ASSERT_EQUAL(mat_solution.at(1, 0), mat_correct.at(1, 0));
  ASSERT_EQUAL(mat_solution.at(1, 1), mat_correct.at(1, 1));
  ASSERT_EQUAL(mat_solution.at(1, 2), mat_correct.at(1, 2));
  ASSERT_EQUAL(mat_solution.at(1, 3), mat_correct.at(1, 3));
  ASSERT_EQUAL(mat_solution.at(2, 0), mat_correct.at(2, 0));
  ASSERT_EQUAL(mat_solution.at(2, 1), mat_correct.at(2, 1));
  ASSERT_EQUAL(mat_solution.at(2, 2), mat_correct.at(2, 2));
  ASSERT_EQUAL(mat_solution.at(2, 3), mat_correct.at(2, 3));
  ASSERT_EQUAL(mat_solution.at(3, 0), mat_correct.at(3, 0));
  ASSERT_EQUAL(mat_solution.at(3, 1), mat_correct.at(3, 1));
  ASSERT_EQUAL(mat_solution.at(3, 2), mat_correct.at(3, 2));
  ASSERT_EQUAL(mat_solution.at(3, 3), mat_correct.at(3, 3));

}

TEST(inverse_test) {
  // Initialize matrix
  linalg::Matrix mat(4, 4);

  mat.set(0, 0, 5);
  mat.set(0, 1, -2);
  mat.set(0, 2, 2);
  mat.set(0, 3, 7);
  mat.set(1, 0, 1);
  mat.set(1, 1, 0);
  mat.set(1, 2, 0);
  mat.set(1, 3, 3);
  mat.set(2, 0, -3);
  mat.set(2, 1, 1);
  mat.set(2, 2, 5);
  mat.set(2, 3, 0);
  mat.set(3, 0, 3);
  mat.set(3, 1, -1);
  mat.set(3, 2, -9);
  mat.set(3, 3, 4);

  linalg::Matrix mat_correct(4, 4);
  mat_correct.set(0, 0, -0.136364);
  mat_correct.set(0, 1, 0.863636);
  mat_correct.set(0, 2, -0.681818);
  mat_correct.set(0, 3, -0.409091);
  mat_correct.set(1, 0, -0.636364);
  mat_correct.set(1, 1, 2.36364);
  mat_correct.set(1, 2, -0.931818);
  mat_correct.set(1, 3, -0.659091);
  mat_correct.set(2, 0, 0.0454545);
  mat_correct.set(2, 1, 0.0454545);
  mat_correct.set(2, 2, -0.0227273);
  mat_correct.set(2, 3, -0.113636);
  mat_correct.set(3, 0, 0.0454545);
  mat_correct.set(3, 1, 0.0454545);
  mat_correct.set(3, 2, 0.227273);
  mat_correct.set(3, 3, 0.136364);

  linalg::Matrix mat_solution = mat.getInverse();

  ASSERT_ALMOST_EQUAL(mat_solution.at(0, 0), mat_correct.at(0, 0), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(0, 1), mat_correct.at(0, 1), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(0, 2), mat_correct.at(0, 2), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(0, 3), mat_correct.at(0, 3), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(1, 0), mat_correct.at(1, 0), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(1, 1), mat_correct.at(1, 1), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(1, 2), mat_correct.at(1, 2), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(1, 3), mat_correct.at(1, 3), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(2, 0), mat_correct.at(2, 0), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(2, 1), mat_correct.at(2, 1), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(2, 2), mat_correct.at(2, 2), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(2, 3), mat_correct.at(2, 3), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(3, 0), mat_correct.at(3, 0), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(3, 1), mat_correct.at(3, 1), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(3, 2), mat_correct.at(3, 2), 0.00001);
  ASSERT_ALMOST_EQUAL(mat_solution.at(3, 3), mat_correct.at(3, 3), 0.00001);
}

TEST(vector_matrix_1) {
  // Initialize matrix
  linalg::Matrix mat(3, 4);
  std::vector<double> vet{3, 4, 2};

  mat.set(0, 0, 13);
  mat.set(0, 1, 9);
  mat.set(0, 2, 7);
  mat.set(0, 3, 15);
  mat.set(1, 0, 8);
  mat.set(1, 1, 7);
  mat.set(1, 2, 4);
  mat.set(1, 3, 6);
  mat.set(2, 0, 6);
  mat.set(2, 1, 4);
  mat.set(2, 2, 0);
  mat.set(2, 3, 3);

  std::vector<double> vet_correct{83, 63, 37, 75};
  int correct_size = 4; 

  std::vector<double> vect_solution = mat.vecmat(vet) ;


  ASSERT_EQUAL(vect_solution.size(), correct_size);

  // Assertions
  for (size_t i = 0; i < correct_size; i++) {
    ASSERT_ALMOST_EQUAL(vect_solution[i], vet_correct[i], 0.00001);
  }
}

TEST_MAIN();
