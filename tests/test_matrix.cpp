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

  int val_solution = mat.determinantOfMatrix(mat, 3);

  std::cout << "determinant_test val is: " << val_solution;

  ASSERT_ALMOST_EQUAL(val_solution, correct_val, 0.00001);

  
  
}

TEST_MAIN();
