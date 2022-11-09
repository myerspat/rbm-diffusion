#include "unit_test_framework.hpp"

#include "../src/matrix.hpp"
#include <cstdio>

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

TEST(matrix_vector_based_constructor) {
  std::vector<size_t> I = {0, 0, 0, 1, 1, 1, 2, 2, 2};
  std::vector<size_t> J = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  std::vector<double> vals = {1.2, 3.1, 6.0, 5.5, 7.8, 9.5, 6.3, 3.5, 1.8};

  // Construct matrix
  linalg::Matrix mat(I, J, vals, 3, 3);

  // Check each value in the matrix
  for (size_t i : I) {
    for (size_t j : J) {
      ASSERT_EQUAL(mat.at(i, j), vals[i * 3 + j]);
    }
  }
}

TEST(test_matmat_multiply) {
  // Initialize two 3x3 matricies
  std::vector<size_t> I = {0, 0, 0, 1, 1, 1, 2, 2, 2};
  std::vector<size_t> J = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  std::vector<double> vals1 = {1.2, 3.1, 6.0, 5.5, 7.8, 9.5, 6.3, 3.5, 1.8};
  std::vector<double> vals2 = {0.0, 6.7, 4.1, 2.2, 0.3, 6.0, 7.1, 8.1, 3.3};

  linalg::Matrix mat1(I, J, vals1, 3, 3);
  linalg::Matrix mat2(I, J, vals2, 3, 3);

  double correct_mat[3][3] = {
      {49.42, 57.57, 43.32}, {84.61, 116.14, 100.70}, {20.48, 57.84, 52.77}};

  linalg::Matrix result = mat1.matmat(mat2);

  for (size_t i : I) {
    for (size_t j : J) {
      ASSERT_ALMOST_EQUAL(result.at(i, j), correct_mat[i][j], 0.000001);
    }
  }
}

TEST(test_sparse_matmat_multiply) {
  // Initialize two 3x3 matricies
  std::vector<size_t> I = {0, 0, 0, 1, 1, 1, 2, 2, 2};
  std::vector<size_t> J = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  std::vector<double> vals1 = {0.0, 3.1, 0.0, 5.5, 0.0, 9.5, 0.0, 3.5, 0.0};
  std::vector<double> vals2 = {0.0, 0.0, 4.1, 0.0, 0.0, 6.0, 7.1, 8.1, 3.3};

  linalg::Matrix mat1(I, J, vals1, 3, 3);
  linalg::Matrix mat2(I, J, vals2, 3, 3);

  double correct_mat[3][3] = {
      {0.0, 0.0, 18.6}, {67.45, 76.95, 53.90}, {0.0, 0.0, 21}};

  linalg::Matrix result = mat1.matmat(mat2);

  for (size_t i : I) {
    for (size_t j : J) {
      ASSERT_ALMOST_EQUAL(result.at(i, j), correct_mat[i][j], 0.000001);
    }
  }
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
