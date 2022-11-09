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
