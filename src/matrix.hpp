#ifndef _MATRIX_
#define _MATRIX_

#include <cstddef>
#include <utility>
#include <vector>

namespace linalg {

// Coordinate format matrix
struct coo {
  size_t I;
  size_t J;
  double val;
};

class Matrix {
private:
  // Data
  std::vector<coo> _data; // vector of I, J, val

public:
  // Constructors / Destructor
  Matrix(){};
  Matrix(std::vector<size_t> &I, std::vector<size_t> &J,
         std::vector<double> &val);

  // Methods
  // Matrix vector multiplication
  std::vector<double> matvec(const std::vector<double> &col_vec);

  // Vector matrix multiplication
  std::vector<double> vecmat(const std::vector<double> &row_vec);

  // Matris matrix multiplication
  Matrix matmat(Matrix &mat);

  // At i,j set the value to val
  void set(size_t i, size_t j, double val);

  // Get value at i,j
  double at(size_t i, size_t j);

  // Get ith row
  std::vector<double> getRow(const size_t i);

  // Get jth column
  std::vector<double> getColumn(const size_t j);

  // Get matrix inverse
  Matrix getInverse();
};

} // namespace linalg

#endif // _MATRIX_
