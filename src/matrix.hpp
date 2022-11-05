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
  std::vector<coo> _data;          // vector of I, J, val
  std::pair<size_t, size_t> _size; // Holds the size of the matrix (MxN)

public:
  // Constructors / Destructor
  Matrix(const size_t m, const size_t n) : _size(std::make_pair(m, n)){};
  Matrix(std::vector<size_t> &I, std::vector<size_t> &J,
         std::vector<double> &val, const size_t m, const size_t n);

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

  // Get the size of the matrix
  std::pair<size_t, size_t> size() { return _size; };

  // Scalar Operations
  Matrix operator*(double);
};

} // namespace linalg

#endif // _MATRIX_
