#include "matrix.hpp"
#include <bits/stdc++.h>
#include <stdexcept>

linalg::Matrix::Matrix(std::vector<size_t> &I,
                       std::vector<size_t> &J,
                       std::vector<double> &val,
                       const size_t m,
                       const size_t n)
    : _size(std::make_pair(m, n)) {
  // Assert the vectors are the same length
  assert(I.size() == J.size() && I.size() == val.size());

  _data.reserve(val.size());
  for (size_t i = 0; i < val.size(); i++) {
    coo new_value;
    new_value.I = I[i];
    new_value.J = J[i];
    new_value.val = val[i];
    _data[i] = new_value;
  }
}

linalg::Matrix linalg::Matrix::matmat(Matrix &mat) {
  Matrix result(_size.first, _size.second);
  return result;
}

std::vector<double> linalg::Matrix::matvec(const std::vector<double> &col_vec) {
  // Initialize result vector
  std::vector<double> result(_size.first, 0);

  // Assert cols of first mat == rows of vec
  assert(_size.second == col_vec.size());

  // For each element in _data iterate through the row of mat at element.I and col at element.J
  for (auto &element : _data) {
    result[element.I] += element.val * col_vec[element.J];
  }

  return result;
}

std::vector<double> linalg::Matrix::vecmat(const std::vector<double> &row_vec) {
  std::vector<double> result(_size.second, 0);

  // add code for vecmat function here

  return result;
}

void linalg::Matrix::set(size_t i, size_t j, double val) {
  // Find the iterator at i and j
  auto it = std::find_if(_data.begin(), _data.end(), [&](const auto &a) { return i == a.I && j == a.J; });

  // If found update to new value and if not add to the end of the vector
  if (it != _data.end()) {
    it->val = val;
  }
  else {
    linalg::coo new_value;
    new_value.I = i;
    new_value.J = j;
    new_value.val = val;
    _data.push_back(new_value);
  }
}

double linalg::Matrix::at(size_t i, size_t j) {
  // Find position i,j
  auto it = std::find_if(_data.cbegin(), _data.cend(), [&](const auto a) { return i == a.I && j == a.J; });

  // If that position exists return the value, if not return 0
  if (it != _data.cend()) {
    return it->val;
  }
  return 0.0;
}

std::vector<double> linalg::Matrix::getRow(const size_t i) {
  // Initialize row vector
  std::vector<double> row;

  // Grab all in row i
  for (auto it = _data.cbegin(); it < _data.cend(); it++) {
    if (it->I == i) {
      row.push_back(it->val);
    }
  }
  return row;
}

std::vector<double> linalg::Matrix::getColumn(const size_t j) {
  // Initialize row vector
  std::vector<double> col;

  // Grab all in row i
  for (auto it = _data.cbegin(); it < _data.cend(); it++) {
    if (it->J == j) {
      col.push_back(it->val);
    }
  }
  return col;
}

linalg::Matrix linalg::Matrix::getInverse() {
  Matrix result(_size.first, _size.second);
  double determinant = 0;
  determinant = this->determinantOfMatrix();
  if(determinant == 0){
    throw std::runtime_error("Matrix not invertable"); 
  }

  Matrix ADJ = this->adjoint();

  result= ADJ*(1/determinant);


  return result;
}

linalg::Matrix linalg::Matrix::operator*(double scalar) {
  Matrix mat(_size.first, _size.second);

  // Multiply each element by the scalar
  for (size_t i = 0; i < _data.size(); i++) {
    mat.set(_data[i].I, _data[i].J, _data[i].val * scalar);
  }

  return mat;
}

linalg::Matrix linalg::Matrix::subMatrix(int p, int q, int n) {
  // Initialize temp matrix
  Matrix temp(n, n);

  // For each element in the matrix
  for (const auto &element : _data) {
    if (element.I != p && element.J != q) {
      // Calculate new column position to the left of q it will 
      // remain the same to the right we must subtract one
      // because we are removing a column
      size_t new_col = element.J;
      if (element.J > q)
        new_col--;

      size_t new_row = element.I;
      if (element.I > p)
        new_row--;

      // Set the value of the temp matrix
      temp.set(new_row, new_col, element.val);

    }
  }
  return temp;
}

double linalg::Matrix::determinantOfMatrix() {
  // Assert  matrix is square
  assert(_size.second == _size.first);
  int n = _size.first;

  double determinant = 0;
  if (n == 1) {
    return this->at(0, 0);
  }
  else if (n == 2) {
    return (at(0, 0) * at(1, 1)) - (at(0, 1) * at(1, 0));
  }
  int sign = 1;
  for (int i = 0; i < n; i++) {
    Matrix temp = subMatrix(0, i, n - 1);
    determinant += sign * this->at(0, i) * temp.determinantOfMatrix();
    sign = -sign;
  }
  return determinant;
}

linalg::Matrix linalg::Matrix::adjoint() {
  // Initialize temp matrix
  Matrix ADJ(_size.first, _size.second);

  int n = _size.first;
  if (n == 1) {
    ADJ.set(0,0,1);
    return ADJ;
  }
  int sign = 1;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Matrix temp = subMatrix(i, j, n-1);
      sign = ((i + j) % 2 == 0) ? 1 : -1;
      ADJ.set(j,i, (sign) * temp.determinantOfMatrix());
      printf("%d %d\n", (int)i, (int)j);
    }
  }  
  return ADJ;
}
