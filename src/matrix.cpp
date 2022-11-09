#include "matrix.hpp"
#include <bits/stdc++.h>

linalg::Matrix::Matrix(std::vector<size_t> &I, std::vector<size_t> &J,
                       std::vector<double> &val, const size_t m, const size_t n)
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
  for (auto& element : _data) {
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
  auto it = std::find_if(_data.begin(), _data.end(),
                         [&](const auto &a) { return i == a.I && j == a.J; });

  // If found update to new value and if not add to the end of the vector
  if (it != _data.end()) {
    it->val = val;
  } else {
    linalg::coo new_value;
    new_value.I = i;
    new_value.J = j;
    new_value.val = val;
    _data.push_back(new_value);
  }
}

double linalg::Matrix::at(size_t i, size_t j) {
  // Find position i,j
  auto it = std::find_if(_data.cbegin(), _data.cend(),
                         [&](const auto a) { return i == a.I && j == a.J; });

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

/// --------------------------------------------------
//void linalg::Matrix::set(size_t i, size_t j, double val) {
// void linalg::Matrix::subMatrix(int mat[N][N], int temp[N][N], int p, int q, int n) {
void linalg::Matrix::subMatrix(Matrix &temp, int p, int q, int n) {
   int i = 0, j = 0;
   // filling the sub matrix
   for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {
         // skipping if the current row or column is not equal to the current
         // element row and column
         if (row != p && col != q) {
            temp.set(i,j++,this->at(row,col));
            // temp[i][j++] = mat[row][col];
            if (j == n - 1) {
               j = 0;
               i++;
            }
         }
      }
   }
}
//double linalg::Matrix::at(size_t i, size_t j) {
//int linalg::Matrix::determinantOfMatrix(int matrix[N][N], int n) {
// matrix can be found with _data 
// N = _size.first
// n = _size.first
double linalg::Matrix::determinantOfMatrix(Matrix &mat,int n) {
  // Assert  matrix is square
  assert(_size.second == _size.first );
  //Matrix mat(_size.first, _size.second);

   double determinant = 0;
   if (n == 1) {
      return this->at(0,0);
   }
  //  if (n == 2) {
  //     return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
  //  }
   //int temp[N][N];
   int sign = 1;
   Matrix temp(_size.first, _size.second);
   for (int i = 0; i < n; i++) {
      subMatrix(temp, 0, i, n);
      determinant += sign * this->at(0,i) * determinantOfMatrix(temp, n - 1);
      sign = -sign;
   }
   return determinant;
}
