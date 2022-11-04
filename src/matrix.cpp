#include "matrix.hpp"
#include <bits/stdc++.h>

linalg::Matrix::Matrix(std::vector<size_t> &I, std::vector<size_t> &J,
                       std::vector<double> &val) {
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

// std::vector<double> linalg::Matrix::matvec(const std::vector<double> &col_vec) {
// }

// std::vector<double> linalg::Matrix::vecmat(const std::vector<double> &row_vec) {
// }

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

// linalg::Matrix linalg::Matrix::getInverse() {}
