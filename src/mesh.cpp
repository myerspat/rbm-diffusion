#include "./mesh.hpp"
#include "matrix.hpp"
#include <bits/stdc++.h>

void Mesh::run(double fissionXS, double absorptionXS, double D) {
  //--------------------------
  // Checking boundary Conditions
  // Stores a/b in a pair (a,b)
  // -------------------------
  //
  assert(_left_bound.first != 0.0 || _right_bound.second != 0.0);
  double rL, rR, Do, Di, Dn;
  bool Do_flag, Dn_flag;
  double dx = _distance / static_cast<double>(_N);
  // Left boundary conditions: branch statements setting Do tilda
  // If a=0  b!= 0  then we get 1/ (0/1) = inf so Do = 0
  // If a!=0 b = 0 then 1/ (1/0) = 0 so 1/rL = 0
  // else we do normal conditions
  if (_left_bound.first == 0.0) {
    Do = 0;
  }
  else if (_left_bound.second == 0.0) {
    Do = 1.0 / (dx / (2.0 * D));
  }
  else {
    Do = 1 / (dx / (2.0 * D) - _left_bound.second / _left_bound.first);
  }
  // Right boundary conditions: branch statements
  if (_right_bound.first == 0.0) {
    Dn = 0; // This term is 1 / r
  }
  else if (_right_bound.second == 0.0) {
    Dn = 1 / (dx / (2 * D));
  }
  else {
    Dn = 1 / (dx / (2 * D) - _right_bound.second / _right_bound.first);
  }
  // For loop iterating of mesh elements
  for (int i = 0; i < _N; i++) {
    //----------------------
    // Creating Migration Matrix
    //  if left boundary case
    if (i == 0) {
      _M.set(i, i + 1, -1.0 / (dx / (2.0 * D) + dx / (2.0 * D)));
      _M.set(i, i, absorptionXS * dx + Do - _M.at(i, i + 1));
      // right boundary case
    }
    else if (i == (_N - 1)) {
      Di = 1.0 / (dx / (2.0 * D) + dx / (2.0 * D));
      _M.set(i, i - 1, -Di);
      _M.set(i, i, absorptionXS * dx + Dn - _M.at(i, i - 1));
      // General Case
    }
    else {
      _M.set(i, i + 1, -1.0 / (dx / (2.0 * D) + dx / (2.0 * D)));
      _M.set(i, i - 1, -1.0 / (dx / (2.0 * D) + dx / (2.0 * D)));
      _M.set(i, i, absorptionXS * dx - _M.at(i, i - 1) - _M.at(i, i + 1));
    }
    //----------------------
    // Creating Fission Matrix
    //-----------------------
    _F.set(i, i, fissionXS * dx);
  }
}
// --------------------------
// Method for outside class information to access the Fission Matrix, _F.
// --------------------------
linalg::Matrix Mesh::getF() { return _F; }
// --------------------------
// Method for ouside class information to access the Migration Matrix, _M.
linalg::Matrix Mesh::getM() { return _M; }

// Mesh::~Mesh() {}
