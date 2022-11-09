#include "./mesh.hpp"
#include "matrix.hpp"
void Mesh::run(double fissionXS, double absorptionXS, double D) {
  //--------------------------
  // Checking boundary Conditions
  // Stores a/b in a pair (a,b)
  // -------------------------
  double rL, rR, Do, Di, Dn;
  bool Do_flag, Dn_flag;
  // Left boundary conditions: branch statements setting Do tilda
  // If a or b = 0 then we get 1/0 = inf so D0 = 0
  // else we do normal conditions
  if (_left_bound.second == 0.0 || _left_bound.first == 0.0) {
    Do = 0;
    Do_flag = false;
  } else {
    rL = 1.0 / (_left_bound.first / _left_bound.second);
    Do_flag = true;
  }
  // Right boundary conditions: branch statements
  if (_right_bound.second == 0.0 || _right_bound.first == 0.0) {
    Dn = 0; // This term is 1 / r
    Dn_flag = false;
  } else {
    rR = 1.0 /
         (_right_bound.first / _right_bound.second); // this term is 1/ (a/b)
    Dn_flag = true;
  }

  double dx = _distance / static_cast<double>(_N);
  for (int i = 0; i < _N; i++) {
    //----------------------
    // Creating Migration Matrix
    //  if left boundary case
    if (i == 0) {
      if (Do_flag) {
        Do = 1.0 / (dx / (2 * D) + rL);
      }
      _M.set(i, i + 1, -1.0 / (dx / (2.0 * D) + dx / (2.0 * D)));
      _M.set(i, i, absorptionXS * dx + Do - _M.at(i, i + 1));
      // right boundary case
    } else if (i == (_N - 1)) {
      Di = 1.0 / (dx / (2.0 * D) + dx / (2.0 * D));
      if (Dn_flag) {
        Dn = 1.0 / (dx / (2.0 * D) + rR);
      }
      _M.set(i, i - 1, -Di);
      _M.set(i, i, absorptionXS * dx + Dn - _M.at(i, i - 1));

      // General Case
    } else {
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
