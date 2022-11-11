#include "rbm/mesh.hpp"
#include <bits/stdc++.h>

void Mesh::run(double fissionXS, double absorptionXS, double D)
{
  //--------------------------
  // Checking boundary Conditions
  // Stores a/b in a pair (a,b)
  // -------------------------
  //
  double Do, Dn;
  double dx = _distance / static_cast<double>(_N);

  // Left boundary conditions: branch statements setting Do tilda
  // If a=0  b!= 0  then we get 1/ (0/1) = inf so Do = 0
  // If a!=0 b = 0 then 1/ (1/0) = 0 so 1/rL = 0
  // else we do normal conditions
  if (_left_bound.first == 0.0) {
    Do = 0;
  } else if (_left_bound.second == 0.0) {
    Do = 1.0 / (dx / (2.0 * D));
  } else {
    Do = 1 / (dx / (2.0 * D) - _left_bound.second / _left_bound.first);
  }

  // Right boundary conditions: branch statements
  if (_right_bound.first == 0.0) {
    Dn = 0; // This term is 1 / r
  } else if (_right_bound.second == 0.0) {
    Dn = 1 / (dx / (2 * D));
  } else {
    Dn = 1 / (dx / (2 * D) - _right_bound.second / _right_bound.first);
  }

  // Left boundary (top row)
  _M(0, 1) = -1.0 / (dx / (2.0 * D) + dx / (2.0 * D));
  _M(0, 0) = absorptionXS * dx + Do - _M(0, 1);
  _F(0, 0) = fissionXS * dx;

  // Right boundary (bottom row)
  int right_idx = _N - 1;
  _M(right_idx, right_idx - 1) = -1.0 / (dx / (2.0 * D) + dx / (2.0 * D));
  _M(right_idx, right_idx) =
    absorptionXS * dx + Dn - _M(right_idx, right_idx - 1);
  _F(right_idx, right_idx) = fissionXS * dx;

  // For loop iterating of mesh elements
  for (int i = 1; i < _N - 1; i++) {
    // Build the rest of the migration operator
    _M(i, i - 1) = -1.0 / (dx / (2.0 * D) + dx / (2.0 * D));
    _M(i, i + 1) = -1.0 / (dx / (2.0 * D) + dx / (2.0 * D));
    _M(i, i) = absorptionXS * dx - _M(i, i - 1) - _M(i, i + 1);

    // Build fission operator
    _F(i, i) = fissionXS * dx;
  }
}
