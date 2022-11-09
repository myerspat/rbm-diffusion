#ifndef _MESH_
#define _MESH_
#include "matrix.hpp"
#include <utility>
#include <vector>

class Mesh {
private:
  linalg::Matrix _F;
  linalg::Matrix _M;
  int _N;           // Finite sizes of Mesh
  double _distance; // Length of Mesh
  std::pair<double, double> _left_bound;
  std::pair<double, double> _right_bound;
  // double _absorptionXS;
  // double _D;
  // std::pair<int, int> _left_boundary_condition;
  // std::pair<int, int> _right_boundary_condition;
  //
public:
  Mesh(int N, double distance, std::pair<double, double> left_bound,
       std::pair<double, double> right_bound) {
    _N = N;
    _distance = distance;
    _left_bound = left_bound;
    _right_bound = right_bound;
    linalg::Matrix mat(_N, _N);
    _F = mat;
    _M = mat;
  }

  void run(double fissionXS, double absorptionXS, double D);
  linalg::Matrix getF();
  linalg::Matrix getM();

  // ~Mesh();
};

#endif // !_MESH_
