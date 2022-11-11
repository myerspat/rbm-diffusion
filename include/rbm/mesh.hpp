#ifndef _MESH_
#define _MESH_

#include <cassert>
#include <utility>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>

class Mesh {
private:
  // Private variables
  xt::xarray<double> _F;                  // Fission operator
  xt::xarray<double> _M;                  // Migration operator
  int _N;                                 // Finite sizes of Mesh
  double _distance;                       // Length of Mesh
  std::pair<double, double> _left_bound;  // Left boundary condition
  std::pair<double, double> _right_bound; // Right boundary condition

public:
  //-----------------
  // Base Constructor
  // Inputs:
  //   distance = Mesh length
  //   left_bound = left boundary condition (a,b)
  //   right_bound = right_boundary condition (a,b)
  // ----------------
  Mesh(int N, double distance, std::pair<double, double> left_bound,
    std::pair<double, double> right_bound)
  {
    // asserting a and b are both not equal to zero for both boundary conditions
    // This is not possible (I think)
    assert(!(left_bound.first == 0.0 && left_bound.second == 0.0));
    assert(!(right_bound.first == 0.0 && right_bound.second == 0.0));

    _N = N;
    _distance = distance;
    _left_bound = left_bound;
    _right_bound = right_bound;
    _F = xt::zeros<double>({_N, _N});
    _M = xt::zeros<double>({_N, _N});
  }

  //---------------
  // This function "runs" the mesh to populate the Migration "_M" and Fission
  // "_F" Matrix Inputs:
  //   fissionXS = the fission cross section
  //   absorptionXS = the absorption cross section
  //   D = diffusion coefficient
  // ----------------
  void run(double fissionXS, double absorptionXS, double D);

  //------------
  // returns the Fission Matrix
  //----------------
  xt::xarray<double> getF() { return _F; }

  //------------------
  // Returns the Migration Matrix
  // --------------------
  xt::xarray<double> getM() { return _M; };
};

#endif // !_MESH_
