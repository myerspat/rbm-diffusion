#ifndef _MESH_
#define _MESH_
#include "matrix.hpp"
#include <cassert>
#include <utility>
#include <vector>

class Mesh {
  private:
    // Private variables
    linalg::Matrix _F;
    linalg::Matrix _M;
    int _N;           // Finite sizes of Mesh
    double _distance; // Length of Mesh
    std::pair<double, double> _left_bound;
    std::pair<double, double> _right_bound;

  public:
    //-----------------
    // Base Constructor
    // Inputs:
    //   distance = Mesh length
    //   left_bound = left boundary condition (a,b)
    //   right_bound = right_boundary condition (a,b)
    // ----------------
    Mesh(int N, double distance, std::pair<double, double> left_bound, std::pair<double, double> right_bound) {
      _N = N;
      _distance = distance;
      // asserting a and b are both not equal to zero for both boundary conditions
      // This is not possible (I think)
      assert(!(left_bound.first == 0.0 && left_bound.second == 0.0));
      assert(!(right_bound.first == 0.0 && right_bound.second == 0.0));
      _left_bound = left_bound;
      _right_bound = right_bound;
      linalg::Matrix mat(_N, _N);
      _F = mat;
      _M = mat;
    }

    //---------------
    // This function "runs" the mesh to populate the Migration "_M" and Fission "_F" Matrix
    // Inputs:
    //   fissionXS = the fission cross section
    //   absorptionXS = the absorption cross section
    //   D = diffusion coefficient
    // ----------------
    void run(double fissionXS, double absorptionXS, double D);

    //------------
    // returns the Fission Matrix
    //----------------
    linalg::Matrix getF() { return _F; }

    //------------------
    // Returns the Migration Matrix
    // --------------------
    linalg::Matrix getM() { return _M; };
};

#endif // !_MESH_
