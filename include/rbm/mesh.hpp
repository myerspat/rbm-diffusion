#ifndef _MESH_
#define _MESH_

#include "rbm/cell.hpp"
#include "rbm/meshElement.hpp"
#include <cassert>
#include <utility>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>

namespace mesh {

class Mesh {
private:
  //========================================================================
  // Data
  xt::xarray<MeshElement> _mesh;           // Mesh centered elements
  std::pair<double, double> _left_bound;   // Left boundary condition
  std::pair<double, double> _right_bound;  // Right boundary condition
  std::pair<double, double> _top_bound;    // Top boundary condition
  std::pair<double, double> _bottom_bound; // Bottom boundary condition
  size_t _xN; // Number of mesh elements in the x direction
  size_t _yN; // Number of mesh elements in the y direction

public:
  //========================================================================
  // Constructors
  Mesh() {};
  Mesh(size_t xN, size_t yN, std::pair<double, double>& left_bound,
    std::pair<double, double>& right_bound, std::pair<double, double>& top_bound,
    std::pair<double, double>& bottom_bound)
    : _xN(xN), _yN(yN), _left_bound(left_bound), _right_bound(right_bound),
      _top_bound(top_bound), _bottom_bound(bottom_bound)
  {
    // Asserting a and b are both not equal to zero for both boundary conditions
    assert(!(left_bound.first == 0.0 && left_bound.second == 0.0));
    assert(!(right_bound.first == 0.0 && right_bound.second == 0.0));
    assert(!(top_bound.first == 0.0 && top_bound.second == 0.0));
    assert(!(bottom_bound.first == 0.0 && bottom_bound.second == 0.0));
    _mesh = xt::xarray<MeshElement>::from_shape({_yN, _xN});
  }

  //========================================================================
  // Methods
  // Initialize x-axis
  void initXAxis(
    std::vector<double>& section_lengths, std::vector<size_t>& x_bins);

  // Initialize y-axis
  void initYAxis(
    std::vector<double>& section_lengths, std::vector<size_t>& y_bins);

  // Build a rectangular section of the mesh and fill respective mesh elements
  void buildCells(std::vector<Cell>& cells);

  // Change the target parameter to new_value in every MeshElement with cell_id
  void changeCell(int cell_id, std::string& target, double new_value);

  xt::xarray<double> constructF(); // Construct the fission operator matrix
  xt::xarray<double> constructM(); // Construct the migration operator matrix

  //========================================================================
  // Getters
  // Get the number of bins along the x-axis
  size_t getXN() const { return _xN; };

  // Get the number of bins along the y-axis
  size_t getYN() const { return _yN; };

  // Get the number of total bins in the mesh
  size_t getSize() { return _xN * _yN; };
};

} // namespace mesh

#endif // !_MESH_
