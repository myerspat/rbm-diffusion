#ifndef _MESH_
#define _MESH_

#include "rbm/material.hpp"
#include "rbm/meshElement.hpp"
#include "rbm/parameter.hpp"
#include "rbm/region.hpp"
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
  xt::xarray<MeshElement> _fine_grid;   // Mesh centered elements
  xt::xarray<MeshElement> _course_grid; // Mesh centered elements
  std::vector<Region> _regions;
  xt::xarray<double> _D_matrix;
  std::pair<double, double> _left_bound;   // Left boundary condition
  std::pair<double, double> _right_bound;  // Right boundary condition
  std::pair<double, double> _top_bound;    // Top boundary condition
  std::pair<double, double> _bottom_bound; // Bottom boundary condition
  size_t _xN_fine;
  size_t _yN_fine;
  size_t _xN_course; // Number of mesh elements in the x direction
  size_t _yN_course; // Number of mesh elements in the y direction

public:
  //========================================================================
  // Constructors
  Mesh() {};
  Mesh(const size_t& xN_fine, const size_t& yN_fine, const size_t& xN_course,
    const size_t& yN_course, const std::pair<double, double>& left_bound,
    const std::pair<double, double>& right_bound,
    const std::pair<double, double>& top_bound,
    const std::pair<double, double>& bottom_bound)
    : _xN_fine(xN_fine), _yN_fine(yN_fine), _xN_course(xN_course),
      _yN_course(yN_course), _left_bound(left_bound), _right_bound(right_bound),
      _top_bound(top_bound), _bottom_bound(bottom_bound)
  {
    // Asserting a and b are both not equal to zero for all boundary conditions
    assert(!(left_bound.first == 0.0 && left_bound.second == 0.0));
    assert(!(right_bound.first == 0.0 && right_bound.second == 0.0));
    assert(!(top_bound.first == 0.0 && top_bound.second == 0.0));
    assert(!(bottom_bound.first == 0.0 && bottom_bound.second == 0.0));
  }

  //========================================================================
  // Methods
  // Construct mesh
  void constructMesh(const std::vector<MeshElement>& elements);

  // Construct course grid
  xt::xarray<MeshElement> constructCourseGrid(
    const std::vector<MeshElement>& elements);

  // Print course grid
  void printCourseGrid(const xt::xarray<MeshElement>& course_grid);

  // Construct fine grid
  xt::xarray<MeshElement> constructFineGrid(
    const xt::xarray<MeshElement>& course_grid);

  // Check lengths of adjacent elements match for all elements in _course_grid
  bool checkSharedLengths(const xt::xarray<MeshElement>& course_grid);

  // Change target mesh elements' material property
  void changeMaterial(const std::size_t& id, const double& new_value,
    const material::Property& target_property);

  xt::xarray<double> constructF(); // Construct the fission operator matrix
  xt::xarray<double> constructM(); // Construct the migration operator matrix
  xt::xarray<double>
  constructD(); // Construct the migration shell consisting of only D

  void updateD(std::vector<rbm::Parameter>& parameters, std::size_t& idx);

  // Assuming row major ordering, returns 1D idx given 2D idx
  size_t ravelIDX(const size_t& i, const size_t& j);

  //========================================================================
  // Getters
  // Get the number of bins along the x-axis
  size_t getXN() const { return _xN_fine * _xN_course; };

  // Get the number of bins along the y-axis
  size_t getYN() const { return _yN_fine * _yN_course; };

  // Get the number of total bins in the mesh
  size_t getSize() const { return getXN() * getYN(); };

  // Gets the course grid
  xt::xarray<MeshElement> getCourseGrid() const { return _course_grid; }

  // Gets the fine grid
  xt::xarray<MeshElement> getFineGrid() const { return _fine_grid; }

  // Gets the regions vector
  std::vector<Region> getRegions() const { return _regions; }

  // Gets the migration matrix shell
  xt::xarray<double> getDMatrix() const { return _D_matrix; }
};

} // namespace mesh

#endif // !_MESH_
