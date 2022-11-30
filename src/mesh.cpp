#include "rbm/mesh.hpp"
#include "rbm/meshElement.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xslice.hpp"
#include "xtensor/xview.hpp"
#include <bits/stdc++.h>
#include <cstdio>
#include <xtensor/xarray.hpp>

namespace mesh {

void Mesh::constructMesh(const std::vector<MeshElement>& elements)
{
  // Construct course grid
  _course_grid = constructCourseGrid(elements);

  // Construct fine grid
  _fine_grid = constructFineGrid(_course_grid);
}

xt::xarray<MeshElement> Mesh::constructCourseGrid(
  const std::vector<MeshElement>& elements)
{
  // Allocate space
  xt::xarray<MeshElement> course_grid({_yN_course, _xN_course});

  // Iterate through elements
  for (const auto& element : elements) {
    // Iterate through row and column indicies and assign respective elements
    for (size_t i = element.getRowIdx().first; i <= element.getRowIdx().second;
         i++) {
      for (size_t j = element.getColIdx().first;
           j <= element.getColIdx().second; j++) {
        course_grid(i, j) = element;
      }
    }
  }

  // Check the adjacent lengths of MeshElements are the same
  assert(checkSharedLengths(course_grid));

  return course_grid;
}

xt::xarray<MeshElement> Mesh::constructFineGrid(
  const xt::xarray<MeshElement>& course_grid)
{

  // Allocate space
  xt::xarray<MeshElement> fine_grid({getYN(), getXN()});

  // Initialize position on y-axis of fine grid
  size_t fine_i = 0;

  for (size_t course_i = 0; course_i < course_grid.shape(0); course_i++) {
    // Find the next y position for each time we go up a row in course grid
    size_t next_fine_i = fine_i + _yN_fine;

    // Initialize position on x-axis of fine grid
    size_t fine_j = 0;

    // For each MeshElement along the x-axis at row course_i
    for (size_t course_j = 0; course_j < course_grid.shape(1); course_j++) {
      // Get element at course_i and course_j
      const MeshElement& element = course_grid(course_i, course_j);

      // Find the next x position in the fine grid given the elements number of
      // bins
      size_t next_fine_j = fine_j + _xN_fine;

      // View the square in the fine_grid that corresponds to the course_grid
      auto positions = xt::view(fine_grid, xt::range(fine_i, next_fine_i),
        xt::range(fine_j, next_fine_j));

      // Assign those positions to element
      positions = element;

      // Increment the x position
      fine_j = next_fine_j;
    }

    // Increment the y position
    fine_i = next_fine_i;
  }

  return fine_grid;
}

bool Mesh::checkSharedLengths(const xt::xarray<MeshElement>& course_grid)
{

  // Checking if ly is correct in each row
  for (size_t i = 0; i < course_grid.shape(0); i++) {
    auto row = xt::row(course_grid, i);
    for (int j = 1; j < course_grid.shape(0); j++) {
      if (row(j - 1).getLY() != row(j).getLY()) {
        return false;
      }
    }
  }
  // checking if lx is correct in each column
  for (size_t i = 0; i < course_grid.shape(1); i++) {
    auto col = xt::col(course_grid, i);
    for (int j = 1; j < course_grid.shape(1); j++) {
      if (col(j - 1).getLX() != col(j).getLX()) {
        return false;
      }
    }
  }
  // if each column and row have the same lx and ly respectively the check is
  // complete
  return true;
}

void Mesh::changeMaterial(const std::size_t& id, const double& new_value,
  const rbm::Parameter& target_parameter)
{
  // Iterating over the course mesh and changingg the value needed for a
  // specific id
  for (size_t i = 0; i < _course_grid.shape(0); i++) {
    for (size_t j = 0; j < _course_grid.shape(1); j++) {
      if (_course_grid.at(i, j).getID() == id) {
        _course_grid.at(i, j).setParameter(new_value, target_parameter);
      }
    }
  }
  // After reconstructing the course grid the fine mesh must be reconstructed.
  _fine_grid = constructFineGrid(_course_grid);
}

xt::xarray<double> Mesh::constructF()
{
  // Allocate space for fission operator
  std::vector<size_t> mesh_shape = {getSize(), getSize()};
  xt::xarray<double> F(mesh_shape);

  // Fill diagonal array
  for (size_t i = 0; i < _fine_grid.shape(0); i++) {
    for (size_t j = 0; j < _fine_grid.shape(1); j++) {
      // Convert 2D idx to 1D idx assuming row major ordering
      size_t position = ravelIDX(i, j);
      // Fill diagonal
      F(position, position) = _fine_grid(i, j).getMaterial().getNuFission();
    }
  }

  return F;
}

size_t Mesh::ravelIDX(const size_t& i, const size_t& j)
{
  return j + i * getXN();
}

xt::xarray<double> Mesh::constructM()
{
  // Allocate space
  std::vector<size_t> mesh_shape = {getSize(), getSize()};
  xt::xarray<double> M(mesh_shape);

  return M;
}

} // namespace mesh
