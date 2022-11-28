#include "rbm/mesh.hpp"
#include "xtensor/xbuilder.hpp"
#include <bits/stdc++.h>

namespace mesh {

void Mesh::constructMesh(const std::vector<MeshElement>& elements) {
  // Construct course grid
  _course_grid = constructCourseGrid(elements);

  // Construct fine grid
  _fine_grid = constructFineGrid();
}

xt::xarray<MeshElement> Mesh::constructCourseGrid(
  const std::vector<MeshElement>& elements)
{
  // Allocate space
  xt::xarray<MeshElement> course_grid({_yN_course, _xN_course});

  return course_grid;
}

xt::xarray<MeshElement> Mesh::constructFineGrid()
{
  // Allocate space
  xt::xarray<MeshElement> fine_grid({_yN_fine, _xN_fine});

  return fine_grid;
}

bool Mesh::checkSharedLengths(const xt::xarray<MeshElement>& course_grid)
{
  bool shared_lengths = true;

  return shared_lengths;
}

void Mesh::changeMaterail(
  std::size_t id, double new_value, std::string& target_parameter)
{}

xt::xarray<double> Mesh::constructF()
{
  // Allocate space
  std::vector<size_t> mesh_shape = {getSize(), getSize()};
  xt::xarray<double> F(mesh_shape);

  return F;
}

xt::xarray<double> Mesh::constructM()
{
  // Allocate space
  std::vector<size_t> mesh_shape = {getSize(), getSize()};
  xt::xarray<double> M(mesh_shape);

  return M;
}

} // namespace mesh
