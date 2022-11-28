#include "rbm/mesh.hpp"
#include "rbm/meshElement.hpp"
#include "xtensor/xbuilder.hpp"
#include <bits/stdc++.h>
#include <cstdio>
#include <xtensor/xarray.hpp>

namespace mesh {

void Mesh::constructMesh(const std::vector<MeshElement>& elements) {}

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

void Mesh::changeMaterail(const std::size_t& id, const double& new_value,
  const std::string& target_parameter)
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
