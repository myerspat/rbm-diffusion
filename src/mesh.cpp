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

  return course_grid;
}

xt::xarray<MeshElement> Mesh::constructFineGrid(
  const xt::xarray<MeshElement>& course_grid)
{
  // Check the adjacent lengths of MeshElements are the same
  assert(checkSharedLengths(course_grid));

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
  bool shared_lengths = true;

  return shared_lengths;
}

void Mesh::changeMaterail(const std::size_t& id, const double& new_value,
  const rbm::Parameter& target_parameter)
{}

xt::xarray<double> Mesh::constructF()
{
  // Allocate space for fission operator
  xt::xarray<double> F = xt::zeros<double>({getSize(), getSize()});

  // Fill diagonal array
  for (size_t i = 0; i < _fine_grid.shape(0); i++) {
    for (size_t j = 0; j < _fine_grid.shape(1); j++) {
      // Convert 2D idx to 1D idx assuming row major ordering
      size_t position = ravelIDX(i, j);

      // Fill diagonal
      F(position, position) = _fine_grid(i, j).getMaterial().getNuFission() *
                              _fine_grid(i, j).getLX() / _xN_fine *
                              _fine_grid(i, j).getLY() / _yN_fine;
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
  xt::xarray<double> M = xt::zeros<double>({getSize(), getSize()});

  // Lambda functions for calculating dx and dy of MeshElement i, j
  auto dx = [&](size_t i, size_t j) {
    return _fine_grid(i, j).getLX() / _xN_fine;
  };

  auto dy = [&](size_t i, size_t j) {
    return _fine_grid(i, j).getLY() / _yN_fine;
  };

  // Lambda function for boundary coupling coefficient for interface between i,
  // j and i, j + 1
  auto a_x = [&](size_t i, size_t j) {
    return dy(i, j) /
           (dx(i, j) / (2 * _fine_grid(i, j).getMaterial().getD()) +
             dx(i, j + 1) / (2 * _fine_grid(i, j + 1).getMaterial().getD()));
  };

  // Lambda function for boundary coupling coefficient for interface between i,
  // j and i + 1, j
  auto a_y = [&](size_t i, size_t j) {
    return dx(i, j) /
           (dy(i, j) / (2 * _fine_grid(i, j).getMaterial().getD()) +
             dy(i + 1, j) / (2 * _fine_grid(i + 1, j).getMaterial().getD()));
  };

  // Lambda function for boundary coupling coefficient for interface between i,
  // j and boundary along y
  auto a_xb = [&](size_t i, size_t j, std::pair<double, double>& bound) {
    const double& a = bound.first;
    const double& b = bound.second;

    return a == 0 ? 0
                  : dy(i, j) /
                      (dx(i, j) / (2 * _fine_grid(i, j).getMaterial().getD()) -
                        b / a);
  };

  // Lambda function for boundary coupling coefficient for interface between i,
  // j and boundary along x
  auto a_yb = [&](size_t i, size_t j, std::pair<double, double>& bound) {
    const double& a = bound.first;
    const double& b = bound.second;

    return a == 0 ? 0
                  : dx(i, j) /
                      (dy(i, j) / (2 * _fine_grid(i, j).getMaterial().getD()) -
                        b / a);
  };

  // Last indicies along x and y of _fine_grid
  size_t x_last = getXN() - 1;
  size_t y_last = getYN() - 1;

  // If y_last == 0 then the problem is 1D else 2D
  if (y_last == 0) {
    // Left boundary element
    M(ravelIDX(0, 0), ravelIDX(0, 1)) = -a_x(0, 0);
    M(ravelIDX(0, 0), ravelIDX(0, 0)) =
      _fine_grid(0, 0).getMaterial().getAbsorption() * dx(0, 0) -
      M(ravelIDX(0, 0), ravelIDX(0, 1)) + a_xb(0, 0, _left_bound);

    // Length (no edges)
    for (size_t j = 1; j < x_last; j++) {
      M(ravelIDX(0, j), ravelIDX(0, j + 1)) = -a_x(0, j);
      M(ravelIDX(0, j), ravelIDX(0, j - 1)) = -a_x(0, j - 1);
      M(ravelIDX(0, j), ravelIDX(0, j)) =
        _fine_grid(0, j).getMaterial().getAbsorption() * dx(0, j) * dy(0, j) -
        M(ravelIDX(0, j), ravelIDX(0, j + 1)) -
        M(ravelIDX(0, j), ravelIDX(0, j - 1));
    }

    // Right boundary element
    M(ravelIDX(0, x_last), ravelIDX(0, x_last - 1)) = -a_x(0, x_last - 1);
    M(ravelIDX(0, x_last), ravelIDX(0, x_last)) =
      _fine_grid(0, x_last).getMaterial().getAbsorption() * dx(0, x_last) *
        dy(0, x_last) -
      M(ravelIDX(0, x_last), ravelIDX(0, x_last - 1)) +
      a_xb(0, x_last, _right_bound);

  } else {
    // Bottom left corner
    M(ravelIDX(0, 0), ravelIDX(0, 1)) = -a_x(0, 0);
    M(ravelIDX(0, 0), ravelIDX(1, 0)) = -a_y(0, 0);
    M(ravelIDX(0, 0), ravelIDX(0, 0)) =
      _fine_grid(0, 0).getMaterial().getAbsorption() * dx(0, 0) * dy(0, 0) -
      M(ravelIDX(0, 0), ravelIDX(0, 1)) - M(ravelIDX(0, 0), ravelIDX(1, 0)) +
      a_xb(0, 0, _left_bound) + a_yb(0, 0, _bottom_bound);

    // Bottom edge (no corners)
    for (size_t j = 1; j < x_last; j++) {
      M(ravelIDX(0, j), ravelIDX(0, j + 1)) = -a_x(0, j);
      M(ravelIDX(0, j), ravelIDX(1, j)) = -a_y(0, j);
      M(ravelIDX(0, j), ravelIDX(0, j - 1)) = -a_x(0, j - 1);
      M(ravelIDX(0, j), ravelIDX(0, j)) =
        _fine_grid(0, j).getMaterial().getAbsorption() * dx(0, j) * dy(0, j) -
        M(ravelIDX(0, j), ravelIDX(0, j + 1)) -
        M(ravelIDX(0, j), ravelIDX(1, j)) -
        M(ravelIDX(0, j), ravelIDX(0, j - 1)) + a_yb(0, j, _bottom_bound);
    }

    // Bottom right corner
    M(ravelIDX(0, x_last), ravelIDX(0, x_last - 1)) = -a_x(0, x_last - 1);
    M(ravelIDX(0, x_last), ravelIDX(1, x_last)) = -a_y(0, x_last);
    M(ravelIDX(0, x_last), ravelIDX(0, x_last)) =
      _fine_grid(0, x_last).getMaterial().getAbsorption() * dx(0, x_last) *
        dy(0, x_last) -
      M(ravelIDX(0, x_last), ravelIDX(0, x_last - 1)) -
      M(ravelIDX(0, x_last), ravelIDX(1, x_last)) +
      a_xb(0, x_last, _right_bound) + a_yb(0, x_last, _bottom_bound);

    // Left and right edge (no corners)
    for (size_t i = 1; i < y_last; i++) {
      M(ravelIDX(i, 0), ravelIDX(i, 1)) = -a_x(i, 0);
      M(ravelIDX(i, 0), ravelIDX(i + 1, 0)) = -a_y(i, 0);
      M(ravelIDX(i, 0), ravelIDX(i - 1, 0)) = -a_y(i - 1, 0);
      M(ravelIDX(i, 0), ravelIDX(i, 0)) =
        _fine_grid(i, 0).getMaterial().getAbsorption() * dx(i, 0) * dy(i, 0) -
        M(ravelIDX(i, 0), ravelIDX(i, 1)) -
        M(ravelIDX(i, 0), ravelIDX(i + 1, 0)) -
        M(ravelIDX(i, 0), ravelIDX(i - 1, 0)) + a_xb(i, 0, _left_bound);

      M(ravelIDX(i, x_last), ravelIDX(i, x_last - 1)) = -a_x(i, x_last - 1);
      M(ravelIDX(i, x_last), ravelIDX(i + 1, x_last)) = -a_y(i, x_last);
      M(ravelIDX(i, x_last), ravelIDX(i - 1, x_last)) = -a_y(i - 1, x_last);
      M(ravelIDX(i, x_last), ravelIDX(i, x_last)) =
        _fine_grid(i, x_last).getMaterial().getAbsorption() * dx(i, x_last) *
          dy(i, x_last) -
        M(ravelIDX(i, x_last), ravelIDX(i, x_last - 1)) -
        M(ravelIDX(i, x_last), ravelIDX(i + 1, x_last)) -
        M(ravelIDX(i, x_last), ravelIDX(i - 1, x_last)) +
        a_xb(i, x_last, _right_bound);
    }

    // Top left corner
    M(ravelIDX(y_last, 0), ravelIDX(y_last, 1)) = -a_x(y_last, 0);
    M(ravelIDX(y_last, 0), ravelIDX(y_last - 1, 0)) = -a_y(y_last - 1, 0);
    M(ravelIDX(y_last, 0), ravelIDX(y_last, 0)) =
      _fine_grid(y_last, 0).getMaterial().getAbsorption() * dx(y_last, 0) *
        dy(y_last, 0) -
      M(ravelIDX(y_last, 0), ravelIDX(y_last, 1)) -
      M(ravelIDX(y_last, 0), ravelIDX(y_last - 1, 0)) +
      a_xb(y_last, 0, _left_bound) + a_yb(y_last, 0, _top_bound);

    // Top edge (no corners)
    for (size_t j = 1; j < x_last; j++) {
      M(ravelIDX(y_last, j), ravelIDX(y_last, j + 1)) = -a_x(y_last, j);
      M(ravelIDX(y_last, j), ravelIDX(y_last - 1, j)) = -a_y(y_last - 1, j);
      M(ravelIDX(y_last, j), ravelIDX(y_last, j - 1)) = -a_x(y_last, j - 1);
      M(ravelIDX(y_last, j), ravelIDX(y_last, j)) =
        _fine_grid(y_last, j).getMaterial().getAbsorption() * dx(y_last, j) *
          dy(y_last, j) -
        M(ravelIDX(y_last, j), ravelIDX(y_last, j + 1)) -
        M(ravelIDX(y_last, j), ravelIDX(y_last - 1, j)) -
        M(ravelIDX(y_last, j), ravelIDX(y_last, j - 1)) +
        a_yb(y_last, 0, _top_bound);
    }

    // Top right corner
    M(ravelIDX(y_last, x_last), ravelIDX(y_last, x_last - 1)) =
      -a_x(y_last, x_last - 1);
    M(ravelIDX(y_last, x_last), ravelIDX(y_last - 1, x_last)) =
      -a_y(y_last - 1, x_last);
    M(ravelIDX(y_last, x_last), ravelIDX(y_last, x_last)) =
      _fine_grid(y_last, x_last).getMaterial().getAbsorption() *
        dx(y_last, x_last) * dy(y_last, x_last) -
      M(ravelIDX(y_last, x_last), ravelIDX(y_last, x_last - 1)) -
      M(ravelIDX(y_last, x_last), ravelIDX(y_last - 1, x_last)) +
      a_xb(y_last, x_last, _right_bound) + a_yb(y_last, x_last, _top_bound);

    // Inner elements
    for (size_t i = 1; i < getYN() - 1; i++) {
      for (size_t j = 1; j < getXN() - 1; j++) {
        M(ravelIDX(i, j), ravelIDX(i, j - 1)) = -a_x(i, j - 1);
        M(ravelIDX(i, j), ravelIDX(i - 1, j)) = -a_y(i - 1, j);
        M(ravelIDX(i, j), ravelIDX(i, j + 1)) = -a_x(i, j);
        M(ravelIDX(i, j), ravelIDX(i + 1, j)) = -a_y(i, j);
        M(ravelIDX(i, j), ravelIDX(i, j)) =
          _fine_grid(i, j).getMaterial().getAbsorption() * dx(i, j) * dy(i, j) -
          M(ravelIDX(i, j), ravelIDX(i, j - 1)) -
          M(ravelIDX(i, j), ravelIDX(i - 1, j)) -
          M(ravelIDX(i, j), ravelIDX(i, j + 1)) -
          M(ravelIDX(i, j), ravelIDX(i + 1, j));
      }
    }
  }

  return M;
}

} // namespace mesh
