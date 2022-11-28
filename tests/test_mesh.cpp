#include "rbm/material.hpp"
#include "rbm/mesh.hpp"
#include "rbm/meshElement.hpp"
#include "unit_test_framework.hpp"
#include <cstddef>
#include <utility>
#include <xtensor/xarray.hpp>

TEST(test_constructCourseGrid_1)
{
  // Material 1
  Material mat_1("fuel", 0.10, 0.11, 2.0);

  // Element 1
  double lx_1 = 30;
  double ly_1 = 10;
  size_t id_1 = 0;
  std::pair<size_t, size_t> idx_row_1 = std::make_pair(0, 1);
  std::pair<size_t, size_t> idx_col_1 = std::make_pair(0, 0);
  mesh::MeshElement element_1(mat_1, lx_1, ly_1, id_1, idx_row_1, idx_col_1);

  // Material 2
  Material mat_2("reflector", 0.01, 0.0, 1.5);

  // Element 2
  double lx_2 = 10;
  double ly_2 = 10;
  size_t id_2 = 0;
  std::pair<size_t, size_t> idx_row_2 = std::make_pair(0, 1);
  std::pair<size_t, size_t> idx_col_2 = std::make_pair(1, 1);
  mesh::MeshElement element_2(mat_2, lx_2, ly_2, id_2, idx_row_2, idx_col_2);

  // Elements vector
  std::vector<mesh::MeshElement> elements = {element_1, element_2};

  // Initialize mesh
  size_t xN_fine, yN_fine;
  size_t xN_course = 2, yN_course = 2;
  std::pair<double, double> bounds = std::make_pair(0, 1);
  mesh::Mesh mesh(
    xN_fine, yN_fine, xN_course, yN_course, bounds, bounds, bounds, bounds);

  // Create xarray of MeshElements using the vector
  auto course_grid = mesh.constructCourseGrid(elements);

  // Assertions
  ASSERT_EQUAL("reflector", course_grid(1, 1).getMaterial().getName());
  ASSERT_EQUAL(lx_1, course_grid(0, 0).getLX());
  ASSERT_EQUAL(lx_2, course_grid(0, 1).getLX());
}

TEST_MAIN();
