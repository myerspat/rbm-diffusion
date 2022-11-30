#include "rbm/material.hpp"
#include "rbm/mesh.hpp"
#include "rbm/meshElement.hpp"
#include "unit_test_framework.hpp"
#include <cstddef>
#include <utility>
#include <xtensor/xarray.hpp>

TEST(test_constructM_1)
{
  // This problem tests the 1D implementation of constructM

  // Initialize materials
  Material fuel("fuel", 0.10, 0.11, 2.0);

  // Lengths
  double lx_0 = 30.0;
  double ly_0 = 1.0;

  // Create mesh elements
  // Fuel region
  size_t id_0 = 0;
  mesh::MeshElement e_1(
    fuel, lx_0, ly_0, id_0, std::make_pair(0, 0), std::make_pair(0, 0));

  // Create course grid
  std::vector<mesh::MeshElement> course_grid = {e_1};

  // Initialize mesh
  size_t xN_fine = 5;
  size_t yN_fine = 1;
  size_t xN_course = 1;
  size_t yN_course = 1;
  std::pair<double, double> bounds = std::make_pair(0, 1);
  std::pair<double, double> left_bound = std::make_pair(0, 1);
  std::pair<double, double> right_bound = std::make_pair(1, -2);

  mesh::Mesh mesh(xN_fine, yN_fine, xN_course, yN_course, left_bound,
    right_bound, bounds, bounds);

  // Construct fine grid
  mesh.constructMesh(course_grid);
  auto M = mesh.constructM();

  // Expected result
  xt::xarray<double> M_expected = {{0.9333, -0.3333, 0, 0, 0},
    {-0.3333, 1.2666, -0.3333, 0, 0}, {0, -0.3333, 1.2666, -0.3333, 0},
    {0, 0, -0.3333, 1.2666, -0.3333}, {0, 0, 0, -0.3333, 1.2190}};

  // Check each value
  for (size_t i = 0; i < M_expected.shape(0); i++) {
    for (size_t j = 0; j < M_expected.shape(1); j++) {
      ASSERT_ALMOST_EQUAL(M(i, j), M_expected(i, j), 0.0001);
    }
  }
}

TEST(test_constructFineGrid_1)
{
  // Initialize materials
  Material fuel("fuel", 0.10, 0.11, 2.0);
  Material reflector("reflector", 0.01, 0.0, 1.5);

  // Lengths
  double lx_0 = 30.0;
  double ly_0 = 30.0;
  double lx_1 = 10.0;
  double ly_1 = 10.0;

  // Create mesh elements
  // Fuel region
  size_t id_0 = 0;
  mesh::MeshElement e_1(
    fuel, lx_0, ly_0, id_0, std::make_pair(0, 0), std::make_pair(0, 0));

  // Reflector region
  size_t id_1 = 1;
  mesh::MeshElement e_2(
    reflector, lx_1, ly_0, id_1, std::make_pair(0, 0), std::make_pair(1, 1));
  mesh::MeshElement e_3(
    reflector, lx_0, ly_1, id_1, std::make_pair(1, 1), std::make_pair(0, 0));
  mesh::MeshElement e_4(
    reflector, lx_1, ly_1, id_1, std::make_pair(1, 1), std::make_pair(1, 1));

  // Create course grid
  xt::xarray<mesh::MeshElement> course_grid = {{e_1, e_2}, {e_3, e_4}};

  // Initialize mesh
  size_t xN_fine = 10;
  size_t yN_fine = 10;
  size_t xN_course = 2;
  size_t yN_course = 2;
  std::pair<double, double> bounds = std::make_pair(0.0, 1.0);
  mesh::Mesh mesh(
    xN_fine, yN_fine, xN_course, yN_course, bounds, bounds, bounds, bounds);

  // Construct fine grid
  auto fine_grid = mesh.constructFineGrid(course_grid);

  // Assertions
  ASSERT_EQUAL(mesh.getSize(), 400);
  for (size_t i = 0; i < 10; i++) {
    // Fuel region
    for (size_t j = 0; j < 10; j++) {
      ASSERT_EQUAL(
        fine_grid(i, j).getMaterial().getName(), e_1.getMaterial().getName());
    }

    // Right reflector region
    for (size_t j = 10; j < 20; j++) {
      ASSERT_EQUAL(
        fine_grid(i, j).getMaterial().getName(), e_2.getMaterial().getName());
      ASSERT_EQUAL(fine_grid(i, j).getLX(), e_2.getLX());
      ASSERT_EQUAL(fine_grid(i, j).getLY(), e_2.getLY());
    }
  }

  // Upper reflector region
  for (size_t i = 10; i < 20; i++) {
    for (size_t j = 0; j < 10; j++) {
      ASSERT_EQUAL(
        fine_grid(i, j).getMaterial().getName(), e_3.getMaterial().getName());
      ASSERT_EQUAL(fine_grid(i, j).getLX(), e_3.getLX());
      ASSERT_EQUAL(fine_grid(i, j).getLY(), e_3.getLY());
    }

    // Right reflector region
    for (size_t j = 10; j < 20; j++) {
      ASSERT_EQUAL(
        fine_grid(i, j).getMaterial().getName(), e_4.getMaterial().getName());
      ASSERT_EQUAL(fine_grid(i, j).getLX(), e_4.getLX());
      ASSERT_EQUAL(fine_grid(i, j).getLY(), e_4.getLY());
    }
  }
}

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

TEST(test_constructF_1)
{
  // Initialize materials
  Material fuel("fuel", 0.10, 0.11, 2.0);
  Material reflector("reflector", 0.01, 0.0, 1.5);

  // Lengths
  double lx_0 = 30.0;
  double ly_0 = 30.0;
  double lx_1 = 10.0;
  double ly_1 = 10.0;

  // Create mesh elements
  // Fuel region
  size_t id_0 = 0;
  mesh::MeshElement e_1(
    fuel, lx_0, ly_0, id_0, std::make_pair(0, 0), std::make_pair(0, 0));

  // Reflector region
  size_t id_1 = 1;
  mesh::MeshElement e_2(
    reflector, lx_1, ly_0, id_1, std::make_pair(0, 0), std::make_pair(1, 1));
  mesh::MeshElement e_3(
    reflector, lx_0, ly_1, id_1, std::make_pair(1, 1), std::make_pair(0, 0));
  mesh::MeshElement e_4(
    reflector, lx_1, ly_1, id_1, std::make_pair(1, 1), std::make_pair(1, 1));

  // Create course grid
  std::vector<mesh::MeshElement> elements = {e_1, e_2, e_3, e_4};

  // Initialize mesh
  size_t xN_fine = 10;
  size_t yN_fine = 10;
  size_t xN_course = 2;
  size_t yN_course = 2;
  std::pair<double, double> bounds = std::make_pair(0.0, 1.0);
  mesh::Mesh mesh(
    xN_fine, yN_fine, xN_course, yN_course, bounds, bounds, bounds, bounds);

  // Construct mesh given vector of elements
  mesh.constructMesh(elements);

  // dx and dy
  double dx_fuel = lx_0 / xN_fine;
  double dy_fuel = ly_0 / yN_fine;
  double dx_reflector = lx_1 / xN_fine;
  double dy_reflector = ly_1 / yN_fine;

  auto F = mesh.constructF();
  for (size_t i = 0; i < 10; i++) {
    ASSERT_ALMOST_EQUAL(
      F(i, i), fuel.getNuFission() * dx_fuel * dy_fuel, 0.0000001);
  }
  for (size_t i = 10; i < 20; i++) {
    ASSERT_ALMOST_EQUAL(F(i, i),
      reflector.getNuFission() * dx_reflector * dy_reflector, 0.0000001);
  }
  for (size_t i = 20; i < 30; i++) {
    ASSERT_ALMOST_EQUAL(
      F(i, i), fuel.getNuFission() * dx_fuel * dy_fuel, 0.0000001);
  }
  for (size_t i = mesh.getSize() - 20; i < mesh.getSize(); i++) {
    ASSERT_ALMOST_EQUAL(F(i, i),
      reflector.getNuFission() * dx_reflector * dy_reflector, 0.0000001);
  }
}

TEST(test_CheckSharedLengths)
{
  // Initialize materials
  Material fuel("fuel", 0.10, 0.11, 2.0);
  Material reflector("reflector", 0.01, 0.0, 1.5);

  // Lengths
  double lx_0 = 30.0;
  double ly_0 = 30.0;
  double lx_1 = 10.0;
  double ly_1 = 10.0;

  // Create mesh elements
  // Fuel region
  size_t id_0 = 0;
  mesh::MeshElement e_1(
    fuel, lx_0, ly_0, id_0, std::make_pair(0, 0), std::make_pair(0, 0));

  // Reflector region
  size_t id_1 = 1;
  mesh::MeshElement e_2(
    reflector, lx_1, ly_0, id_1, std::make_pair(0, 0), std::make_pair(1, 1));
  mesh::MeshElement e_3(
    reflector, lx_0, ly_1, id_1, std::make_pair(1, 1), std::make_pair(0, 0));
  mesh::MeshElement e_4(
    reflector, lx_1, ly_1, id_1, std::make_pair(1, 1), std::make_pair(1, 1));

  // Create course grid
  xt::xarray<mesh::MeshElement> course_grid = {{e_1, e_2}, {e_3, e_4}};

  // Initialize mesh
  size_t xN_fine = 10;
  size_t yN_fine = 10;
  size_t xN_course = 2;
  size_t yN_course = 2;
  std::pair<double, double> bounds = std::make_pair(0.0, 1.0);
  mesh::Mesh mesh(
    xN_fine, yN_fine, xN_course, yN_course, bounds, bounds, bounds, bounds);

  // Image of course_grid -> Which is correct each columns lx are the same and
  // each rows ly are the same
  //------------------------------
  // e3 = (1,1) (0,0) | e4 = (1,1) (1,1)
  // xl = 30          | lx = 10
  // ly = 10          | ly = 10
  // ---------------------------
  // e1 = (0,0) (0,0) | e2 = (0,0) (1,1)
  // xl = 30          | lx = 10
  // ly = 30          | ly = 30
  ASSERT_TRUE(mesh.checkSharedLengths(course_grid) == true);

  // constructing bad course grid where lx's are wrong
  size_t lx_2 = 100;
  mesh::MeshElement e_5(
    reflector, lx_2, ly_1, id_1, std::make_pair(1, 1), std::make_pair(0, 0));
  // Image of course_grid_bad_in_xl -> Which is incorrect in the lx e1 and e5
  // meshElements
  //------------------------------
  // e5 = (1,1) (0,0)      | e4 = (1,1) (1,1)
  // xl = 100              | lx = 10
  // ly = 10               | ly = 10
  // ---------------------------
  // e1 = (0,0) (0,0) | e2 = (0,0) (1,1)
  // xl = 30          | lx = 10
  // ly = 30          | ly = 30
  xt::xarray<mesh::MeshElement> course_grid_bad_in_xl = {
    {e_1, e_2}, {e_5, e_4}};
  ASSERT_TRUE(mesh.checkSharedLengths(course_grid_bad_in_xl) == false);

  // Image of course_grid_bad_in_ly -> Which is incorrect in the ly e6 and e4
  // meshElements
  //------------------------------
  // e5 = (1,1) (0,0)      | e4 = (1,1) (1,1)
  // xl = 30               | lx = 10
  // ly = 100              | ly = 10
  //                      | -----------------------
  //                      |
  //                      |
  // ---------------------------
  // e1 = (0,0) (0,0)     | e2 = (0,0) (1,1)
  // xl = 30              | lx = 10
  // ly = 30              | ly = 30

  lx_2 = 30;
  size_t ly_2 = 100;
  mesh::MeshElement e_6(
    reflector, lx_2, ly_2, id_1, std::make_pair(1, 1), std::make_pair(0, 0));
  xt::xarray<mesh::MeshElement> course_grid_bad_in_ly = {
    {e_1, e_2}, {e_6, e_4}};
  ASSERT_TRUE(mesh.checkSharedLengths(course_grid_bad_in_ly) == false);
}

TEST_MAIN();
