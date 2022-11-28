#include "rbm/mesh.hpp"
#include "unit_test_framework.hpp"
#include <utility>
#include <xtensor/xarray.hpp>

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

TEST_MAIN();
