#include "../unit_test_framework.hpp"
#include "rbm/mesh.hpp"
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

TEST(regression_test_1D_mesh)
{
  using namespace material;

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

  // Number of course and fine mesh elements
  size_t xN_fine = 5;
  size_t yN_fine = 1;
  size_t xN_course = 1;
  size_t yN_course = 1;

  // Boundary conditions
  std::pair<double, double> reflective_bc = std::make_pair(0, 1);
  std::pair<double, double> vacuum_bc = std::make_pair(1, -2);

  // Initialize mesh
  mesh::Mesh mesh(xN_fine, yN_fine, xN_course, yN_course, reflective_bc,
    vacuum_bc, reflective_bc, reflective_bc);

  // Construct fine grid
  mesh.constructMesh(course_grid);

  // Construct M and F operators
  auto M = mesh.constructM();
  auto F = mesh.constructF();

  // Calculate eigenvalues/vectors
  auto A = xt::linalg::dot(xt::linalg::inv(M), F);
  auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

  // Get fundamental
  double k = eigenvalues(0).real();
  auto fluxes = xt::abs(xt::real(xt::col(eigenvectors, 0)));

  // Expected fundamental
  double k_expected = 1.05519149103;
  xt::xarray<double> fluxes_expected = {
    0.58889, 0.54388, 0.45729, 0.33576, 0.18855};

  // Assertions
  ASSERT_ALMOST_EQUAL(k, k_expected, 0.00000000001);
  for (size_t i = 0; i < fluxes.shape(0); i++) {
    ASSERT_ALMOST_EQUAL(fluxes(i), fluxes_expected(i), 0.00001);
  }
}

TEST_MAIN();
