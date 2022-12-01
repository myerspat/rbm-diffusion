#include "../unit_test_framework.hpp"
#include "rbm/mesh.hpp"
#include "xtensor/xutils.hpp"
#include <cmath>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnorm.hpp>

TEST(regression_test_2D_mesh)
{
  // Initialize materials
  Material fuel("fuel", 0.10, 0.1109662271, 2.0);

  // Lengths
  double lx_0 = 30.0;
  double ly_0 = 30.0;

  // Create mesh elements
  // Fuel region
  size_t id_0 = 0;
  mesh::MeshElement e_1(
    fuel, lx_0, ly_0, id_0, std::make_pair(0, 0), std::make_pair(0, 0));

  // Create course grid
  std::vector<mesh::MeshElement> course_grid = {e_1};

  // Number of course and fine mesh elements
  size_t xN_fine = 20;
  size_t yN_fine = 20;
  size_t xN_course = 1;
  size_t yN_course = 1;

  // Boundary conditions
  std::pair<double, double> reflective_bc = std::make_pair(0, 1);
  std::pair<double, double> zero_flux_bc = std::make_pair(1, 0);

  // Initialize mesh
  mesh::Mesh mesh(xN_fine, yN_fine, xN_course, yN_course, reflective_bc,
    zero_flux_bc, zero_flux_bc, reflective_bc);

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

  // Initialize analytical solution
  xt::xarray<double> analytical_flux = xt::zeros<double>({mesh.getSize()});

  // Calculate the bare slab reactor analytical solution for BCs
  for (size_t i = 0; i < mesh.getYN(); i++) {
    for (size_t j = 0; j < mesh.getXN(); j++) {
      double pi = xt::numeric_constants<double>::PI;
      double x = (j + 1.0 / 2.0) * lx_0 / xN_fine;
      double y = (i + 1.0 / 2.0) * ly_0 / yN_fine;
      analytical_flux(mesh.ravelIDX(i, j)) =
        std::cos(pi * x / (lx_0 * 2)) * std::cos(pi * y / (ly_0 * 2));
    }
  }

  // Normalize analytical solution to a unit vector
  analytical_flux /= xt::norm_l2(analytical_flux);

  // Check that the eigenvalue is critical
  ASSERT_ALMOST_EQUAL(k, 1.0, 0.0001);

  // Check each flux solution
  for (size_t i = 0; i < fluxes.shape(0); i++) {
    ASSERT_ALMOST_EQUAL(fluxes(i), analytical_flux(i), 0.000000001);
  }
}

TEST_MAIN();
