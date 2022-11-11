#include "rbm/mesh.hpp"
#include "xtensor/xslice.hpp"
#include "xtensor/xview.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>

// 1D one-speed example with uniform slab reactor with reflective boundary
// condition at x=0 and vacuum boundary condition at x=a/2. The parameters are
// as follows: a/2 = 30 cm D = 2.0 cm nu_Sigma_f = 0.11 cm^-1
//
// The parameter of interest for this problem is the macroscopic cross section
// of absorption (Sigma_a) and the manually constructed training space consists
// of the following: 0.02, 0.04, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18 cm^-1
//
// The user can choose the target Sigma_a value and the number of nodes in the
// finite differences. The code will produce the resultant RBM approximated
// neutron multiplication factor (k) as well as the one calcualted from 1D
// one-speed finite difference. The relative error will be computed and printed
// as well.

int main(int argc, char* argv[])
{
  // Check that two arguments are given in the command line
  if (argc < 3) {
    throw std::runtime_error("Usage: " + std::string(argv[0]) +
                             " [number-of-nodes] [target-parameter-value]\n");
  }

  // Number of nodes in the mesh
  int N = std::stoi(argv[1]);

  // Target macroscopic xs for absorption in centemeters
  double Sigma_a_t = std::stod(argv[2]);

  // Training points
  // std::vector<double> Sigma_a_train = {0.02, 0.04, 0.08, 0.10, 0.12, 0.14,
  // 0.16, 0.18};
  std::vector<double> Sigma_a_train = {0.08, 0.10, 0.12, 0.14};
  int num_training_points = Sigma_a_train.size();

  xt::xarray<double> flux_train = xt::zeros<double>({num_training_points, N});
  std::vector<double> k_train(Sigma_a_train.size(), 0.0);

  printf("For the following target parameter: Absorption XS\n");
  printf("The training points are: ");
  for (size_t i = 0; i < Sigma_a_train.size(); i++) {
    printf("%lg ", Sigma_a_train[i]);
  }
  printf("\n\n");

  // Constants
  double half_thickness = 30.0; // cm
  double D = 2.0;               // cm
  double nu_Sigma_f = 0.11;     // cm^-1

  // Boundaries
  std::pair<double, double> left_bound =
    std::make_pair(0.0, 1.0); // a/b left boundary condition
  std::pair<double, double> right_bound =
    std::make_pair(1.0, -2.0); // a/b right boundary condition

  printf("Calculating training eigenvectors/values\n");

  // Iterate over each training point
  for (size_t i = 0; i < Sigma_a_train.size(); i++) {
    // Create finite difference mesh
    // Making Mesh object and getting M, F matricies
    Mesh my_mesh(N, half_thickness, left_bound, right_bound);
    my_mesh.run(nu_Sigma_f, Sigma_a_train[i], D);

    // Initializing M, F into variable (could be done in powerIteration)
    xt::xarray<double> M = my_mesh.getM();
    xt::xarray<double> F = my_mesh.getF();

    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto eigenfunction = xt::linalg::eig(A);

    k_train[i] = std::get<0>(eigenfunction)(0, 0).real();

    auto flux = xt::view(flux_train, i, xt::all());
    flux = xt::real(xt::col(std::get<1>(eigenfunction), 0)) * -1;

    // Run power iteration to solve eigenvalue problem
    printf("k = %5lg\n", k_train[i]);
  }

  printf("\nCalculating target value %lg\n", Sigma_a_t);

  // Initialize target matricies
  xt::xarray<double> M_t =
    xt::zeros<double>({num_training_points, num_training_points});
  xt::xarray<double> F_t =
    xt::zeros<double>({num_training_points, num_training_points});

  // Making Mesh object and getting M, F matricies for target parameter
  Mesh my_mesh(N, half_thickness, left_bound, right_bound);
  my_mesh.run(nu_Sigma_f, Sigma_a_t, D);
  xt::xarray<double> M = my_mesh.getM();
  xt::xarray<double> F = my_mesh.getF();

  // Run finite difference to get "exact" solution
  xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
  auto eigenfunction = xt::linalg::eig(A);
  double k_t_exact = std::get<0>(eigenfunction)(0, 0).real();

  // Build target matricies
  for (size_t i = 0; i < Sigma_a_train.size(); i++) {
    for (size_t j = 0; j < Sigma_a_train.size(); j++) {
      // Set the dot of previous result and forward flux
      auto flux_i = xt::view(flux_train, i, xt::all());
      auto flux_j = xt::view(flux_train, i, xt::all());
      M_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(M, flux_j))(0, 0);
      F_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(F, flux_j))(0, 0);
    }
  }

  xt::xarray<double> A_t = xt::linalg::dot(xt::linalg::inv(M_t), F_t);
  auto eigenfunction_t = xt::linalg::eig(A_t);

  // Run power iteration to solve generalized eigenvalue problem
  double k_t = std::get<0>(eigenfunction_t)(0, 0).real();

  // Calculate relative error
  double rel_err = std::abs(k_t - k_t_exact) / k_t_exact;

  // Print results
  printf("For target Sigma_a = %5lg, k_t = %6lg, with %6lg relative error.",
    Sigma_a_t, k_t, rel_err);

  return 0;
}
