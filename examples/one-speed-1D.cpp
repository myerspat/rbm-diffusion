#include "../src/linalg.hpp"
#include "../src/matrix.hpp"
#include "../src/mesh.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

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

int main(int argc, char *argv[]) {
  // Check that two arguments are given in the command line
  if (argc < 3) {
    throw std::runtime_error("Usage: " + std::string(argv[0]) + " [number-of-nodes] [target-parameter-value]\n");
  }

  // Number of nodes in the mesh
  int N = std::stoi(argv[1]);

  // Target macroscopic xs for absorption in centemeters
  double Sigma_a_t = std::stod(argv[2]);

  // Training points
  std::vector<double> Sigma_a_train = {0.02, 0.04, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18};
  std::vector<std::vector<double> > flux_train;
  std::vector<double> k_train;

  flux_train.reserve(Sigma_a_train.size());
  k_train.reserve(Sigma_a_train.size());

  // Constants
  double half_thickness = 30.0; // cm
  double D = 2.0;               // cm
  double nu_Sigma_f = 0.11;     // cm^-1

  // Boundaries
  std::pair<double, double> left_bound;  // a/b left boundary condition
  std::pair<double, double> right_bound; // a/b right boundary condition
  left_bound.first = 0.0;                // a
  left_bound.second = 1.0;               // b
  right_bound.first = 1.0;               // a
  right_bound.second = -2.0;              // b

  // Iterate over each training point
  for (double &Sigma_a : Sigma_a_train) {
    // Create finite difference mesh
    // Making Mesh object and getting M, F matricies
    Mesh my_mesh(N, half_thickness, left_bound, right_bound);
    my_mesh.run(nu_Sigma_f, Sigma_a, D);

    // Initializing M, F into variable (could be done in powerIteration)
    linalg::Matrix M = my_mesh.getF();
    linalg::Matrix F = my_mesh.getM();

    // Run power iteration to solve eigenvalue problem
    auto eigen_result = linalg::powerIteration(M, F);
    flux_train.push_back(eigen_result.first);
    k_train.push_back(eigen_result.second);
  }

  // Initialize target matricies
  linalg::Matrix M_t(Sigma_a_train.size(), Sigma_a_train.size());
  linalg::Matrix F_t(Sigma_a_train.size(), Sigma_a_train.size());

  // Build target matricies
  for (size_t i = 0; i < Sigma_a_train.size(); i++) {
    for (size_t j = 0; j < Sigma_a_train.size(); j++) {
      // Intermediate vector-matrix multiplication for adjoint flux
      std::vector<double> M_t_inter = M_t.vecmat(flux_train[i]);
      std::vector<double> F_t_inter = F_t.vecmat(flux_train[i]);

      // Set the dot of previous result and forward flux
      M_t.set(i, j, linalg::dot(M_t_inter, flux_train[j]));
      F_t.set(i, j, linalg::dot(F_t_inter, flux_train[j]));
    }
  }

  // Run power iteration to solve generalized eigenvalue problem
  auto target = linalg::powerIteration(M_t, F_t);
  double k_t = target.second;

  // Run finite difference to get "exact" solution
  // USE MESH CONSTROCTOR HERE

  // Making Mesh object and getting M, F matricies for target parameter
  Mesh my_mesh(N, half_thickness, left_bound, right_bound);
  my_mesh.run(nu_Sigma_f, Sigma_a_t, D);
  linalg::Matrix M = my_mesh.getF();
  linalg::Matrix F = my_mesh.getM();

  auto exact = linalg::powerIteration(M_t, F_t);
  double k_t_exact = exact.second;

  // Calculate relative error
  double rel_err = std::abs(k_t - k_t_exact) / k_t_exact;

  // Print results
  printf("For target Sigma_a = %5lg, k_t = %6lg, with %6lg relative error.", Sigma_a_t, k_t, rel_err);

  return 0;
}
