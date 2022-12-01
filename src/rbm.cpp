#include "rbm/rbm.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <cstdio>
#include <tuple>
#include <utility>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xview.hpp>

namespace rbm {

void Perturb::pcaReduce(xt::xarray<double>& training_fluxes)
{
  // Center flux data
  xt::xarray<double> col_means = xt::mean(training_fluxes, 0);
  for (size_t i = 0; i < training_fluxes.shape(0); i++) {
    auto center = xt::view(training_fluxes, i, xt::all());
    center -= col_means;
  }

  // Rank of flux data
  size_t rank = xt::linalg::matrix_rank(training_fluxes);

  // Singular value decomposition
  auto [U, L, At] = xt::linalg::svd(training_fluxes);

  // Extract U, L, and A' (At) from X = ULA'
  // l is given as a row vector but is really a diagonal metrix that is
  // rank x rank
  U = xt::view(U, xt::all(), xt::range(0, rank));

  // Calculate total variance
  double total_variance = xt::sum(xt::square(L))(0);

  // Reduce subspace to the first PCs to _num_pcs
  double var = 0.0;
  for (size_t i = 0; i < _num_pcs; i++) {
    var += (L(i) * L(i));
  }

  // Throw warning if the number of PCs preserved is less than 90% of the total
  // variance
  if (var < 0.9 * total_variance) {
    printf("Warning: Subspace was reduced to %lu PCs which has only %lg "
           "percent of the total variance\n",
      _num_pcs, var / total_variance * 100);
  }

  // Reduced U, L, and At
  xt::xarray<double> U_r = xt::view(U, xt::all(), xt::range(0, _num_pcs));
  xt::xarray<double> L_r = xt::view(L, xt::range(0, _num_pcs));
  xt::xarray<double> At_r =
    xt::view(At, xt::range(0, _num_pcs), xt::range(0, _num_pcs));

  // Calculate reduced training_fluxes
  training_fluxes = xt::linalg::dot(U_r, xt::linalg::dot(xt::diag(L_r), At_r));

  // Uncenter training_fluxes with col_means from 0 to num_pcs
  col_means = xt::view(col_means, xt::range(0, _num_pcs));
  for (size_t i = 0; i < training_fluxes.shape(0); i++) {
    auto center = xt::view(training_fluxes, i, xt::all());
    center += col_means;
  }
}

xt::xarray<double> Perturb::constructF_t(
  const xt::xarray<double>& F, const xt::xarray<double>& training_fluxes)
{
  // Initialize F_t and allocate space
  xt::xarray<double> F_t =
    xt::zeros<double>({training_fluxes.shape(1), training_fluxes.shape(1)});

  for (size_t i = 0; i < training_fluxes.shape(1); i++) {
    // Get flux at i
    xt::xarray<double> flux_i = xt::view(training_fluxes, xt::all(), i);
    for (size_t j = 0; j < training_fluxes.shape(1); j++) {
      // Get flux at j
      xt::xarray<double> flux_j = xt::view(training_fluxes, xt::all(), j);
      // Calculate F_t
      F_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(F, flux_j))(0);
    }
  }

  return F_t;
}

xt::xarray<double> Perturb::constructM_t(
  const xt::xarray<double>& M, const xt::xarray<double>& training_fluxes)
{
  // Initialize M_t and allocate space
  xt::xarray<double> M_t =
    xt::zeros<double>({training_fluxes.shape(1), training_fluxes.shape(1)});

  for (size_t i = 0; i < training_fluxes.shape(1); i++) {
    // Get flux at i
    xt::xarray<double> flux_i = xt::view(training_fluxes, xt::all(), i);
    for (size_t j = 0; j < training_fluxes.shape(1); j++) {
      // Get flux at j
      xt::xarray<double> flux_j = xt::view(training_fluxes, xt::all(), j);
      // Calculate F_t
      M_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(M, flux_j))(0);
    }
  }

  return M_t;
}

void Perturb::initialize(
  xt::xarray<double>& training_points, mesh::Mesh& mesh, size_t& element_id)
{}

std::pair<double, xt::xarray<double>> Perturb::findMaxEigen(
  xt::xarray<double> eigenvalues, xt::xarray<double> eigenvectors)
{
  size_t idx = 0;
  double max_eigenvalue = 0.0;
  for (size_t i = 0; i < eigenvalues.size(); i++) {
    if (eigenvalues(i) > max_eigenvalue) {
      max_eigenvalue = eigenvalues(i);
      idx = i;
    }
  }

  return std::make_pair(eigenvalues(idx), xt::abs(xt::col(eigenvectors, idx)));
}

void Perturb::train()
{
  // Prompt user of training begun
  printf("\n Begin training\n");

  // full training_fluxes and training_k
  _training_fluxes = xt::xarray<double>::from_shape(
    {_mesh.getSize(), _training_points.shape(0)});
  _training_k = xt::xarray<double>::from_shape({_training_points.shape(0)});

  for (size_t i = 0; i < _training_points.shape(0); i++) {
    // set parameter
    _mesh.changeMaterial(_element_id, _training_points(i), _target_parameter);

    // find F and M matrices
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // find eigenvalues and eigenvectors
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

    // Find max eigenvalue corresponding eigenvector
    auto fundamental =
      findMaxEigen(xt::real(eigenvalues), xt::real(eigenvectors));

    _training_k(i) = fundamental.first;
    xt::col(_training_fluxes, i) = fundamental.second;

    // Update user
    printf("   Point %lu = %3lg => k = %6lg\n", i + 1, _training_points(i),
      _training_k(i));
  }

  // reduce to PxP
  pcaReduce(_training_fluxes);
}

void Perturb::calcTargets()
{
  // Prompt user of online mode
  printf("\n Calculating Targets\n");

  _target_fluxes =
    xt::zeros<double>({_training_fluxes.shape(0), _target_points.size()});
  _target_k = xt::xarray<double>::from_shape({_target_points.size()});

  for (size_t i = 0; i < _target_points.size(); i++) {
    // change with specific parameter
    _mesh.changeMaterial(_element_id, _target_points(i), _target_parameter);

    // get F and M matricies
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // get F_t and M_t
    xt::xarray<double> F_t = constructF_t(F, _training_fluxes);
    xt::xarray<double> M_t = constructM_t(M, _training_fluxes);

    // calculate the eigenvlue and eigenvector for target
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M_t), F_t);
    auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

    // Find max eigenvalue corresponding eigenvector
    auto fundamental =
      findMaxEigen(xt::real(eigenvalues), xt::real(eigenvectors));

    // Get target eigenvalue
    _target_k(i) = fundamental.first;

    // Calculate approximate target flux using the sum of ritz vector *
    // _training_fluxes
    for (size_t j = 0; j < _training_points.size(); j++) {
      xt::col(_target_fluxes, i) +=
        fundamental.second(j) * xt::col(_training_fluxes, j);
    }

    // Normlize target fluxes
    xt::col(_target_fluxes, i) /= xt::norm_l2(xt::col(_target_fluxes, i));

    // Update user
    printf("   Point %lu = %3lg => k = %6lg\n", i + 1, _target_points(i),
      _target_k(i));
  }
}

void Perturb::checkError()
{
  printf("\n Calculating Errors\n");

  auto exact_target_fluxes =
    xt::xarray<double>::from_shape({_mesh.getSize(), _target_points.size()});
  auto exact_target_k = xt::xarray<double>::from_shape({_target_points.size()});
  auto error = xt::xarray<double>::from_shape({_target_points.size()});

  for (size_t i = 0; i < _target_points.shape(0); i++) {
    // set parameter
    _mesh.changeMaterial(_element_id, _target_points(i), _target_parameter);

    // find F and M matrices
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // find eigenvalues and eigenvectors
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

    // Find max eigenvalue corresponding eigenvector
    auto fundamental =
      findMaxEigen(xt::real(eigenvalues), xt::real(eigenvectors));

    exact_target_k(i) = fundamental.first;
    xt::col(exact_target_fluxes, i) = fundamental.second;

    // Calculate relative error
    error(i) = (exact_target_k(i) - _target_k(i)) / exact_target_k(i);

    // Update user
    printf(
      "   Point %lu = %3lg => k_exact = %6lg, k_rbm = %6lg, error = %6lg\n",
      i + 1, _target_points(i), exact_target_k(i), _target_k(i), error(i));
  }
}

} // namespace rbm
