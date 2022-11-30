#include "rbm/rbm.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <cstdio>
#include <tuple>
#include <utility>
#include <xtensor-blas/xlinalg.hpp>
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
  xt::xarray<double> F_t = xt::xarray<double>::from_shape(
    {training_fluxes.shape(1), training_fluxes.shape(1)});

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
  xt::xarray<double> M_t = xt::xarray<double>::from_shape(
    {training_fluxes.shape(1), training_fluxes.shape(1)});

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

void Perturb::train()
{
  // full training_fluxes and training_k (which is 1/egienvalue)
  xt::xarray<double> training_fluxes = xt::xarray<double>::from_shape(
    {_mesh.getSize(), _training_points.shape(0)});

  xt::xarray<double> training_k =
    xt::xarray<double>::from_shape({_training_points.shape(0)});

  for (size_t i = 0; i < _training_points.shape(0); i++) {
    // set parameter
    _mesh.changeMaterial(_element_id, _training_points(i), _target_parameter);
    // find F and M matrices
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)
    // find eigenvalues and eigenvectors
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto eigenfunction = xt::linalg::eig(A);
    training_k(i) = 1 / (std::get<0>(eigenfunction)(0).real());
    xt::col(training_fluxes, i) = std::get<1>(eigenfunction)(0).real();
  }

  // reduce to PxP
  pcaReduce(training_fluxes);
}

std::pair<xt::xarray<double>, double> Perturb::calcTarget(double target_value)
{
  xt::xarray<double> target_flux =
    xt::xarray<double>::from_shape({_mesh.getSize()});
  double target_k;
  // change with specific parameter
  _mesh.changeMaterail(_element_id, target_value, _target_parameter);
  // get F and M matricies
  xt::xarray<double> F = _mesh.constructF(); //(nxn)
  xt::xarray<double> M = _mesh.constructM(); //(nxn)
  rbm::Perturb object;
  // get F_t and M_t
  xt::xarray<double> F_t = object.constructF_t(F, _training_fluxes);
  xt::xarray<double> M_t = object.constructF_t(M, _training_fluxes);
  // calculate the eigenvlue and eigenvector for target
  xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M_t), F_t);
  auto eigenfunction = xt::linalg::eig(A);
  target_k = 1 / (std::get<0>(eigenfunction)(0).real());
  xt::col(target_flux, 0) = std::get<1>(eigenfunction)(0).real();

  return std::make_pair(target_flux, target_k);
}

} // namespace rbm
