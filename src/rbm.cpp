#include "rbm/rbm.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <string>
#include <tuple>
#include <utility>

namespace rbm {

void Perturb::pcaReduce(
  xt::xarray<double>& training_fluxes, xt::xarray<double>& training_k)
{}

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
    std::string absor = "absorption";
    _mesh.changeCell(_cell_id, absor, _training_points(i));
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
  pcaReduce(training_fluxes, training_k);
}

std::pair<xt::xarray<double>, double> Perturb::calcTarget(double target_value)
{
  xt::xarray<double> target_flux =
    xt::xarray<double>::from_shape({_mesh.getSize()});
  double target_k;

  return std::make_pair(target_flux, target_k);
}

} // namespace rbm
