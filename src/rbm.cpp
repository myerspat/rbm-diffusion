#include "rbm/rbm.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <utility>
#include "xtensor-blas/xlinalg.hpp"

namespace rbm {

void RBM::pcaReduce(
  xt::xarray<double>& training_fluxes, xt::xarray<double>& training_k)
{}

xt::xarray<double> RBM::constructF_t(
  const xt::xarray<double>& F, const xt::xarray<double>& training_fluxes)
{
  // Initialize F_t and allocate space
  xt::xarray<double> F_t = xt::xarray<double>::from_shape(
    {training_fluxes.shape(0), training_fluxes.shape(0)});

  return F_t;
}

xt::xarray<double> RBM::constructM_t(
  const xt::xarray<double>& M, const xt::xarray<double>& training_fluxes)
{
  // Initialize M_t and allocate space
  xt::xarray<double> M_t = xt::xarray<double>::from_shape(
    {training_fluxes.shape(0), training_fluxes.shape(0)});

    for (size_t i = 0; i < training_fluxes.shape(1); i++) {
      // Get flux at i
      auto flux_i = xt::view(training_fluxes, xt::all(), i);
      for (size_t j = 0; j < training_fluxes.shape(1); j++) {
        // Get flux at j
        auto flux_j = xt::view(training_fluxes, xt::all(), j);
        // Calculate F_t
        M_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(M, flux_j))(0);
      }
    }

  return M_t;
}

void PerturbAbsorption::initialize(
  xt::xarray<double>& training_points, mesh::Mesh& mesh, int& cell_id)
{}

void PerturbAbsorption::train() {}

std::pair<xt::xarray<double>, double> PerturbAbsorption::calcTarget(
  double target_value)
{
  xt::xarray<double> target_flux = xt::xarray<double>::from_shape({_mesh.getSize()});
  double target_k;

  return std::make_pair(target_flux, target_k);
}

} // namespace rbm
