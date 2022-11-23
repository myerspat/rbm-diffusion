#include "rbm/rbm.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <utility>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xview.hpp>

namespace rbm {

void RBM::pcaReduce(xt::xarray<double>& training_fluxes)
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

  // Reduce subspace to the first PCs with 80% of the original variance
  // TODO: Replace 0.8 with a user choice between 0 and 1 for variance they want
  // to preserve
  size_t num_pcs = 1;
  double var = 0.0;
  while (var < 0.8 * total_variance && num_pcs <= rank) {
    var += (L(num_pcs - 1) * L(num_pcs - 1));
    num_pcs++;
  }

  // Reduced U, L, and At
  xt::xarray<double> U_r = xt::view(U, xt::all(), xt::range(0, num_pcs));
  xt::xarray<double> L_r = xt::view(L, xt::range(0, num_pcs));
  xt::xarray<double> At_r =
    xt::view(At, xt::range(0, num_pcs), xt::range(0, num_pcs));

  // Calculate reduced training_fluxes
  training_fluxes = xt::linalg::dot(U_r, xt::linalg::dot(xt::diag(L_r), At_r));

  // Uncenter training_fluxes with col_means from 0 to num_pcs
  col_means = xt::view(col_means, xt::range(0, num_pcs));
  for (size_t i = 0; i < training_fluxes.shape(0); i++) {
    auto center = xt::view(training_fluxes, i, xt::all());
    center += col_means;
  }
}

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

  return M_t;
}

void PerturbAbsorption::initialize(
  xt::xarray<double>& training_points, mesh::Mesh& mesh, int& cell_id)
{}

void PerturbAbsorption::train() {}

std::pair<xt::xarray<double>, double> PerturbAbsorption::calcTarget(
  double target_value)
{
  xt::xarray<double> target_flux =
    xt::xarray<double>::from_shape({_mesh.getSize()});
  double target_k;

  return std::make_pair(target_flux, target_k);
}

} // namespace rbm
