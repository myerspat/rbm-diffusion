#include "rbm/rbm.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <string>
#include <tuple>
#include <utility>

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
  // change with specific parameter
  std::string absor = "absorption";
  _mesh.changeCell(_cell_id, absor, target_value);
  // get F and M matricies
  xt::xarray<double> F = _mesh.constructF(); //(nxn)
  xt::xarray<double> M = _mesh.constructM(); //(nxn)
  rbm::PerturbAbsorption object;
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
