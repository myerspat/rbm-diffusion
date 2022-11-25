#include "rbm/rbm.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <utility>
#include "xtensor-blas/xlinalg.hpp"
#include <tuple>
#include <string>

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
  xt::xarray<double> target_flux = xt::xarray<double>::from_shape({_mesh.getSize()});
  double target_k;
	std::string absor = "absorption";
     _mesh.changeCell(_cell_id, absor,  target_value);
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M =  _mesh.constructM(); //(nxn)
    rbm::PerturbAbsorption object;
    xt::xarray<double> F_t = object.constructF_t(F, _training_fluxes);
    xt::xarray<double> M_t = object.constructF_t(M, _training_fluxes);
  //xt::xarray<double> M = my_mesh.getM();
  //  xt::xarray<double> F = my_mesh.getF();
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M_t), F_t);
    auto eigenfunction = xt::linalg::eig(A);
    target_k = 1/(std::get<0>(eigenfunction)(0).real());
    xt::col(target_flux, 0) = std::get<1>(eigenfunction)(0).real();

  return std::make_pair(target_flux, target_k);
}

} // namespace rbm
