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

  return M_t;
}

void PerturbAbsorption::initialize(
  xt::xarray<double>& training_points, mesh::Mesh& mesh, int& cell_id)
{}

void PerturbAbsorption::train() {
// full training_fluxes and training_k (which is 1/egienvalue)
  xt::xarray<double> training_fluxes = xt::xarray<double>::from_shape({_mesh.getSize(), _training_points.shape(0)});

  xt::xarray<double> training_k = xt::xarray<double>::from_shape({_training_points.shape(0)});

  for (size_t i = 0; i < _training_points.shape(0)-1; i++){
    _mesh.changeCell(_cell_id, "absorption",  _training_points(i))
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M =  _mesh.constructM(); //(nxn)
  //xt::xarray<double> M = my_mesh.getM();
  //  xt::xarray<double> F = my_mesh.getF();
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto eigenfunction = xt::linalg::eig(A);
    training_k[i] = 1/(std::get<0>(eigenfunction)(i).real());
    xt::row(training_fluxes, i) = std::get<1>(eigenfunction)(i).real();

  }


  pcaReduce();

}

std::pair<xt::xarray<double>, double> PerturbAbsorption::calcTarget(
  double target_value)
{
  xt::xarray<double> target_flux = xt::xarray<double>::from_shape({_mesh.getSize()});
  double target_k;

  return std::make_pair(target_flux, target_k);
}

} // namespace rbm
