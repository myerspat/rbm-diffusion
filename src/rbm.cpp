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

    

    for (size_t i = 0; i < training_fluxes.shape(0); i++) {
      // get row and coln fluxes 
      auto flux_i = xt::view(training_fluxes, i, xt::all());
      auto flux_j = xt::view(training_fluxes, i, xt::all());
      // find F_t
      xt::row(F_t, i) = xt::linalg::dot(flux_i, xt::linalg::dot(F, flux_j));
    }

    // for 0, M-1
    // Ft(i :)=fluxs_row(i) *[ F * fluxs_coln(i) ] 

  return F_t;
}

xt::xarray<double> RBM::constructM_t(
  const xt::xarray<double>& M, const xt::xarray<double>& training_fluxes)
{
  // Initialize M_t and allocate space
  xt::xarray<double> M_t = xt::xarray<double>::from_shape(
    {training_fluxes.shape(0), training_fluxes.shape(0)});

    for (size_t i = 0; i < training_fluxes.shape(0); i++) {
      // get row and coln fluxes 
      auto flux_i = xt::view(training_fluxes, i, xt::all());
      auto flux_j = xt::view(training_fluxes, i, xt::all());
      // find M_t
      xt::row(M_t, i) = xt::linalg::dot(flux_i, xt::linalg::dot(M, flux_j));
    }

    

  return M_t;
}

void PerturbAbsorption::initialize(
  xt::xarray<double>& training_points, mesh::Mesh& mesh, int& cell_id)
{}

void PerturbAbsorption::train() {

// full training_fluxes and training_k (which is 1/egienvalue)

//auto flux = xt::view(flux_train, i, xt::all());
//flux = xt::real(xt::col(std::get<1>(eigenfunction), 0)) * -1;

// for i = 0, training points.shape(0)-1
// change Cell( "abosoption", _traning points(i),cellid)
// F=constructF() (nxn)
// M=constructM() (nxn)
//xt::xarray<double> M = my_mesh.getM();
//  xt::xarray<double> F = my_mesh.getF();
//eig(M^(-1) *F  //tuple eigenvecotors and eigencetor (want one evigitveotor and cospiong eigvalue)

//fill traing flux / k
//pcareduce() (dont use for unit test)

//xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
// auto eigenfunction = xt::linalg::eig(A);
//k_train[i] = std::get<0>(eigenfunction)(0, 0).real();

// // Build target matricies
//   for (size_t i = 0; i < Sigma_a_train.size(); i++) {
//     for (size_t j = 0; j < Sigma_a_train.size(); j++) {
//       // Set the dot of previous result and forward flux
//       auto flux_i = xt::view(flux_train, i, xt::all());
//       auto flux_j = xt::view(flux_train, i, xt::all());
//       M_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(M, flux_j))(0, 0);
//       F_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(F, flux_j))(0, 0);
//     }
//   }

}

std::pair<xt::xarray<double>, double> PerturbAbsorption::calcTarget(
  double target_value)
{
  xt::xarray<double> target_flux = xt::xarray<double>::from_shape({_mesh.getSize()});
  double target_k;

  return std::make_pair(target_flux, target_k);
}

} // namespace rbm
