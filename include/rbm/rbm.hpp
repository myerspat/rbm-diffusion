#ifndef _RBM_
#define _RBM_

#include "rbm/mesh.hpp"
#include "rbm/rbmEnums.hpp"

namespace rbm {

class Perturb {
public:
  //=============================================================
  // Constructors / Destructor
  Perturb() {};
  Perturb(xt::xarray<double>& training_points,
    xt::xarray<double>& target_points, mesh::Mesh& mesh, size_t& element_id,
    Parameter& target_parameter)
    : _training_points(training_points), _target_points(target_points),
      _mesh(mesh), _element_id(element_id),
      _target_parameter(target_parameter) {};

  //=============================================================
  // Methods
  // Initialize rbm with the training points
  void initialize(
    xt::xarray<double>& training_points, mesh::Mesh& mesh, size_t& element_id);

  // Create the training subspace
  void train();

  // Calculate the eigenvalue and eigenvector for the target value
  void calcTargets();

  // Reduce the subspace using principle component analysis
  void pcaReduce(xt::xarray<double>& training_fluxes);

  std::pair<double, xt::xarray<double>> findMaxEigen(
    xt::xarray<double> eigenvalues, xt::xarray<double> eigenvectors);

  // Check errors
  void checkError();

  // Construct the target fission matrix
  xt::xarray<double> constructF_t(
    const xt::xarray<double>& F, const xt::xarray<double>& training_fluxes);

  // Construct the target migration matrix
  xt::xarray<double> constructM_t(
    const xt::xarray<double>& M, const xt::xarray<double>& training_fluxes);

  // Set the number of PCs to keep in reduction
  void setNumPCs(const size_t num_pcs) { _num_pcs = num_pcs; }

  //=============================================================
  // Getters
  const size_t getNumTraining() const { return _training_points.size(); };
  const size_t getNumTarget() const { return _target_points.size(); };

private:
  //=============================================================
  // Data
  xt::xarray<double> _training_points; // 1D array of training paramater values
  xt::xarray<double> _training_fluxes; // 2D array of training fluxes
  xt::xarray<double> _training_k;      // 1D array of training k
  xt::xarray<double> _target_points;   // 1D array of target parameter values
  xt::xarray<double> _target_fluxes;   // 2D array of target fluxes
  xt::xarray<double> _target_k;        // 2D array of target k
  mesh::Mesh _mesh;                    // Mesh for the problem
  size_t _num_pcs = 3; // number of PCs to keep when reducing, defualts to 3
  size_t _element_id;  // Element within mesh that will be perturbed
  Parameter _target_parameter; // Perturbed parameter type (absorption, D,
                               // nu_fission)
};

} // namespace rbm

#endif // !_RBM_
