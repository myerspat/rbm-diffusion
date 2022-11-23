#ifndef _RBM_
#define _RBM_

#include "rbm/mesh.hpp"

namespace rbm {

class RBM {
public:
  //=============================================================
  // Constructors / Destructor
  RBM() {};

  //=============================================================
  // Methods
  // Initialize rbm with the training points
  virtual void initialize(
    xt::xarray<double>& training_points, mesh::Mesh& mesh, int& cell_id) = 0;

  // Create the training subspace
  virtual void train() = 0;

  // Calculate the eigenvalue and eigenvector for the target value
  virtual std::pair<xt::xarray<double>, double> calcTarget(
    double target_value) = 0;

  // Reduce the subspace using principle component analysis
  void pcaReduce(xt::xarray<double>& training_fluxes);

  // Construct the target fission matrix
  xt::xarray<double> constructF_t(
    const xt::xarray<double>& F, const xt::xarray<double>& training_fluxes);

  // Construct the target migration matrix
  xt::xarray<double> constructM_t(
    const xt::xarray<double>& M, const xt::xarray<double>& training_fluxes);

  // Set the number of PCs to keep in reduction
  void setNumPCs(const size_t num_pcs) { _num_pcs = num_pcs; }

private:
  //=============================================================
  // Data
  // TODO: Make _num_pcs editable for the user in the XML interface
  size_t _num_pcs = 3; // number of PCs to keep when reducing, defualts to 3
};

class PerturbAbsorption : public RBM {
public:
  //=============================================================
  // Constructors / Destructor
  PerturbAbsorption() : RBM() {};
  PerturbAbsorption(
    xt::xarray<double>& training_points, mesh::Mesh& mesh, int& cell_id)
    : RBM(), _training_points(training_points), _mesh(mesh),
      _cell_id(cell_id) {};

  //=============================================================
  // Methods
  // Initialize rbm with the training points
  void initialize(
    xt::xarray<double>& training_points, mesh::Mesh& mesh, int& cell_id);

  // Create the training subspace
  void train();

  // Calculate the eigenvalue and eigenvector for the target value
  std::pair<xt::xarray<double>, double> calcTarget(double target_value);

private:
  //=============================================================
  // Data
  xt::xarray<double> _training_points; // 1D array of training paramaters
  xt::xarray<double> _training_fluxes; // 2D array of training fluxes
  xt::xarray<double> _training_k;      // 1D array of training k
  mesh::Mesh _mesh;                    // Mesh for the problem
  int _cell_id; // Cell within mesh that will be perturbed
};

} // namespace rbm

#endif // !_RBM_
