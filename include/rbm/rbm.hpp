#ifndef _RBM_
#define _RBM_

#include "rbm/mesh.hpp"
#include "rbm/parameter.hpp"
#include <vector>

namespace rbm {

class Perturb {
public:
  //=============================================================
  // Constructors / Destructor
  Perturb() {};
  Perturb(const std::vector<Parameter>& parameters, const mesh::Mesh& mesh,
    const bool& constructD)
    : _parameters(parameters), _mesh(mesh), _constructD(constructD) {};

  //=============================================================
  // Methods
  // Create the training subspace
  void train();

  // Load target data if already calculated
  void loadFluxes(const std::string& path);

  // Calculate the eigenvalue and eigenvector for the target value
  void calcTargets();

  // Reduce the subspace using principle component analysis
  void pcaReduce(xt::xarray<double>& training_fluxes);

  std::pair<double, xt::xarray<double>> findMaxEigen(
    xt::xarray<double> eigenvalues, xt::xarray<double> eigenvectors);

  // Check errors
  void checkError();

  // Write data to file_name
  void writeTrainingData(const std::string& file_name,
    const xt::xarray<double>& k, const xt::xarray<double>& fluxes);

  // Write data to file_name
  void writeTargetData(const std::string& file_name,
    const xt::xarray<double>& k, const xt::xarray<double>& fluxes);

  // Multiply fluxes by matrix
  xt::xarray<double> constructMatrix(
    const xt::xarray<double>& matrix, const xt::xarray<double> training_fluxes);

  // Write PCA data to reduced.csv
  void writePCAData();

  // FOPT verification function
  void foptCalcTargets();

  //=============================================================
  // Setters
  // Set the number of PCs to keep in reduction
  void setNumPCs(const size_t num_pcs) { _num_pcs = num_pcs; }

  // Set bool for EIMs use
  void setEIM(const bool& use_eim) { _use_eim = use_eim; }

  void setPrecompute(const bool& precompute) { _precompute = precompute; }

  void setFOPT(const bool& fopt) { _fopt = fopt; }

  std::size_t getNumTraining()
  {
    return _parameters[0].getTrainingPoints().size();
  }
  std::size_t getNumTarget() { return _parameters[0].getTargetPoints().size(); }

  bool getFOPT() { return _fopt; }

private:
  //=============================================================
  // Data
  std::vector<Parameter> _parameters; // Vector of perturbed parameters
  mesh::Mesh _mesh;                   // Mesh of the problem

  size_t _num_pcs = 3; // number of PCs to keep when reducing, defualts to 3
  xt::xarray<double> _variance; // Variance of PCA data

  bool _fopt = false;
  bool _constructD;
  bool _precompute = true;
  bool _use_eim = true;
  xt::xarray<double> _A_t;
  xt::xarray<double> _D_t;

  xt::xarray<double> _training_fluxes; // 2D array of training fluxes
  xt::xarray<double> _training_k;      // 1D array of training k
  xt::xarray<double> _target_fluxes;   // 2D array of target fluxes
  xt::xarray<double> _target_k;        // 2D array of target k
};

} // namespace rbm

#endif // !_RBM_
