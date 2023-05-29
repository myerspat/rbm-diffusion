#include "rbm/rbm.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xmanipulation.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xslice.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xview.hpp>

namespace rbm {

void Perturb::pcaReduce(xt::xarray<double>& training_fluxes)
{
  // Rank of flux data
  size_t rank = xt::linalg::matrix_rank(training_fluxes);
  if (rank < _num_pcs) {
    printf("\n Warning: Rank of training fluxes matrix is less than the number "
           "of PCs, finding number of PCs based on variance\n");
    _num_pcs = 0;
  }

  // Singular value decomposition
  auto [U, L, At] = xt::linalg::svd(training_fluxes, false);

  // Flip eigenvectors' sign to ensure deterministic output
  auto max_abs_col = xt::argmax(xt::abs(U), 0);
  xt::xarray<double> signs(max_abs_col.size());
  for (size_t i = 0; i < max_abs_col.size(); i++) {
    auto sign = U(max_abs_col(i), i) / std::abs(U(max_abs_col(i), i));

    xt::col(U, i) *= sign;
    xt::row(At, i) *= sign;
  }

  // Variance
  _variance = xt::square(L) / (training_fluxes.shape(0) - 1);

  // Calculate total variance
  double total_variance = xt::sum(_variance)(0);

  // If _num_pcs is zero use variance to determine subspace size
  // else use the number defined by the user
  double var = 0.0;
  if (_num_pcs == 0) {
    while (var < 0.9 * total_variance) {
      var += _variance(_num_pcs);
      _num_pcs++;
    }

    std::cout << "\n Preserved " << var / total_variance * 100
              << " % of the variance with " << _num_pcs << " PCs\n\n";
  } else {
    // If the _num_pcs was not zero reduced to user defined number
    for (size_t i = 0; i < _num_pcs; i++) {
      var += _variance(i);
    }

    // Throw warning if the number of PCs preserved is less than 90% of the
    // total variance
    if (var < 0.9 * total_variance) {
      printf("\n Warning: Subspace was reduced to %lu PCs which has only %lg "
             "percent of the total variance\n",
        _num_pcs, var / total_variance * 100);
    }
  }

  training_fluxes = xt::view(U, xt::all(), xt::range(0, _num_pcs));
}

void Perturb::loadFluxes(const std::string& path)
{
  std::ifstream in_file;
  in_file.open(path);
  auto data = xt::load_csv<double>(in_file, ',', _parameters.size());
  in_file.close();

  _training_k = xt::view(data, 0, xt::all());
  _training_fluxes =
    xt::view(data, xt::range(1, xt::placeholders::_), xt::all());

  // reduce to PxP
  if (_parameters[0].getTrainingPoints().size() != 1) {
    pcaReduce(_training_fluxes);
  }

  // Ensure fluxes are normalized to 1
  for (size_t i = 0; i < _training_fluxes.shape(1); i++) {
    xt::col(_training_fluxes, i) /= xt::norm_l2(xt::col(_training_fluxes, i));
  }

  // Write PCA data
  writePCAData();

  // Precompute
  if (_precompute) {
    // Allocate space for matrix
    _A_t = xt::xarray<double>::from_shape({_training_fluxes.shape(1),
      _training_fluxes.shape(1), _mesh.getRegions().size()});

    // For each region precompute flux weighting
    for (size_t i = 0; i < _mesh.getRegions().size(); i++) {
      auto A_t_i = xt::view(_A_t, xt::all(), xt::all(), i);
      A_t_i = constructMatrix(
        xt::diag(_mesh.getRegions()[i].getMask()), _training_fluxes);
    }

    // Precompute diffusion region if diffusion coefficients aren't perturbed
    if (!_constructD) {
      _D_t = constructMatrix(_mesh.getDMatrix(), _training_fluxes);
    }
  } else {
    // Allocate space for matrix
    _A_t = xt::xarray<double>::from_shape(
      {_mesh.getSize(), _mesh.getSize(), _mesh.getRegions().size()});

    // For each region precompute flux weighting
    for (size_t i = 0; i < _mesh.getRegions().size(); i++) {
      auto A_t_i = xt::view(_A_t, xt::all(), xt::all(), i);
      A_t_i = xt::diag(_mesh.getRegions()[i].getMask());
    }

    // Get diffusion matrix
    if (!_constructD) {
      _D_t = _mesh.getDMatrix();
    }
  }
}

std::pair<double, xt::xarray<double>> Perturb::findMaxEigen(
  xt::xarray<double> eigenvalues, xt::xarray<double> eigenvectors)
{
  size_t idx = 0;
  double max_eigenvalue = 0.0;
  for (size_t i = 0; i < eigenvalues.size(); i++) {
    if (eigenvalues(i) > max_eigenvalue) {
      max_eigenvalue = eigenvalues(i);
      idx = i;
    }
  }

  return std::make_pair(eigenvalues(idx), xt::col(eigenvectors, idx));
}

void Perturb::train()
{
  // Prompt user of training begun
  printf("\n Begin training\n");

  // Initialize training_fluxes and training_k
  _training_fluxes = xt::xarray<double>::from_shape(
    {_mesh.getSize(), _parameters[0].getTrainingPoints().size()});
  _training_k =
    xt::xarray<double>::from_shape({_parameters[0].getTrainingPoints().size()});

  for (size_t i = 0; i < _parameters[0].getTrainingPoints().size(); i++) {
    // Change material properties for each region
    printf("   Point %lu = ", i + 1);
    for (Parameter& parameter : _parameters) {
      _mesh.changeMaterial(parameter.getID(), parameter.getTrainingPoints()(i),
        parameter.getMaterialProperty());
      printf("%3lg ", parameter.getTrainingPoints()(i));
    }

    if (_constructD == true) {
      _mesh.updateD(_parameters, i);
    }

    // Find F and M matrices
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // Find eigenvalues and eigenvectors
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

    // Find max eigenvalue corresponding eigenvector
    auto fundamental =
      findMaxEigen(xt::real(eigenvalues), xt::real(eigenvectors));

    _training_k(i) = fundamental.first;
    xt::col(_training_fluxes, i) = xt::abs(fundamental.second);

    // Update user
    printf("=> k = %6lg\n", _training_k(i));
  }

  // Write training data to csv file
  writeTrainingData("training.csv", _training_k, _training_fluxes);

  // reduce to PxP
  if (_parameters[0].getTrainingPoints().size() != 1) {
    pcaReduce(_training_fluxes);
  }

  // Ensure fluxes are normalized to 1
  for (size_t i = 0; i < _training_fluxes.shape(1); i++) {
    xt::col(_training_fluxes, i) /= xt::norm_l2(xt::col(_training_fluxes, i));
  }

  // Write PCA data
  writePCAData();

  // Precompute
  if (_precompute) {
    // Allocate space for matrix
    _A_t = xt::xarray<double>::from_shape({_training_fluxes.shape(1),
      _training_fluxes.shape(1), _mesh.getRegions().size()});

    // For each region precompute flux weighting
    for (size_t i = 0; i < _mesh.getRegions().size(); i++) {
      auto A_t_i = xt::view(_A_t, xt::all(), xt::all(), i);
      A_t_i = constructMatrix(
        xt::diag(_mesh.getRegions()[i].getMask()), _training_fluxes);
    }

    // Precompute diffusion region if diffusion coefficients aren't perturbed
    if (!_constructD) {
      _D_t = constructMatrix(_mesh.getDMatrix(), _training_fluxes);
    }
  } else {
    // Allocate space for matrix
    _A_t = xt::xarray<double>::from_shape(
      {_mesh.getSize(), _mesh.getSize(), _mesh.getRegions().size()});

    // For each region precompute flux weighting
    for (size_t i = 0; i < _mesh.getRegions().size(); i++) {
      auto A_t_i = xt::view(_A_t, xt::all(), xt::all(), i);
      A_t_i = xt::diag(_mesh.getRegions()[i].getMask());
    }

    // Get diffusion matrix
    if (!_constructD) {
      _D_t = _mesh.getDMatrix();
    }
  }
}

xt::xarray<double> Perturb::constructMatrix(
  const xt::xarray<double>& matrix, const xt::xarray<double> training_fluxes)
{
  xt::xarray<double> results_matrix =
    xt::zeros<double>({training_fluxes.shape(1), training_fluxes.shape(1)});

  for (size_t i = 0; i < training_fluxes.shape(1); i++) {
    // Get flux at i
    xt::xarray<double> flux_i = xt::view(training_fluxes, xt::all(), i);
    for (size_t j = 0; j < training_fluxes.shape(1); j++) {
      // Get flux at j
      xt::xarray<double> flux_j = xt::view(training_fluxes, xt::all(), j);
      // Calculate F_t
      results_matrix(i, j) =
        xt::linalg::dot(flux_i, xt::linalg::dot(matrix, flux_j))(0);
    }
  }

  return results_matrix;
}

void Perturb::calcTargets()
{
  // Prompt user of online mode
  printf("\n Calculating Targets\n");

  _target_fluxes = xt::zeros<double>(
    {_mesh.getSize(), _parameters[0].getTargetPoints().size()});
  _target_k =
    xt::xarray<double>::from_shape({_parameters[0].getTargetPoints().size()});

  for (size_t i = 0; i < _parameters[0].getTargetPoints().size(); i++) {
    // Change with specific parameter
    printf("   Point %lu = ", i + 1);
    for (Parameter& parameter : _parameters) {
      _mesh.changeMaterial(parameter.getID(), parameter.getTargetPoints()(i),
        parameter.getMaterialProperty());
      printf("%3lg ", parameter.getTargetPoints()(i));
    }

    // If a diffusion coefficient is perturbed use EIM or don't use affine
    // decomposition
    if (_constructD) {
      if (_use_eim) {
        throw std::runtime_error("EIM is not yet supported");
      } else {
        _mesh.updateD(_parameters, i);
        if (_precompute) {
          _D_t = constructMatrix(_mesh.getDMatrix(), _training_fluxes);
        }
      }
    }

    xt::xarray<double> F_t, M_t;
    if (_precompute) {
      F_t = xt::zeros<double>({_A_t.shape(0), _A_t.shape(1)});
      M_t = xt::zeros<double>({_A_t.shape(0), _A_t.shape(1)});

      for (size_t i = 0; i < _mesh.getRegions().size(); i++) {
        const auto& material = _mesh.getRegions()[i].getMat();
        const auto& A_t = xt::view(_A_t, xt::all(), xt::all(), i);
        F_t += material.getNuFission() * A_t;
        M_t += material.getAbsorption() * A_t;
      }
      M_t += _D_t;
    } else {
      F_t = constructMatrix(_mesh.constructF(), _training_fluxes);
      M_t = constructMatrix(_mesh.constructM(), _training_fluxes);
    }

    // calculate the eigenvlue and eigenvector for target
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M_t), F_t);
    auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

    // Find max eigenvalue corresponding eigenvector
    auto fundamental =
      findMaxEigen(xt::real(eigenvalues), xt::real(eigenvectors));

    // Get target eigenvalue
    _target_k(i) = fundamental.first;

    // Calculate approximate target flux using the sum of ritz vector *
    // _training_fluxes
    for (size_t j = 0; j < _training_fluxes.shape(1); j++) {
      xt::col(_target_fluxes, i) +=
        (fundamental.second(j) * xt::col(_training_fluxes, j));
    }

    // Normlize target fluxes
    xt::col(_target_fluxes, i) = xt::abs(xt::col(_target_fluxes, i)) /
                                 xt::norm_l2(xt::col(_target_fluxes, i));

    // Update user
    printf("=> k = %6lg\n", _target_k(i));
  }

  // Write target data
  writeTargetData("target.csv", _target_k, _target_fluxes);
}

void Perturb::foptCalcTargets()
{
  assert(_training_fluxes.shape(1) == 1 && _training_k.size() == 1);

  // Prompt user of online mode
  printf("\n Calculating Targets\n");

  _target_fluxes = xt::zeros<double>(
    {_mesh.getSize(), _parameters[0].getTargetPoints().size()});
  _target_k =
    xt::xarray<double>::from_shape({_parameters[0].getTargetPoints().size()});

  const double k1 = _training_k(0);
  const xt::xarray<double> flux1 = xt::view(_training_fluxes, xt::all(), 0);

  const xt::xarray<double> F1 = _mesh.constructF(); //(nxn)
  const xt::xarray<double> M1 = _mesh.constructM(); //(nxn)

  for (size_t i = 0; i < _parameters[0].getTargetPoints().size(); i++) {
    printf("   Point %lu = ", i + 1);
    for (Parameter& parameter : _parameters) {
      _mesh.changeMaterial(parameter.getID(), parameter.getTargetPoints()(i),
        parameter.getMaterialProperty());
      printf("%3lg ", parameter.getTargetPoints()(i));
    }

    if (_constructD == true) {
      _mesh.updateD(_parameters, i);
    }

    // Find F and M matrices
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // Get target eigenvalue using FOPT equation
    _target_k(i) =
      1 / ((xt::linalg::dot(
              flux1, xt::linalg::dot((M - M1) - k1 * (F - F1), flux1)) /
             xt::linalg::dot(flux1, xt::linalg::dot(F, flux1)))(0) +
            1 / k1);

    // Update user
    printf("=> k = %6lg\n", _target_k(i));
  }

  // Write target data
  writeTargetData("target.csv", _target_k, _target_fluxes);
}

void Perturb::checkError()
{
  printf("\n Calculating Errors\n");

  auto exact_target_fluxes = xt::xarray<double>::from_shape(
    {_mesh.getSize(), _parameters[0].getTargetPoints().size()});
  auto exact_target_k =
    xt::xarray<double>::from_shape({_parameters[0].getTargetPoints().size()});
  auto error =
    xt::xarray<double>::from_shape({_parameters[0].getTargetPoints().size()});

  for (size_t i = 0; i < _parameters[0].getTargetPoints().size(); i++) {
    printf("   Point %lu = ", i + 1);
    // Set parameter
    for (Parameter& parameter : _parameters) {
      _mesh.changeMaterial(parameter.getID(), parameter.getTargetPoints()(i),
        parameter.getMaterialProperty());
      printf("%3lg ", parameter.getTargetPoints()(i));
    }

    if (_constructD == true) {
      _mesh.updateD(_parameters, i);
    }

    // find F and M matrices
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // find eigenvalues and eigenvectors
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

    // Find max eigenvalue corresponding eigenvector
    auto fundamental =
      findMaxEigen(xt::real(eigenvalues), xt::real(eigenvectors));

    exact_target_k(i) = fundamental.first;
    xt::col(exact_target_fluxes, i) = xt::abs(fundamental.second);

    // Calculate relative error
    error(i) = (exact_target_k(i) - _target_k(i)) / exact_target_k(i);

    // Update user
    printf("=> k_exact = %6lg, k_rbm = %6lg, error = %6lg\n", exact_target_k(i),
      _target_k(i), error(i));
  }

  // Print data
  writeTargetData("error.csv", exact_target_k, exact_target_fluxes);
}

void Perturb::writeTrainingData(const std::string& file_name,
  const xt::xarray<double>& k, const xt::xarray<double>& fluxes)
{
  std::cout << "\n Writing data to " + file_name + "\n";

  // Open file
  std::ofstream file;
  file.open(file_name);

  // Printing point data
  for (Parameter& parameter : _parameters) {
    xt::xarray<double> points = parameter.getTrainingPoints();
    file << std::setprecision(12) << points(0);
    for (size_t i = 1; i < points.size(); i++) {
      file << "," << std::setprecision(12) << points(i);
    }
    file << std::endl;
  }

  // Printing eigenvalue data
  file << std::setprecision(12) << k(0);
  for (size_t i = 1; i < k.size(); i++) {
    file << "," << std::setprecision(12) << k(i);
  }
  file << std::endl;

  // Printing flux data
  for (size_t i = 0; i < fluxes.shape(0); i++) {
    file << std::setprecision(12) << fluxes(i, 0);
    for (size_t j = 1; j < fluxes.shape(1); j++) {
      file << "," << std::setprecision(12) << fluxes(i, j);
    }
    file << std::endl;
  }

  file.close();
}

void Perturb::writeTargetData(const std::string& file_name,
  const xt::xarray<double>& k, const xt::xarray<double>& fluxes)
{
  std::cout << "\n Writing data to " + file_name + "\n";

  // Open file
  std::ofstream file;
  file.open(file_name);

  // Printing point data
  for (Parameter& parameter : _parameters) {
    xt::xarray<double> points = parameter.getTargetPoints();
    file << std::setprecision(12) << points(0);
    for (size_t i = 1; i < points.size(); i++) {
      file << "," << std::setprecision(12) << points(i);
    }
    file << std::endl;
  }

  // Printing eigenvalue data
  file << std::setprecision(12) << k(0);
  for (size_t i = 1; i < k.size(); i++) {
    file << "," << std::setprecision(12) << k(i);
  }
  file << std::endl;

  // Printing flux data
  for (size_t i = 0; i < fluxes.shape(0); i++) {
    file << std::setprecision(12) << fluxes(i, 0);
    for (size_t j = 1; j < fluxes.shape(1); j++) {
      file << "," << std::setprecision(12) << fluxes(i, j);
    }
    file << std::endl;
  }

  file.close();
}

void Perturb::writePCAData()
{
  std::cout << " Writing data to reduced.csv\n";

  // Open file
  std::ofstream file;
  file.open("reduced.csv");

  // Printing variance data
  file << std::setprecision(12) << _variance(0);
  for (size_t i = 1; i < _variance.size(); i++) {
    file << "," << std::setprecision(12) << _variance(i);
  }
  file << std::endl;

  // Printing flux data
  for (size_t i = 0; i < _training_fluxes.shape(0); i++) {
    file << std::setprecision(12) << _training_fluxes(i, 0);
    for (size_t j = 1; j < _training_fluxes.shape(1); j++) {
      file << "," << std::setprecision(12) << _training_fluxes(i, j);
    }
    file << std::endl;
  }

  file.close();
}

} // namespace rbm
