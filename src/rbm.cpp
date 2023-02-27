#include "rbm/rbm.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xmanipulation.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xslice.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <cstdio>
#include <fstream>
#include <iomanip>
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

xt::xarray<double> Perturb::constructF_t(
  const xt::xarray<double>& F, const xt::xarray<double>& training_fluxes)
{
  // Initialize F_t and allocate space
  xt::xarray<double> F_t =
    xt::zeros<double>({training_fluxes.shape(1), training_fluxes.shape(1)});

  for (size_t i = 0; i < training_fluxes.shape(1); i++) {
    // Get flux at i
    xt::xarray<double> flux_i = xt::view(training_fluxes, xt::all(), i);
    for (size_t j = 0; j < training_fluxes.shape(1); j++) {
      // Get flux at j
      xt::xarray<double> flux_j = xt::view(training_fluxes, xt::all(), j);
      // Calculate F_t
      F_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(F, flux_j))(0);
    }
  }

  return F_t;
}

xt::xarray<double> Perturb::constructM_t(
  const xt::xarray<double>& M, const xt::xarray<double>& training_fluxes)
{
  // Initialize M_t and allocate space
  xt::xarray<double> M_t =
    xt::zeros<double>({training_fluxes.shape(1), training_fluxes.shape(1)});

  for (size_t i = 0; i < training_fluxes.shape(1); i++) {
    // Get flux at i
    xt::xarray<double> flux_i = xt::view(training_fluxes, xt::all(), i);
    for (size_t j = 0; j < training_fluxes.shape(1); j++) {
      // Get flux at j
      xt::xarray<double> flux_j = xt::view(training_fluxes, xt::all(), j);
      // Calculate F_t
      M_t(i, j) = xt::linalg::dot(flux_i, xt::linalg::dot(M, flux_j))(0);
    }
  }

  return M_t;
}

void Perturb::initialize(
  xt::xarray<double>& training_points, mesh::Mesh& mesh, size_t& element_id)
{}

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
    {_mesh.getSize(), _training_points.shape(0)});
  _training_k = xt::xarray<double>::from_shape({_training_points.shape(0)});

  for (size_t i = 0; i < _training_points.shape(0); i++) {
    // set parameter
    _mesh.changeMaterial(_element_id, _training_points(i), _target_parameter);

    // find F and M matrices
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // find eigenvalues and eigenvectors
    xt::xarray<double> A = xt::linalg::dot(xt::linalg::inv(M), F);
    auto [eigenvalues, eigenvectors] = xt::linalg::eig(A);

    // Find max eigenvalue corresponding eigenvector
    auto fundamental =
      findMaxEigen(xt::real(eigenvalues), xt::real(eigenvectors));

    _training_k(i) = fundamental.first;
    xt::col(_training_fluxes, i) = xt::abs(fundamental.second);

    // Update user
    printf("   Point %lu = %3lg => k = %6lg\n", i + 1, _training_points(i),
      _training_k(i));
  }

  // Write training data to csv file
  writePointData(
    "training.csv", _training_points, _training_k, _training_fluxes);

  // reduce to PxP
  if (_training_points.size() != 1) {
    pcaReduce(_training_fluxes);
  }

  // Ensure fluxes are normalized to 1
  for (size_t i = 0; i < _training_fluxes.shape(1); i++) {
    xt::col(_training_fluxes, i) /= xt::norm_l2(xt::col(_training_fluxes, i));
  }

  // Write PCA data
  writePCAData();
}

void Perturb::calcTargets()
{
  // Prompt user of online mode
  printf("\n Calculating Targets\n");

  _target_fluxes =
    xt::zeros<double>({_training_fluxes.shape(0), _target_points.size()});
  _target_k = xt::xarray<double>::from_shape({_target_points.size()});

  for (size_t i = 0; i < _target_points.size(); i++) {
    // change with specific parameter
    _mesh.changeMaterial(_element_id, _target_points(i), _target_parameter);

    // get F and M matricies
    xt::xarray<double> F = _mesh.constructF(); //(nxn)
    xt::xarray<double> M = _mesh.constructM(); //(nxn)

    // get F_t and M_t
    xt::xarray<double> F_t = constructF_t(F, _training_fluxes);
    xt::xarray<double> M_t = constructM_t(M, _training_fluxes);

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
    printf("   Point %lu = %3lg => k = %6lg\n", i + 1, _target_points(i),
      _target_k(i));
  }

  // Write target data
  writePointData("target.csv", _target_points, _target_k, _target_fluxes);
}

void Perturb::checkError()
{
  printf("\n Calculating Errors\n");

  auto exact_target_fluxes =
    xt::xarray<double>::from_shape({_mesh.getSize(), _target_points.size()});
  auto exact_target_k = xt::xarray<double>::from_shape({_target_points.size()});
  auto error = xt::xarray<double>::from_shape({_target_points.size()});

  for (size_t i = 0; i < _target_points.shape(0); i++) {
    // set parameter
    _mesh.changeMaterial(_element_id, _target_points(i), _target_parameter);

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
    printf(
      "   Point %lu = %3lg => k_exact = %6lg, k_rbm = %6lg, error = %6lg\n",
      i + 1, _target_points(i), exact_target_k(i), _target_k(i), error(i));
  }

  // Print data
  writePointData("error.csv", error, exact_target_k, exact_target_fluxes);
}

void Perturb::writePointData(const std::string& file_name,
  const xt::xarray<double>& points, const xt::xarray<double>& k,
  const xt::xarray<double>& fluxes)
{
  std::cout << "\n Writing data to " + file_name + "\n";

  // Open file
  std::ofstream file;
  file.open(file_name);

  // Printing point data
  file << std::setprecision(12) << points(0);
  for (size_t i = 1; i < points.size(); i++) {
    file << "," << std::setprecision(12) << points(i);
  }
  file << std::endl;

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
