#include "rbm/rbm.hpp"
#include "unit_test_framework.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xview.hpp"
#include <fstream>
#include <tuple>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xcsv.hpp>

TEST(test_pcaReduce_1)
{
  // Read in Iris data
  std::ifstream file_in;
  file_in.open("../../tests/supporting/iris.csv");
  xt::xarray<double> data = xt::load_csv<double>(file_in);
  file_in.close();

  // Initialize object for PCA
  rbm::Perturb object;

  // Set size of reduced subspace to 2 PCs
  object.setNumPCs(2);

  // Run pcaReduce on Iris data
  object.pcaReduce(data);

  // Read in the expected uncentered data from sklearn.decomposition.PCA
  file_in.open("../../tests/supporting/iris_pcas.csv");
  xt::xarray<double> expected_data = xt::load_csv<double>(file_in);
  file_in.close();

  // Calculate rank of reduced data
  size_t rank = xt::linalg::matrix_rank(data);
  size_t expected_rank = xt::linalg::matrix_rank(expected_data);
  ASSERT_EQUAL(rank, expected_rank);

  // Assert the local_pcs are the same generated by sklearn.decomposition.PCA
  for (size_t i = 0; i < data.shape(0); i++) {
    for (size_t j = 0; j < data.shape(1); j++) {
      ASSERT_ALMOST_EQUAL(data(i, j), expected_data(i, j), 1e-13);
    }
  }
}

TEST(ContructM_t_test)
{
  xt::xarray<double> M = xt::xarray<double>::from_shape({3, 3});

  M(0, 0) = 2.2;
  M(1, 0) = 0.0;
  M(2, 0) = 0.0;
  M(0, 1) = 0.0;
  M(1, 1) = 1.3;
  M(2, 1) = 0.0;
  M(0, 2) = 0.0;
  M(1, 2) = 0.0;
  M(2, 2) = 3.5;

  xt::xarray<double> Flux = xt::xarray<double>::from_shape({3, 2});

  Flux(0, 0) = 0.2;
  Flux(1, 0) = 0.94;
  Flux(2, 0) = 0.03;
  Flux(0, 1) = 0.3;
  Flux(1, 1) = 0.230;
  Flux(2, 1) = 0.540;

  rbm::Perturb object;
  xt::xarray<double> M_t = object.constructM_t(M, Flux);

  ASSERT_ALMOST_EQUAL(M_t(0, 0), 1.2398, 0.0001);
  ASSERT_ALMOST_EQUAL(M_t(1, 0), 0.4698, 0.0001);
  ASSERT_ALMOST_EQUAL(M_t(0, 1), 0.4698, 0.0001);
  ASSERT_ALMOST_EQUAL(M_t(1, 1), 1.2874, 0.0001);
}

TEST(ContructF_t_test)
{
  xt::xarray<double> F = xt::xarray<double>::from_shape({3, 3});

  F(0, 0) = 2.2;
  F(1, 0) = 0.0;
  F(2, 0) = 0.0;
  F(0, 1) = 0.0;
  F(1, 1) = 1.3;
  F(2, 1) = 0.0;
  F(0, 2) = 0.0;
  F(1, 2) = 0.0;
  F(2, 2) = 3.5;

  xt::xarray<double> Flux = xt::xarray<double>::from_shape({3, 2});

  Flux(0, 0) = 0.2;
  Flux(1, 0) = 0.94;
  Flux(2, 0) = 0.03;
  Flux(0, 1) = 0.3;
  Flux(1, 1) = 0.230;
  Flux(2, 1) = 0.540;

  rbm::Perturb object;
  xt::xarray<double> F_t = object.constructF_t(F, Flux);

  ASSERT_ALMOST_EQUAL(F_t(0, 0), 1.2398, 0.0001);
  ASSERT_ALMOST_EQUAL(F_t(1, 0), 0.4698, 0.0001);
  ASSERT_ALMOST_EQUAL(F_t(0, 1), 0.4698, 0.0001);
  ASSERT_ALMOST_EQUAL(F_t(1, 1), 1.2874, 0.0001);
}

TEST_MAIN();
