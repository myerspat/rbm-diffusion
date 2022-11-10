#include "linalg.hpp"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace linalg {

  std::pair<std::vector<double>, double> powerIteration(Matrix &M, Matrix &F) {
    // Get the inverse of M
    Matrix M_inv = M.getInverse();

    // Get A for power iteration
    Matrix A = M_inv.matmat(F);

    // Create initial guess of 0 of size n
    std::vector<double> current_e(A.size().second, 1.0);
    double current_k = 1.0;

    std::vector<double> next_e;
    double next_k;

    // Initialize relative errors
    double rel_err_e;
    double rel_err_k;

    for (size_t num_iter = 0; num_iter < max_iter; num_iter++) {
      // Calculate the new eigenvector using power iteraion
      linalg::Matrix A_inter = A * (1.0 / current_k);
      next_e = A_inter.matvec(current_e);

      // Calculate new eigenvalue
      next_k = current_k * dot(next_e, next_e) / dot(next_e, current_e);

      // Calculate eigenvector error
      rel_err_e = 0.0;
      for (size_t i = 0; i < next_e.size(); i++) {
        double err = std::abs(next_e[i] - current_e[i]) / current_e[i];
        if (rel_err_e < err)
          rel_err_e = err;
      }

      // Calculate eigenvalue error
      rel_err_k = std::abs(next_k - current_k) / current_k;

      // If convergence criterion are met return
      if (rel_err_e < e_crit && rel_err_k < k_crit)
        return std::make_pair(next_e, next_k);

      // Update for next iteration
      current_e = next_e;
      current_k = next_k;
    }

    throw std::runtime_error("Power iteration did not converge in the maximum number of iterations");
  }

  // double norm(double *vec) {}

  double dot(const std::vector<double> &vec1, const std::vector<double> &vec2) {
    assert(vec1.size() == vec2.size());
    double result = 0.0;
    for (size_t i = 0; i < vec1.size(); i++) {
      result += vec1[i] * vec2[i];
    }
    return result;
  }

} // namespace linalg
