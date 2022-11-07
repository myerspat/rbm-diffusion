#ifndef _LINALG_
#define _LINALG_

#include "matrix.hpp"
#include <utility>
#include <vector>

namespace linalg {

constexpr double e_crit = 1.0e-6;   // Eigenvector convergence criterion
constexpr double k_crit = 1.0e-7;   // Eigenvalue convergence criterion
constexpr size_t max_iter = 500000; // Max number of iterations

std::pair<std::vector<double>, double> powerIteration(Matrix &M, Matrix &F);

double norm(std::vector<double> &vec);

double dot(std::vector<double> &vec1, std::vector<double> &vec2);

} // namespace linalg

#endif // !_LINALG_
