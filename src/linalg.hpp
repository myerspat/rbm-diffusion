#ifndef _LINALG_
#define _LINALG_

#include "mesh.hpp"
#include "matrix.hpp"
#include <utility>
#include <vector>

namespace linalg {

std::pair<std::vector<double>, double> powerIteration(Matrix A);

double norm(std::vector<double> vec);

double dot(std::vector<double> vec1, std::vector<double> vec2);

} // namespace linalg

#endif // !_LINALG_
