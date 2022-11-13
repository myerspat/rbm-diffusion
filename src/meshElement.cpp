#include "rbm/meshElement.hpp"

namespace mesh {

void MeshElement::initialize(double dx, double dy, Cell& cell)
{
  _dx = dx;
  _dy = dy;
  _cell = cell;
}

} // namespace mesh
