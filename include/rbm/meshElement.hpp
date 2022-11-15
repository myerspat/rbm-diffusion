#ifndef _MESH_ELEMENT_
#define _MESH_ELEMENT_

#include "rbm/cell.hpp"

namespace mesh {

class MeshElement {
public:
  //=========================================================================================
  // Constructors / Destructor
  MeshElement() {};
  MeshElement(double dx, double dy, Cell& cell)
    : _dx(dx), _dy(dy), _cell(cell) {};

  //=========================================================================================
  // Methods
  void initialize(double dx, double dy, Cell& cell);

  //=========================================================================================
  // Getters
  double dx() { return _dx; };
  double dy() { return _dy; };
  Cell& cell() { return _cell; };

private:
  //=========================================================================================
  // Data
  double _dx;
  double _dy;
  Cell _cell;
};

} // namespace mesh

#endif // !_MESH_ELEMENT_
