#ifndef _MESH_ELEMENT_
#define _MESH_ELEMENT_

#include "rbm/material.hpp"
#include <assert.h>
#include <utility>

namespace mesh {

class MeshElement {
public:
  //=========================================================================================
  // Constructors / Destructor
  MeshElement() {};
  MeshElement(const material::Material& mat, const double& lx, const double& ly,
    const std::size_t& id, const std::pair<std::size_t, std::size_t>& idx_row,
    const std::pair<std::size_t, std::size_t>& idx_col)
    : _idx_row(idx_row), _idx_col(idx_col), _lx(lx), _ly(ly), _id(id), _mat(mat)
  {
    assert(checkIndexing(idx_row, idx_col));
  };

  //=========================================================================================
  // Methods

  bool checkIndexing(const std::pair<std::size_t, std::size_t>& idx_row,
    const std::pair<std::size_t, std::size_t>& idx_col);

  //=========================================================================================

  // Getters
  const material::Material& getMaterial() const { return _mat; };
  const std::size_t& getID() const { return _id; };
  const std::pair<std::size_t, std::size_t>& getRowIdx() const
  {
    return _idx_row;
  };
  const std::pair<std::size_t, std::size_t>& getColIdx() const
  {
    return _idx_col;
  };
  const double& getLX() const { return _lx; };
  const double& getLY() const { return _ly; };

  //=========================================================================================
  // Setters
  void setParameter(
    const double& new_value, const material::Property& target_parameter)
  {
    _mat.setParameter(new_value, target_parameter);
  };

private:
  //=========================================================================================
  // Data
  material::Material _mat;
  std::pair<std::size_t, std::size_t> _idx_row, _idx_col;
  double _lx;
  double _ly;
  std::size_t _id;
};

} // namespace mesh

#endif // !_MESH_ELEMENT_
