#ifndef _REGION_
#define _REGION_

#include "rbm/material.hpp"
#include "rbm/meshElement.hpp"
#include "rbm/parameter.hpp"
#include <xtensor/xarray.hpp>

namespace mesh {

class Region {
private:
  //========================================================================
  // Data
  std::size_t _id;
  material::Material _mat;

  // Affine decomposition masks
  xt::xarray<double> _mask;

public:
  //========================================================================
  // Constructors
  Region() {};
  Region(const xt::xarray<MeshElement>& fine_grid, const size_t& xN_fine,
    const size_t& yN_fine, const std::size_t& id,
    const material::Material& mat);

  //========================================================================
  // Methods
  void changeMaterial(
    const double& new_value, const material::Property& target_property);

  //========================================================================
  // Getters
  std::size_t getID() const { return _id; };
  material::Material getMat() const { return _mat; };
  xt::xarray<double> getMask() const { return _mask; };
};

} // namespace mesh

#endif // !_REGION_
