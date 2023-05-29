#ifndef _PARAMETER_
#define _PARAMETER_

#include "rbm/material.hpp"
#include <xtensor/xarray.hpp>

namespace rbm {

class Parameter {
private:
  //========================================================================
  // Data
  size_t _region_id;
  material::Property _property;
  xt::xarray<double> _training_points;
  xt::xarray<double> _target_points;

public:
  //========================================================================
  // Constructors
  Parameter() {};
  Parameter(const size_t& region_id, const material::Property& property,
    const xt::xarray<double>& training_points,
    const xt::xarray<double>& target_points)
    : _region_id(region_id), _property(property),
      _training_points(training_points), _target_points(target_points) {};

  //========================================================================
  // Getters
  size_t getID() const { return _region_id; }
  material::Property const getMaterialProperty() { return _property; }
  xt::xarray<double> const getTrainingPoints() { return _training_points; }
  xt::xarray<double> const getTargetPoints() { return _target_points; }
};

} // namespace rbm

#endif // !_PARAMETER_
