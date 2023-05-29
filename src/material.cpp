#include "rbm/material.hpp"

namespace material {

void Material::initialize(
  std::string name, double absorption, double nu_fission, double D)
{}

void Material::setParameter(
  const double& new_value, const Property& target_property)
{
  switch (target_property) {

  case Property::absorption:
    _absorption = new_value;
    break;

  case Property::D:
    _D = new_value;
    break;

  case Property::nu_fission:
    _nu_fission = new_value;
    break;
  }
}

} // namespace material
