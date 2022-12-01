#include "rbm/material.hpp"
#include "rbm/rbmEnums.hpp"

void Material::initialize(
  std::string name, double absorption, double nu_fission, double D)
{}

void Material::setParameter(
  const double& new_value, const rbm::Parameter& target_parameter)
{
  switch (target_parameter) {

  case rbm::Parameter::absorption:
    _absorption = new_value;
    break;

  case rbm::Parameter::D:
    _D = new_value;
    break;

  case rbm::Parameter::nu_fission:
    _nu_fission = new_value;
    break;
  }
}
