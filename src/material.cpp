#include "rbm/material.hpp"
#include "rbm/rbmEnums.hpp"

void Material::initialize(
  std::string name, double absorption, double nu_fission, double D)
{}

void Material::setParameter(
  const double& new_value, const rbm::Parameter& target_parameter)
{
  
  switch(target_parameter){

    case rbm::absorption : 
      _absorption = new_value;

    case rbm::D :
      _D = new_value;

    case rbm::nu_fission :
      _nu_fission = new_value;
  }
}
