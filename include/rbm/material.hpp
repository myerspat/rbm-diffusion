#ifndef _MATERIAL_
#define _MATERIAL_

#include "rbm/rbmEnums.hpp"
#include <string>
/// Contains all the properties of a material for a diffusion problem
///
/// This includes the diffusion coefficient, D, Fission Cross Section,
/// nu_fission, absorption cross section, absorption, and  the name of the
/// material, name.
class Material {
public:
  //========================================================================
  // Constructor / Destructor

  Material() {};
  // Constructor of the Material class.
  ///
  /// Initialized the material needed for a diffusion problem
  ///
  /// @param name is the name of the structure being created
  ///             Example: Fuel or Reflector
  /// @param absorption is the absorption macroscopic cross section for a given
  /// structure
  /// @param nu_fission is the fission macroscopic cross sections for a given
  /// structure
  /// @param D is the diffusion coefficient
  /// @returns Reinitializes the private variables with their respective
  /// information above
  Material(const std::string& name, const double& absorption,
    const double& nu_fission, const double& D)
    : _name(name), _absorption(absorption), _nu_fission(nu_fission), _D(D) {};

  //========================================================================
  // Methods

  /// Re-Initialized the material needed for a diffusion problem
  ///
  /// @param name is the name of the structure being created
  ///             Example: Fuel or Reflector
  /// @param absorption is the absorption macroscopic cross section for a given
  /// structure
  /// @param nu_fission is the fission macroscopic cross sections for a given
  /// structure
  /// @param D is the diffusion coefficient
  void initialize(
    std::string name, double absorption, double nu_fission, double D);

  //========================================================================
  // Getters

  /// Getter function to retrieve the private variable _name. Which is the name
  /// of the function
  ///
  /// @return _name
  const std::string& getName() const { return _name; };

  /// Getter function to retrieve the private variable _absorption. This is the
  /// absorption macroscopic cross section of a material.
  ///
  /// @return _absorption
  const double& getAbsorption() const { return _absorption; };

  /// Getter function to retrieve the private variable _nu_fission. This is the
  /// vfission macroscopic cross section.
  ///
  /// @return _nu_fission
  const double& getNuFission() const { return _nu_fission; };

  /// Getter function to retrieve the private variable _D. This is the diffusion
  /// coefficient related the the material.
  ///
  /// @return _D
  const double& getD() const { return _D; };

  //========================================================================
  // Setters
  // Set the value of a target_parameter

  /// Sets a new value for a given private variable.
  ///
  /// @param new_value is the new value or replacement value for a given private
  /// variable.
  /// @param target_parameter is the target parameter or key word to change the
  /// private variable.
  void setParameter(
    const double& new_value, const rbm::Parameter& target_parameter);

private:
  //========================================================================
  // Data
  std::string _name;
  double _absorption;
  double _nu_fission;
  double _D;
};

#endif // !_MATERIAL_
