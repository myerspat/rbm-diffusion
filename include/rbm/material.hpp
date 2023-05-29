#ifndef _MATERIAL_
#define _MATERIAL_

#include <string>

namespace material {

enum Property { absorption, D, nu_fission };

class Material {
public:
  //========================================================================
  // Constructor / Destructor

  Material() {};
  /// Constructor of the Material class.
  Material(const std::string& name, const double& absorption,
    const double& nu_fission, const double& D)
    : _name(name), _absorption(absorption), _nu_fission(nu_fission), _D(D) {};

  //========================================================================
  // Methods
  void initialize(
    std::string name, double absorption, double nu_fission, double D);

  //========================================================================
  // Getters
  const std::string& getName() const { return _name; };
  const double& getAbsorption() const { return _absorption; };
  const double& getNuFission() const { return _nu_fission; };
  const double& getD() const { return _D; };

  //========================================================================
  // Setters
  // Set the value of a target_parameter
  void setParameter(const double& new_value, const Property& target_property);

private:
  //========================================================================
  // Data
  std::string _name;
  double _absorption;
  double _nu_fission;
  double _D;
};

} // namespace material

#endif // !_MATERIAL_
