#ifndef _MATERIAL_
#define _MATERIAL_

#include <string>

class Material {
public:
  //========================================================================
  // Constructor / Destructor
  Material() {};
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
  void setAbsorption(double& value) { _absorption = value; };
  void setNuFission(double& value) { _nu_fission = value; };
  void setD(double& value) { _D = value; };

private:
  //========================================================================
  // Data
  std::string _name;
  double _absorption;
  double _nu_fission;
  double _D;
};

#endif // !_MATERIAL_
