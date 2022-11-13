#ifndef _MATERIAL_
#define _MATERIAL_

#include <string>

class Material {
public:
  //========================================================================
  // Constructor / Destructor
  Material() {};
  Material(std::string name, double absorption, double nu_fission, double D)
    : _name(name), _absorption(absorption), _nu_fission(nu_fission), _D(D) {};

  //========================================================================
  // Methods
  void initialize(
    std::string name, double absorption, double nu_fission, double D);

  //========================================================================
  // Getters
  std::string getName() { return _name; };
  double getAbsorption() { return _absorption; };
  double getNuFission() { return _nu_fission; };
  double getD() { return _D; };

  //========================================================================
  // Setters
  void setAbsorption(double value) { _absorption = value; };
  void setNuFission(double value) { _nu_fission = value; };
  void setD(double value) { _D = value; };

private:
  //========================================================================
  // Data
  std::string _name;
  double _absorption;
  double _nu_fission;
  double _D;
};

#endif // !_MATERIAL_
