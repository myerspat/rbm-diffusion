#include "rbm/material.hpp"
#include "unit_test_framework.hpp"

TEST(test_Material_1)
{
  // Initialize material
  Material mat("fuel", 1, 2, 3);

  // Assertions
  ASSERT_EQUAL("fuel", mat.getName());
  ASSERT_EQUAL(1, mat.getAbsorption());
  ASSERT_EQUAL(2, mat.getNuFission());
  ASSERT_EQUAL(3, mat.getD());
}
TEST(test_setParameters)
{
  //Contructing Material 
  Material mat("fuel", 1, 2, 3);
  // Setting new value for absorption
  double new_value = 10;
  rbm::Parameter target_type = rbm::Parameter::absorption;
  mat.setParameter(new_value, target_type);
  // checking if value changed is correct
  ASSERT_EQUAL(mat.getAbsorption(), 10);

  // checking new value for nu_fission
  target_type = rbm::Parameter::nu_fission;
  new_value = 11;
  mat.setParameter(new_value, target_type);
  //checking if nu fission cahnged
  ASSERT_EQUAL(mat.getNuFission(), 11);

  // Checking new value for D
  target_type = rbm::Parameter::D;
  new_value = 12;
  mat.setParameter(new_value, target_type);
  ASSERT_EQUAL(mat.getD(), 12);

}

TEST_MAIN();
