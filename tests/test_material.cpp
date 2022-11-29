#include "rbm/material.hpp"
#include "unit_test_framework.hpp"

TEST(test_Material_1) {
  // Initialize material
  Material mat("fuel", 1, 2, 3);

  // Assertions
  ASSERT_EQUAL("fuel", mat.getName());
  ASSERT_EQUAL(1, mat.getAbsorption());
  ASSERT_EQUAL(2, mat.getNuFission());
  ASSERT_EQUAL(3, mat.getD());
}

TEST_MAIN();
