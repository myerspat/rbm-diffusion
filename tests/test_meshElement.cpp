#include "rbm/material.hpp"
#include "rbm/meshElement.hpp"
#include "unit_test_framework.hpp"
#include <utility>

TEST(test_MeshElement_1)
{
  // Private variables
  Material mat("fuel", 1, 2, 3);
  double lx = 4;
  double ly = 5;
  size_t id = 6;
  std::pair<size_t, size_t> idx_row = std::make_pair(7, 8);
  std::pair<size_t, size_t> idx_col = std::make_pair(9, 10);

  // Initialize MeshElement
  mesh::MeshElement element(mat, lx, ly, id, idx_row, idx_col);

  // Assertions
  ASSERT_EQUAL(mat.getName(), element.getMaterial().getName());
  ASSERT_EQUAL(lx, element.getLX());
  ASSERT_EQUAL(ly, element.getLY());
  ASSERT_EQUAL(id, element.getID());
  ASSERT_EQUAL(idx_row.second, element.getRowIdx().second);
  ASSERT_EQUAL(idx_col.first, element.getColIdx().first);
}

TEST_MAIN();
