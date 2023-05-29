#include "rbm/material.hpp"
#include "rbm/meshElement.hpp"
#include "unit_test_framework.hpp"
#include <utility>

TEST(test_MeshElement_1)
{
  // Private variables
  material::Material mat("fuel", 1, 2, 3);
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
TEST(test_checkIndexing)
{

  // Private variables
  material::Material mat("fuel", 3, 6, 9);
  double lx = 4;
  double ly = 5;
  size_t id = 0;
  std::pair<size_t, size_t> idx_row = std::make_pair(7, 8);
  std::pair<size_t, size_t> idx_col = std::make_pair(9, 10);

  // Initialize MeshElement
  mesh::MeshElement element(mat, lx, ly, id, idx_row, idx_col);

  // checking negative indexing for row
  std::pair row_idx = std::make_pair(-1, 1);
  std::pair col_idx = std::make_pair(1, 1);

  ASSERT_FALSE(element.checkIndexing(row_idx, col_idx));

  // checking negative indexing for columns
  row_idx = std::make_pair(1, 1);
  col_idx = std::make_pair(-1, 1);
  ASSERT_FALSE(element.checkIndexing(row_idx, col_idx));

  // checking if first > last for row index
  row_idx = std::make_pair(2, 1);
  ASSERT_FALSE(element.checkIndexing(row_idx, col_idx));

  // checking if first > last for col indexow_idx, col_idx) == false);

  row_idx = std::make_pair(1, 1);
  col_idx = std::make_pair(2, 1);
  ASSERT_FALSE(element.checkIndexing(row_idx, col_idx));
}

TEST_MAIN();
