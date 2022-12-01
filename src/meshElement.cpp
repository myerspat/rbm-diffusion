#include "rbm/meshElement.hpp"

namespace mesh {


bool MeshElement::checkIndexing(const std::pair<std::size_t, std::size_t>& idx_row,
                                const std::pair<std::size_t, std::size_t>& idx_col)
{
  // checking if first < second
   if(_idx_row.first > idx_row.second ||  idx_col.first > idx_col.second){
     return false;
   }
  // checking if both values are non-negative
  if(idx_row.first < 0 || idx_row.second < 0 || idx_col.first < 0 ||
     idx_col.second < 0){
    return false;
  }
  // else indexing should be good

  return true;
}

} // namespace mesh
