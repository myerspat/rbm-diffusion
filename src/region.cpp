#include "rbm/region.hpp"
#include "xtensor/xbuilder.hpp"
#include <xtensor/xio.hpp>

namespace mesh {

Region::Region(const xt::xarray<MeshElement>& fine_grid, const size_t& xN_fine,
  const size_t& yN_fine, const std::size_t& id, const material::Material& mat)
  : _id(id), _mat(mat)
{
  // Initialize masks
  _mask = xt::zeros<double>({fine_grid.size()});

  // Ravel index function to convert from position in space to position in the
  // matrix
  auto ravelIDX = [&](std::size_t i, std::size_t j) {
    return j + i * fine_grid.shape(1);
  };

  for (size_t i = 0; i < fine_grid.shape(0); i++) {
    for (size_t j = 0; j < fine_grid.shape(1); j++) {
      if (fine_grid(i, j).getID() == _id) {
        // Update absorption mask
        _mask(ravelIDX(i, j)) =
          fine_grid(i, j).getLX() / xN_fine * fine_grid(i, j).getLY() / yN_fine;
      }
    }
  }
}

void Region::changeMaterial(
  const double& new_value, const material::Property& target_property)
{
  _mat.setParameter(new_value, target_property);
}

} // namespace mesh
