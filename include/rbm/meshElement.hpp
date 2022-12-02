#ifndef _MESH_ELEMENT_
#define _MESH_ELEMENT_

#include "rbm/material.hpp"
#include <assert.h>
#include <utility>

namespace mesh {

/// MeshElement class allows you to create a piece of a course mesh.
///
/// MeshElements are put together to form a course mesh. These are the physical
/// structures that make up a model of the diffusion model being represented.
/// Each mesh element is given indexes where the position in the rows and
/// columns of a grid where the MeshElements should be orientated along with
/// its material length in the x direction and length in the y direction.
///\n
/// Example: physcial representation MeshElement \n
/// |------------------------ |\n
/// |             |           |carries.. ID =0,\n
/// |             |           |material -> D, absorption, nu_fission, name,\n
/// |________lx_______________|idx_row = [0, 2] inclusive for left and right \n
/// |             ly          |idx_col = [1, 3] inclusive for left and right \n
/// |             |           |\n
/// |-------------------------|\n
class MeshElement {

public:
  //=========================================================================================
  // Constructors / Destructor
  MeshElement() {};
  /// This is the constructor the the MeshElement class
  ///
  /// The constructor takes a material, lengths in the x direction, length
  /// in the y direction, id, and two pairs of indices for rows and columns.
  /// @param mat is the material information for a given position of a spacial
  /// structure.
  /// @param lx is the length in the x direction of the MeshElement's material
  /// information.
  /// @param ly is the length in y direction of the MeshElement's material
  /// information.
  /// @param id is the id that tags the specific mesh element. (This is used to
  /// eventually change multiple MeshElements information all at once)
  /// @param idx_row is the row indexing information for a specfic mesh element.
  /// For a given pair, example (0, 1), the mesh element will fill the indexing
  /// orientation on a grid. Both left and right bounds are inclusive in a pair.
  /// @param idx_row is the column indexing information for a specific mesh
  /// element. For a given pair the mesh element will fill the indexing
  /// orientation on a grid. Both left and right bound are inclusive in a pair.
  MeshElement(const Material& mat, const double& lx, const double& ly,
    const std::size_t& id, const std::pair<std::size_t, std::size_t>& idx_row,
    const std::pair<std::size_t, std::size_t>& idx_col)
    : _idx_row(idx_row), _idx_col(idx_col), _lx(lx), _ly(ly), _id(id), _mat(mat)
  {
    assert(checkIndexing(idx_row, idx_col));
  };

  //=========================================================================================
  // Methods

  /// checkIndexing is a user error checking function
  ///
  /// This function makes sure idx_col and idx_row are always from least to
  /// greates ie pair.first > pair.second for each respective pair. The indexes
  /// also should be non negative.
  ///
  /// @param idx_row is the row indexing range at which the grid will index in
  /// order to be stored. Both these points will be inclusive and values
  /// assosiated to them will also be stored along with values between them.
  /// @param idx_col is the column indexing range at which the grid will index
  /// in order to be stored. Both these points will be inclusive and the values
  /// associated to them will also be stored along with values between the two
  /// points.
  /// @returns A boolian checking non negative and increase order from first to
  /// second.
  bool checkIndexing(const std::pair<std::size_t, std::size_t>& idx_row,
    const std::pair<std::size_t, std::size_t>& idx_col);
  //=========================================================================================

  //=========================================================================================

  // Getters
  /// Getter function for material.
  ///
  /// @returns _mat
  const Material& getMaterial() const { return _mat; };
  /// Getter function for id.
  ///
  /// @returns _id
  const std::size_t& getID() const { return _id; };

  /// Getter function for row indexing information.
  ///
  ///@returns _idx_row
  const std::pair<std::size_t, std::size_t>& getRowIdx() const
  {
    return _idx_row;
  };
  /// Getter function for column indexing information.
  ///
  /// returns _idx_col
  const std::pair<std::size_t, std::size_t>& getColIdx() const
  {
    return _idx_col;
  };
  /// Getter function for x length information.
  ///
  /// @returns _lx
  const double& getLX() const { return _lx; };

  /// Getter function for y length information
  ///
  ///@returns _ly
  const double& getLY() const { return _ly; };

  //=========================================================================================
  // Setters

  /// Sets a new value for a material variable
  ///
  /// If a material needs to be changed this allows the MeshElements material to
  /// be accessed and changed.
  /// @param new_value the new value in the material to be changed.
  /// @param target_parameter the material information key word to be changed.
  /// #returns nothing, bit changes material information.
  void setParameter(
    const double& new_value, const rbm::Parameter& target_parameter)
  {
    _mat.setParameter(new_value, target_parameter);
  };

private:
  //=========================================================================================
  // Data
  Material _mat;
  std::pair<std::size_t, std::size_t> _idx_row, _idx_col;
  double _lx;
  double _ly;
  std::size_t _id;
};

} // namespace mesh

#endif // !_MESH_ELEMENT_
