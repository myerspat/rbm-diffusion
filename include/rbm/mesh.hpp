#ifndef _MESH_
#define _MESH_

#include "rbm/meshElement.hpp"
#include "rbm/rbmEnums.hpp"
#include <cassert>
#include <utility>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>

namespace mesh {

/// Mesh class constructs a mesh based off of MeshElements provided in a vector,
/// 4 boundary conditions and xN's/yN's for both a course and fine grid.
///
/// The mesh class takes input for boundary conditions for a 2D diffusion
/// problem(This problem can also be 1D, but a new constructor is used) along
/// with a vector of MeshElements and iterative chunks to make a grid.
///
/// The MeshElement vector should be filled with MeshElements to make a 2D grid,
/// M by N, in order to make a course grid. A course grid is the physical
/// geometry of the problem at hand. It contains the distance of each N and M
/// chunk in the lx/ly direction along with its material properties.\n General
/// Structure of course Grid: \n 1. For each row of the course grid the
/// MeshElements ly must be same.\n 2. For each column of the course grid the
/// MeshElements lx must be the same.\n 3. For a N x M grid you can have less
/// MeshElements than position as long as row_idx and col_idx cover the rest of
/// the course grid.
///
/// Once the course grid is made then the fine grid can be made. Depending on
/// the users input a single MeshElement is brocken down into a N x M grid. Each
/// MeshElements fine grid is still connected together by their respective (top,
/// bottom, left, right) MeshElements orientation in the course grid.
///
/// Once the fine grid is made the Migration and Fission matricies can be solved
/// for.
class Mesh {
private:
  //========================================================================
  // Data
  xt::xarray<MeshElement> _fine_grid;      // Mesh centered elements
  xt::xarray<MeshElement> _course_grid;    // Mesh centered elements
  std::pair<double, double> _left_bound;   // Left boundary condition
  std::pair<double, double> _right_bound;  // Right boundary condition
  std::pair<double, double> _top_bound;    // Top boundary condition
  std::pair<double, double> _bottom_bound; // Bottom boundary condition
  size_t _xN_fine;
  size_t _yN_fine;
  size_t _xN_course; // Number of mesh elements in the x direction
  size_t _yN_course; // Number of mesh elements in the y direction

public:
  //========================================================================
  // Constructors
  Mesh() {};
  /// Constructor for the mesh
  ///
  /// @param xN_fine is number of bins to break up 1 MeshElement in the x
  /// direction.
  /// @param yN_fine is the number of bins to break up 1 MeshElement in the y
  /// direction.
  /// @param xN_course is the number of bins to break up all MeshElements in the
  /// x direction.
  /// @param yN_course is the number of bins to break up all MeshElements in the
  /// y direction.
  /// @param left_bound is the boundary condition (a,b) for a diffusion problem
  /// on the left surface of a square geometry.
  /// @param right_bound is the boundary condition (a,b) for a diffusion problem
  /// on the right surface of a square geometry.
  /// @param top_bound is the boundary condition (a,b) for a diffusion problem
  /// on the top surface of a square geometry.
  /// @param bottom_bound is the boundary conditions (a,b) for a diffusion
  /// problem on the bottom surface a square geometry.
  ///@warning For ALL boundary conditions a and b cannot both equal 0.
  Mesh(const size_t& xN_fine, const size_t& yN_fine, const size_t& xN_course,
    const size_t& yN_course, const std::pair<double, double>& left_bound,
    const std::pair<double, double>& right_bound,
    const std::pair<double, double>& top_bound,
    const std::pair<double, double>& bottom_bound)
    : _xN_fine(xN_fine), _yN_fine(yN_fine), _xN_course(xN_course),
      _yN_course(yN_course), _left_bound(left_bound), _right_bound(right_bound),
      _top_bound(top_bound), _bottom_bound(bottom_bound)
  {
    // Asserting a and b are both not equal to zero for all boundary conditions
    assert(!(left_bound.first == 0.0 && left_bound.second == 0.0));
    assert(!(right_bound.first == 0.0 && right_bound.second == 0.0));
    assert(!(top_bound.first == 0.0 && top_bound.second == 0.0));
    assert(!(bottom_bound.first == 0.0 && bottom_bound.second == 0.0));
  }

  //========================================================================
  // Methods

  // Construct mesh
  /// constructMesh constructs both the fine and course grids when given a
  /// vector of MeshElements.
  /// @param elements in a vector of mesh elements that makes a N x M course
  /// grid.
  /// @warning if elements is not a full NXM grid and a position is left out
  /// then the mesh will not work properly.
  void constructMesh(const std::vector<MeshElement>& elements);

  // Construct course grid
  /// constructCourseGrid constructs a ( N x M ) course grid of MeshElements
  /// given a vector of mesh elements. In our case the course grid would be
  /// (xN_course by xY_course) in dimensions.
  ///@param elements is a vector of MeshElements that make up a course grid.
  ///@returns void type but private variable _course_grid get modified.
  ///  @warning if elements is not a full NXM grid and a position is left out
  ///  then the mesh will not work properly.
  xt::xarray<MeshElement> constructCourseGrid(
    const std::vector<MeshElement>& elements);

  // Construct fine grid
  /// constructFineGrid construct a finer grid relative to the course grid of
  /// mesh elements. In our case the fine grid would be a total width of
  /// (xN_course * xN_fine by yN_course * yN_fine) of MeshElements.
  /// @returns void type, but _fine_grid get modified;
  xt::xarray<MeshElement> constructFineGrid(
    const xt::xarray<MeshElement>& course_grid);

  // Check lengths of adjacent elements match for all elements in _course_grid
  /// checkSharedLengths makes sure user input is valid for a square/rectangular
  /// geometry. Each MeshElement in a row must have the same ly and each mesh
  /// element in a column must have the same lx for a course grid to be valid.
  ///@param course_grid takes in a course grid that is constructed to check if
  /// it is valid.
  ///@returns true of rows/columns meet the lengths requirements statement else
  /// false.
  ///@warning If a course grid is not valid then the geometry is not a
  /// square/rectangule and the program exits.
  bool checkSharedLengths(const xt::xarray<MeshElement>& course_grid);

  // Change target mesh elements' material property
  /// Changes the material properties of the MeshElement for a course and fine
  /// grid in respect to the MeshElements ID. If multiple MeshElements have the
  /// same id then each ones material properties will be changed respectively.
  /// @param id is the unique tag given to each mesh element.
  /// @param new_value is the new value to change the materal properties of the
  /// MeshElement to.
  /// @param target_parameter is the material properties key word to
  /// know which material property to change. For example if target_parameter=D
  /// then change D.
  /// @returns nothing, but changed _course_grid and _fine_grid.
  /// @warning if given a non-valid id or target_parameter value nothing should
  /// change in the grids.
  void changeMaterial(const std::size_t& id, const double& new_value,
    const rbm::Parameter& target_parameter);

  /// Constructs the fission operator matrix
  ///
  /// Iterates through fine grid and constructs the fission operator matrix.
  /// @returns xt::xarray <doubles>.
  xt::xarray<double> constructF(); // Construct the fission operator matrix

  // Constructs the Migration Operator matrix
  //
  // Iterates the fine grid and constructs the migration operator matrix.
  // @returns xt::xarray<doubles> Migration Matrix
  xt::xarray<double> constructM(); // Construct the migration operator matrix

  // Assuming row major ordering, returns 1D idx given 2D idx
  size_t ravelIDX(const size_t& i, const size_t& j);

  //========================================================================
  // Getters
  // Get the number of bins along the x-axis
  /// Getter function to get the fine grid's number of elements in the x
  /// direction.
  /// @returns size_t number of positions in the x direction of fine mesh.
  size_t getXN() const { return _xN_fine * _xN_course; };

  // Gets the number of bins along the y-axis.
  /// Getter function to get the fine grid's number of elements in the y
  /// direction.
  /// @returns size_t is the number of position in the y direction in the y
  /// direction.
  size_t getYN() const { return _yN_fine * _yN_course; };

  // Get the number of total bins in the mesh
  // Getter Function to get the size of the fine grid ie the total number of
  // elements.
  // @returns size_t if total number of positions in the fine mesh
  size_t getSize() const { return getXN() * getYN(); };

  // Gets the course grid.
  /// Getter function allowing you to return the private course_grid.
  ///@returns xt::xarray<MeshElement> of a course grid
  xt::xarray<MeshElement> getCourseGrid() const { return _course_grid; }

  // Gets the fine grid
  /// Getter function allowing you to return the private fine_grid
  ///@returns xt::xarray<MeshElement> of a fine grid
  xt::xarray<MeshElement> getFineGrid() const { return _fine_grid; }
};

} // namespace mesh

#endif // !_MESH_
