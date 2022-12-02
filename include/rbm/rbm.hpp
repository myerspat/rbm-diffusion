#ifndef _RBM_
#define _RBM_

#include "rbm/mesh.hpp"
#include "rbm/rbmEnums.hpp"

namespace rbm {
/// Perturb class is brocken down into a work flow,\n
/// 1. Takes in training points\n
///     - These training points can be that of the material class -> absorption,
///     diffusion, nu_fission\n
///     - Functions used, Perturb()\n
/// 2. Calculates Fluxes and criticality factor for each each training point\n
///     - This is done by creating a mesh and constructing the fission and
///     migration matricies\n
///     - Each iteration the material infromation is changed for the new flux
///     and criticality factor\n
/// 3. Does a Principle Component Analysis, PCA, over all the fluxes
/// calculated\n
///     - This uses a singular value decomposition, SVA, algorithm to calculate
///     PCA\n
///     - The default principle components return is 3, but the user can set the
///     amount desired\n
/// 4. Constructs the target Fission Matrix, F^t_ij = <flux_i
/// |F(alpha_t)|flux_j>, and the target Migration Matrix ,M^t_ij =
/// <flux_i|M(alpha_t)|flux_j>, based of the training fluxes.\n
/// 5. Solve for the target criticality and target fluxes.\n
///
class Perturb {
public:
  //=============================================================
  // Constructors / Destructor
  Perturb() {};
  /// Constructor is given all the information in order for a perturbation to
  /// happen and then finding the
  ///
  ///@param training_points is a vector of information that allows a mesh's
  /// material information in a grid to be changed. Training points can be
  /// absorption cross section, fission cross section, and D. Both target points
  /// and training points are of the same type.
  /// @param target_points is a vectory of information that allows the class to
  /// solve for target points specified by the user. Target points can be
  /// absorption cross section, fission cross section, and D. Both target points
  /// and training points are of the same type.
  ///@param mesh is a Mesh object which constructs and holds the grid in which
  /// the finite different will take place.
  ///@param element_id is the id given to the mesh object in order to know which
  /// MeshElement's material in the mesh need to be changed.
  ///@param target_parameter is the paramter in which to be changed. This
  /// specified what material information the training_points are (nu_fission,
  /// absorption, D).
  ///
  Perturb(xt::xarray<double>& training_points,
    xt::xarray<double>& target_points, mesh::Mesh& mesh, size_t& element_id,
    Parameter& target_parameter)
    : _training_points(training_points), _target_points(target_points),
      _mesh(mesh), _element_id(element_id),
      _target_parameter(target_parameter) {};

  //=============================================================
  // Methods

  // Initialize rbm with the training points
  /// Initialized Perturb class with training points
  ///
  /// Training points are the point tp purturb to create a subspace in which a
  /// model can be based off of.
  ///@param training_points is a vector of training points (nu_fission,
  /// absorption, D) of the same type.
  ///@param mesh is the information that the finite difference will be taken
  /// plac on and hold the model information.
  ///@param element_id is the id associated the the mesh's MeshElements. The
  /// MeshElement's ids that is equivilant it element_id will be given the
  /// traning points.
  void initialize(
    xt::xarray<double>& training_points, mesh::Mesh& mesh, size_t& element_id);

  // Create the training subspace
  /// train() creates a default of three different training subspace (traning
  /// points, training fluxes, training k) unless specified otherwise by the
  /// user.
  ///
  /// The training fluxes fluxes and training k arrays need to be calculated
  /// by,\n
  ///   1. Constructing a new mesh for each training point\n
  ///   2. Constructing a new fission matrix and migration matrix for each
  ///   training point\n
  ///   3. Finding the eigen values and iegen vectors of the M and F matricies\n
  ///   4. Storing the information in _training_fluxes, _training_k, and
  ///   _training_points respectively.\n
  void train();

  // Calculate the eigenvalue and eigenvector for the target value
  /// calTarget() is given a target value in which to predict criticality and
  /// fluxes on. This is used after the training data structures are populated
  /// with relavent data
  ///   1. Constructing a new mesh for each training point\n
  ///   2. Constructing a new fission matrix and migration matrix for each
  ///   training point\n
  ///   3. Finding the eigen values and iegen vectors of the M and F matricies\n
  ///   4. returns the information of predicted fluxes and predicted
  ///   criticality\n
  ///@returns nothing, but the the predicted flux array and predicted
  /// criticality value in stored in _target_fluxes  and _target_k respectively.
  void calcTargets();

  // Reduce the subspace using principle component analysis
  /// pcaReduce() does a principle component analysis, PCA, on the all training
  /// fluxes and training criticality and creates a smaller subset of training
  /// fluxes and training criticality based off the weight of the eigen vectors
  /// found. The default amount of fluxes/k values to return is 3, but can be
  /// changed by the user. This is done by using an SVA algorithm.
  /// @param training_fluxes is a 1-D array of fluxes (each training flux is
  /// unique row wise).
  /// @returns nothing but filles _training_points, _training_k,
  /// _training_fluxes with information found.
  void pcaReduce(xt::xarray<double>& training_fluxes);

  std::pair<double, xt::xarray<double>> findMaxEigen(
    xt::xarray<double> eigenvalues, xt::xarray<double> eigenvectors);

  // Check errors
  /// Calculates the percent error for both criticality, k, and fluxes of the
  /// target data and training data.
  void checkError();

  // Construct the target fission matrix
  ///
  /// @params training_fluxes is the 1-2D array of training fluxes
  /// @returns target fission matrix F^t_ij = <flux_i |F(alpha_t)|flux_j>
  xt::xarray<double> constructF_t(
    const xt::xarray<double>& F, const xt::xarray<double>& training_fluxes);

  /// Constructs the target migration matrix
  ///
  /// @params training fluxes is the 1-2D array of training fluxes
  /// @returns target migration matrix M^t_ij = <flux_i|M(alpha_t)|flux_j>
  xt::xarray<double> constructM_t(
    const xt::xarray<double>& M, const xt::xarray<double>& training_fluxes);

  /// Sets the number of Principle Components to keep in reduction
  ///
  ///@params num_pcs number of principle components kept after redution
  ///@returns nothing, but sets _numpcs to num_pcs
  void setNumPCs(const size_t num_pcs) { _num_pcs = num_pcs; }

  // Write data to file_name
  // Writes training data, target data, and error data to files
  // training.csv, target.csv, and error,csv respectively.
  // @param file_name
  void writePointData(const std::string& file_name,
    const xt::xarray<double>& points, const xt::xarray<double>& k,
    const xt::xarray<double>& fluxes);

  // Write PCA data to reduced.csv
  /// Writes reduction data from PCA to reduce.csv.
  /// Also calles writePointData to write the data to their respetive files.
  void writePCAData();

  //=============================================================
  // Getters
  /// Getter function to get the number of training points input by the user.
  /// @returns total number of training points
  const size_t getNumTraining() const { return _training_points.size(); };

  /// Getter function to get the number of target points input by the user
  ///@returns total number of target points
  const size_t getNumTarget() const { return _target_points.size(); };

private:
  //=============================================================
  // Data
  xt::xarray<double> _training_points; // 1D array of training paramater values
  xt::xarray<double> _training_fluxes; // 2D array of training fluxes
  xt::xarray<double> _training_k;      // 1D array of training k
  xt::xarray<double> _target_points;   // 1D array of target parameter values
  xt::xarray<double> _target_fluxes;   // 2D array of target fluxes
  xt::xarray<double> _target_k;        // 2D array of target k
  xt::xarray<double> _variance;        // Variance of PCA data
  mesh::Mesh _mesh;                    // Mesh for the problem
  size_t _num_pcs = 3; // number of PCs to keep when reducing, defaults to 3
  size_t _element_id;  // Element within mesh that will be perturbed
  Parameter _target_parameter; // Perturbed parameter type (absorption, D,
                               // nu_fission)
};

} // namespace rbm

#endif // !_RBM
