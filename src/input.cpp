#include "rbm/input.hpp"
#include "rbm/material.hpp"
#include <algorithm>
#include <cstdio>
#include <pugixml.hpp>
#include <stdexcept>
#include <utility>
#include <xtensor/xadapt.hpp>

namespace util {

void loadFile(pugi::xml_document& file, const std::string& file_name)
{
  // Load file
  auto load = file.load_file(file_name.c_str());

  // Check result
  if (!load) {
    throw std::runtime_error("File " + file_name + " failed to load");
  }
}

pugi::xml_node getNode(
  const pugi::xml_node& node, const std::string& child_name)
{
  // Get child node of name child_name
  const auto child = node.child(child_name.c_str());

  // Check if found
  if (!child) {
    throw std::runtime_error(
      "Child " + child_name + " not found in node " + node.name());
  }
  return child;
}

Settings parseSettingsNode(const pugi::xml_node& root)
{
  Settings settings;

  if (auto settings_node = root.child("settings")) {
    if (auto training_path = settings_node.attribute("training")) {
      // Parse training path
      settings.training_path = training_path.as_string();

    } else if (auto pcs_path = settings_node.attribute("pcs")) {
      // Parse pcs path
      settings.pcs_path = pcs_path.as_string();
    }

    if (auto calc_errors = settings_node.attribute("calc_errors")) {
      // Parse calculate errors bool
      settings.calc_errors = calc_errors.as_bool();
    }
  }

  return settings;
}

std::vector<material::Material> parseMaterialsNode(const pugi::xml_node& root)
{
  // Get materials node
  pugi::xml_node materials_node = util::getNode(root, "materials");
  std::vector<material::Material> materials;

  // Iterate through all material nodes and fill materials vector
  for (pugi::xml_node material_node = util::getNode(materials_node, "material");
       material_node; material_node = material_node.next_sibling("material")) {
    materials.emplace_back(
      util::getAttribute<std::string>(material_node, "name"),
      util::getAttribute<double>(material_node, "absorption"),
      util::getAttribute<double>(material_node, "nu_fission"),
      util::getAttribute<double>(material_node, "D"));
  }

  return materials;
}

mesh::Mesh parseMeshNode(const pugi::xml_node& root)
{
  // Notify user of parsing
  std::cout << "\n   Parsing mesh XML node\n";

  // Get material node
  std::vector<material::Material> materials = util::parseMaterialsNode(root);

  // Get mesh node from root
  pugi::xml_node mesh_node = util::getNode(root, "mesh");

  // Get cartesian_mesh node
  pugi::xml_node cartesian_mesh_node =
    util::getNode(mesh_node, "cartesian_mesh");

  // Get total x and y bins in the mesh
  size_t fine_x_bins =
    util::getAttribute<size_t>(cartesian_mesh_node, "fine_x_bins");
  size_t fine_y_bins =
    util::getAttribute<size_t>(cartesian_mesh_node, "fine_y_bins");
  size_t course_x_bins =
    util::getAttribute<size_t>(cartesian_mesh_node, "course_x_bins");
  size_t course_y_bins =
    util::getAttribute<size_t>(cartesian_mesh_node, "course_y_bins");

  std::cout << "     Total x bins = " << fine_x_bins * course_x_bins
            << ", Total y bins = " << fine_y_bins * course_y_bins << std::endl;

  // Get boundary conditions
  std::vector<std::string> sides = {"left", "right", "top", "bottom"};
  std::vector<std::pair<double, double>> bounds;
  bounds.reserve(4);

  // For each side get pair of a and b for boundary conditions
  for (const auto& side : sides) {
    pugi::xml_node bound_node = util::getNode(mesh_node, side + "_condition");
    bounds.emplace_back(util::getAttribute<double>(bound_node, "a"),
      util::getAttribute<double>(bound_node, "b"));

    // Print bound
    std::cout << "     " + side + " bound: a = " << bounds.back().first
              << ", b = " << bounds.back().second << std::endl;
  }

  // Initialize mesh
  mesh::Mesh mesh(fine_x_bins, fine_y_bins, course_x_bins, course_y_bins,
    bounds[0], bounds[1], bounds[2], bounds[3]);

  // Parse cartesian_mesh node for mesh elements
  std::vector<mesh::MeshElement> elements;

  // Reserve the maximum allowed elements
  elements.reserve(course_x_bins * course_y_bins);

  // Iterate through elements
  for (pugi::xml_node element_node =
         util::getNode(cartesian_mesh_node, "element");
       element_node; element_node = element_node.next_sibling("element")) {
    // Get attributes: id, x_length, y_length, and row/col positions
    size_t id = util::getAttribute<size_t>(element_node, "id");
    double x_length = util::getAttribute<double>(element_node, "x_length");
    double y_length = util::getAttribute<double>(element_node, "y_length");
    size_t start_row_idx =
      util::getAttribute<size_t>(element_node, "start_row_idx");
    size_t stop_row_idx =
      util::getAttribute<size_t>(element_node, "stop_row_idx");
    size_t start_col_idx =
      util::getAttribute<size_t>(element_node, "start_col_idx");
    size_t stop_col_idx =
      util::getAttribute<size_t>(element_node, "stop_col_idx");

    // Get material attribute and find the associated material in the vector
    std::string material_str =
      util::getAttribute<std::string>(element_node, "material");
    auto material_itr = std::find_if(materials.cbegin(), materials.cend(),
      [&](const material::Material& a) { return a.getName() == material_str; });

    // Check that the material was found in the vector
    if (material_itr == materials.cend()) {
      throw std::runtime_error("The material " + material_str +
                               " was not given in the materials node");
    }

    // Construct MeshElement in vector
    elements.emplace_back(*material_itr, x_length, y_length, id,
      std::make_pair(start_row_idx, stop_row_idx),
      std::make_pair(start_col_idx, stop_col_idx));
  }
  elements.shrink_to_fit();

  std::cout << "     Number of course mesh elements = " << elements.size()
            << std::endl;

  // Construct mesh with mesh elements
  mesh.constructMesh(elements);

  return mesh;
}

rbm::Perturb parseRBMNode(const pugi::xml_node& root, mesh::Mesh& mesh)
{
  // Notify parsing of rbm node
  std::cout << "\n   Parsing rbm XML node\n";

  // Get rbm node
  pugi::xml_node rbm_node = util::getNode(root, "rbm");

  // Get parameters for multi-parameter perturbation
  std::vector<rbm::Parameter> parameters;
  bool constructD = false;
  for (pugi::xml_node parameter_node = util::getNode(rbm_node, "parameter");
       parameter_node;
       parameter_node = parameter_node.next_sibling("parameter")) {
    // Parse region/meshElement ID and property
    size_t id = util::getAttribute<std::size_t>(parameter_node, "element_id");
    std::string property_str =
      util::getAttribute<std::string>(parameter_node, "property");

    // Convert string to enum
    material::Property property;
    if (property_str == "absorption") {
      property = material::Property::absorption;
    } else if (property_str == "nu_fission") {
      property = material::Property::nu_fission;
    } else if (property_str == "D") {
      property = material::Property::D;
      constructD = true;
    } else {
      throw std::runtime_error("The property provided (" + property_str +
                               ") is not absorption, nu_fission, or D");
    }

    // Parse training and target points from parameter node
    xt::xarray<double> training_points =
      xt::adapt(util::parseString<double>(parameter_node, "training"));
    xt::xarray<double> target_points =
      xt::adapt(util::parseString<double>(parameter_node, "target"));

    parameters.emplace_back(id, property, training_points, target_points);

    // Print target parameter and ID of elements to be perturbed
    std::cout << "     Perturbing elements of ID = " << id << std::endl;
    std::cout << "     Target property: " + property_str << std::endl;
  }

  // Assert the perturbation arrays are the same size
  size_t training_length = parameters[0].getTrainingPoints().size();
  size_t target_length = parameters[0].getTargetPoints().size();
  for (size_t i = 1; i < parameters.size(); i++) {
    assert(training_length == parameters[i].getTrainingPoints().size());
    assert(target_length == parameters[i].getTargetPoints().size());
  }

  // Construct perturb object
  rbm::Perturb perturb(parameters, mesh, constructD);

  // Determine number of PCAs to preserve, default is 3
  size_t num_pcs = 3;
  if (auto attr = rbm_node.attribute("pcs")) {
    num_pcs = attr.as_int();
    perturb.setNumPCs(num_pcs);

  } else if (training_length < 3) {
    num_pcs = training_length;
    perturb.setNumPCs(num_pcs);
  }

  // Determine if EIM will be used default to true (only used on diffusion
  // perturbations)
  if (auto attr = rbm_node.attribute("EIM")) {
    perturb.setEIM(attr.as_bool());
  }

  // Precompute with affine decomposition
  if (auto attr = rbm_node.attribute("precompute")) {
    perturb.setPrecompute(attr.as_bool());
  }

  // Read FOPT attribute if applicable
  if (auto attr = rbm_node.attribute("FOPT")) {
    perturb.setFOPT(attr.as_bool());
  }

  // Print size of subspace
  std::cout << "     Number of PCs = " << num_pcs << std::endl;

  // Print number of target points
  std::cout << "     Number of target points = " << target_length << std::endl;

  return perturb;
}

} // namespace util
