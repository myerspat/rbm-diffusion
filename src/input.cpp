#include "rbm/input.hpp"
#include "rbm/material.hpp"
#include <algorithm>
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

std::vector<Material> parseMaterialsNode(const pugi::xml_node& root)
{
  // Get materials node
  pugi::xml_node materials_node = util::getNode(root, "materials");
  std::vector<Material> materials;

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
  // Get material node
  std::vector<Material> materials = util::parseMaterialsNode(root);

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

  // Get boundary conditions
  std::vector<std::string> sides = {"left", "right", "top", "bottom"};
  std::vector<std::pair<double, double>> bounds;
  bounds.reserve(4);

  // For each side get pair of a and b for boundary conditions
  for (const auto& side : sides) {
    pugi::xml_node bound_node = util::getNode(mesh_node, side + "_condition");
    bounds.emplace_back(util::getAttribute<double>(bound_node, "a"),
      util::getAttribute<double>(bound_node, "b"));
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
      [&](const Material& a) { return a.getName() == material_str; });

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

  // Construct mesh with mesh elements
  mesh.constructMesh(elements);

  return mesh;
}

rbm::Perturb parseRBMNode(const pugi::xml_node& root, mesh::Mesh& mesh)
{
  // Get rbm node
  pugi::xml_node rbm_node = util::getNode(root, "rbm");

  // Get training node
  pugi::xml_node training_node = util::getNode(rbm_node, "training");

  // Get cell id that will be perturbed
  size_t element_id = util::getAttribute<int>(training_node, "element_id");

  // Get target parameter
  std::string parameter =
    util::getAttribute<std::string>(training_node, "parameter");

  // Convert target parameter to enum value
  rbm::Parameter target_parameter;
  if (parameter == "absorption") {
    target_parameter = rbm::Parameter::absorption;
  } else if (parameter == "nu_fission") {
    target_parameter = rbm::Parameter::nu_fission;
  } else if (parameter == "D") {
    target_parameter = rbm::Parameter::D;
  } else {
    throw std::runtime_error("The parameter provided (" + parameter +
                             ") is not absorption, nu_fission, or D");
  }

  // Parse training fluxes
  xt::xarray<double> training_points =
    xt::adapt(util::parseString<double>(training_node, "values"));

  // Create Purturb object
  rbm::Perturb perturb(training_points, mesh, element_id, target_parameter);

  // Determine number of PCAs to preserve, default is 3
  if (auto attr = training_node.attribute("pcas")) {
    perturb.setNumPCs(attr.as_int());
  } else if (training_points.shape(0) < 3) {
    // If the number of training points given is less than the default (3)
    // reduce number of PCs kept
    perturb.setNumPCs(training_points.shape(0));
  }

  return perturb;
}

} // namespace util
