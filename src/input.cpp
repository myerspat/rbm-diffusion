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

// std::map<int, mesh::Section> parseSectionsNode(const pugi::xml_node& root)
// {
//   // Get sections node
//   pugi::xml_node sections_node = util::getNode(root, "sections");
//   std::map<int, mesh::Section> sections;
//
//   // Iterate through all section nodes and fill sections vector
//   for (pugi::xml_node section_node = util::getNode(sections_node, "section");
//        section_node; section_node = section_node.next_sibling("section")) {
//     int section_id = util::getAttribute<size_t>(section_node, "id");
//     mesh::Section section(section_id,
//       std::make_pair(util::getAttribute<double>(section_node, "x0"),
//         util::getAttribute<double>(section_node, "y0")),
//       util::getAttribute<double>(section_node, "l1"),
//       util::getAttribute<double>(section_node, "l2"));
//
//     sections.emplace(section_id, section);
//   }
//
//   return sections;
// }
//
// std::vector<mesh::Cell> parseCellNode(
//   const pugi::xml_node& root, const std::map<int, mesh::Section>& sections)
// {
//   // Parse materials node
//   std::vector<Material> materials = util::parseMaterialsNode(root);
//
//   // Get cells node
//   pugi::xml_node cells_node = util::getNode(root, "cells");
//   std::vector<mesh::Cell> cells;
//
//   // Iterate through all section nodes and fill sections vector
//   for (pugi::xml_node cell_node = util::getNode(cells_node, "cell"); cell_node;
//        cell_node = cell_node.next_sibling("cell")) {
//     // Get cell id
//     int cell_id = util::getAttribute<int>(cell_node, "id");
//
//     // Get cell material name
//     std::string cell_material_name =
//       util::getAttribute<std::string>(cell_node, "material");
//
//     // Get cell material
//     Material cell_material = *std::find_if(materials.begin(), materials.end(),
//       [&](Material& a) { return a.getName() == cell_material_name; });
//
//     // Retrive all the section ids
//     std::vector<int> cell_section_ids =
//       util::parseString<int>(cell_node, "sections");
//
//     // Initialize cell sections
//     std::map<int, mesh::Section> cell_sections;
//
//     // For each id find the correct id in the sections vector and add to the
//     // cell sections vector
//     for (int& id : cell_section_ids) {
//       cell_sections.emplace(id, sections.find(id)->second);
//     }
//
//     // Construct Cell and place in vector
//     cells.emplace_back(cell_id, cell_sections, cell_material);
//   }
//
//   return cells;
// }
//
// mesh::Mesh parseMeshNode(const pugi::xml_node& root)
// {
//   // Parse sections node
//   std::map<int, mesh::Section> sections = util::parseSectionsNode(root);
//
//   // Parse cells node
//   std::vector<mesh::Cell> cells = util::parseCellNode(root, sections);
//
//   // Get mesh node from root
//   pugi::xml_node mesh_node = util::getNode(root, "mesh");
//
//   // Get cartesian_mesh node
//   pugi::xml_node cartesian_mesh_node =
//     util::getNode(mesh_node, "cartesian_mesh");
//
//   // Get total x and y bins in the mesh
//   size_t total_x_bins =
//     util::getAttribute<size_t>(cartesian_mesh_node, "total_x_bins");
//   size_t total_y_bins =
//     util::getAttribute<size_t>(cartesian_mesh_node, "total_y_bins");
//
//   // Get boundary conditions
//   std::vector<std::string> sides = {"left", "right", "top", "bottom"};
//   std::vector<std::pair<double, double>> bounds;
//   bounds.reserve(4);
//
//   // For each side get pair of a and b for boundary conditions
//   for (auto& side : sides) {
//     pugi::xml_node bound = util::getNode(mesh_node, side + "_condition");
//     bounds.emplace_back(util::getAttribute<double>(bound, "a"),
//       util::getAttribute<double>(bound, "b"));
//   }
//
//   // Initialize mesh
//   mesh::Mesh mesh(
//     total_x_bins, total_y_bins, bounds[0], bounds[1], bounds[2], bounds[3]);
//
//   // Initialize axes
//   std::vector<double> x_section_lengths;
//   std::vector<size_t> x_bins;
//
//   // For each xml node named x
//   assert(cartesian_mesh_node.child("x"));
//   for (const pugi::xml_node& axis_node : cartesian_mesh_node.children("x")) {
//     // Find the section
//     mesh::Section& section =
//       sections.find(util::getAttribute<int>(axis_node, "section"))->second;
//
//     // Add to section and bins vectors
//     x_section_lengths.push_back(section.getL1());
//     x_bins.push_back(util::getAttribute<size_t>(axis_node, "bins"));
//   }
//
//   // Fill mesh element dxs
//   mesh.initXAxis(x_section_lengths, x_bins);
//
//   // Initialize y-axis
//   std::vector<double> y_section_lengths;
//   std::vector<size_t> y_bins;
//
//   // For each xml node named y
//   assert(cartesian_mesh_node.child("y"));
//   for (const pugi::xml_node& axis_node : cartesian_mesh_node.children("y")) {
//     // Find the section
//     mesh::Section& section =
//       sections.find(util::getAttribute<int>(axis_node, "section"))->second;
//
//     // Add to section and bins vector
//     y_section_lengths.push_back(section.getL2());
//     y_bins.push_back(util::getAttribute<size_t>(axis_node, "bins"));
//   }
//
//   // Fill mesh element dys
//   mesh.initYAxis(y_section_lengths, y_bins);
//
//   // Build mesh with the cells
//   mesh.buildCells(cells);
//   return mesh;
// }
//
// rbm::PerturbAbsorption parseRBMNode(
//   const pugi::xml_node& root, mesh::Mesh& mesh)
// {
//   // Get rbm node
//   pugi::xml_node rbm_node = util::getNode(root, "rbm");
//
//   // Get training node
//   pugi::xml_node training_node = util::getNode(rbm_node, "training");
//
//   // Get cell id that will be perturbed
//   int cell_id = util::getAttribute<int>(training_node, "cell_id");
//
//   // Parse training fluxes
//   xt::xarray<double> training_points =
//     xt::adapt(util::parseString<double>(training_node, "values"));
//
//   return rbm::PerturbAbsorption(training_points, mesh, cell_id);
// }

} // namespace util
