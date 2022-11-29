#include "rbm/input.hpp"
#include "rbm/material.hpp"
#include "rbm/mesh.hpp"
#include "rbm/rbm.hpp"
#include <algorithm>
#include <pugixml.hpp>
#include <string>
#include <xtensor/xio.hpp>

int main(int argc, char* argv[])
{
  // Default file path
  std::string xml_file_name = "rbm.xml";
  if (argc > 1) {
    xml_file_name = argv[1];
  }

  // Load XML document
  pugi::xml_document xml_file;
  util::loadFile(xml_file, xml_file_name);

  // Get rbm root node
  assert(xml_file.child("rbm"));
  pugi::xml_node rbm_root = util::getNode(xml_file, "rbm");

  // Build mesh from xml file
  mesh::Mesh mesh = util::parseMeshNode(rbm_root);

  // Build target pertubation parameter
  rbm::Perturb parameter = util::parseRBMNode(rbm_root, mesh);

  // Get target node and parse target values
  pugi::xml_node target_node = util::getNode(rbm_root.child("rbm"), "target");
  std::vector<double> target_values =
    util::parseString<double>(target_node, "values");

  // Allocate space for eigenvectors and eigenvalues
  xt::xarray<double> target_fluxes =
    xt::xarray<double>::from_shape({target_values.size(), mesh.getSize()});
  xt::xarray<double> target_k =
    xt::xarray<double>::from_shape({target_values.size()});

  for (size_t i = 0; i < target_values.size(); i++) {
    // Calculate the eigenvector and eigenvalue for parameter target
    std::pair<xt::xarray<double>, double> target =
      parameter.calcTarget(target_values[i]);

    // Fill target eigenvector array
    target_k(i) = target.second;
  }

  // Print target k values
  std::cout << target_k << std::endl;

  return 0;
}
