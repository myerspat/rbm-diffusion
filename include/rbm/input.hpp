#ifndef _INPUT_
#define _INPUT_

#include "rbm/material.hpp"
#include "rbm/mesh.hpp"
#include "rbm/rbm.hpp"
#include <cassert>
#include <map>
#include <pugixml.hpp>
#include <string>
#include <vector>

namespace util {

// Load xml_document
void loadFile(pugi::xml_document& file, const std::string& file_name);

// Get child of pugi node
pugi::xml_node getNode(
  const pugi::xml_node& node, const std::string& child_name);

// Parse materials node
std::vector<Material> parseMaterialsNode(const pugi::xml_node& root);

// Parse mesh node and build Mesh
mesh::Mesh parseMeshNode(const pugi::xml_node& root);

// Parse rbm node and build Perturb class
rbm::Perturb parseRBMNode(
  const pugi::xml_node& root, mesh::Mesh& mesh);

// Get node attribute as type T
template<typename T>
T getAttribute(const pugi::xml_node& node, const std::string& attribute)
{
  auto attr = node.attribute(attribute.c_str());
  if (attr) {
    if constexpr (std::is_same_v<T, double>)
      return attr.as_double();
    else if constexpr (std::is_same_v<T, int> || std::is_same_v<T, size_t>)
      return attr.as_int();
    else if constexpr (std::is_same_v<T, std::string>)
      return attr.as_string();
  }
  throw std::runtime_error("Node " + std::string(node.name()) +
                           " does not have an attribute named " + attribute);
}

// Parse an XML child node
template<typename T>
std::vector<T> parseString(const pugi::xml_node& node, const std::string& value,
  const std::string& delimiter = " ")
{
  // Assert the child node exists
  assert(node.child(value.c_str()));

  // Grab a string of the child values
  std::string node_text = node.child_value(value.c_str());

  // If the values vector is empty
  if (node_text.empty()) {
    // If the vector is empty throw an error
    throw std::runtime_error(
      "No values were found within " + value + " node in " + node.name());
  }

  // Initialize token, node_values, and pos
  std::string token;
  std::vector<T> node_values;
  size_t pos = 0;

  // While node_text is populated
  while ((pos = node_text.find(delimiter)) != std::string::npos) {
    // Get the substring
    token = node_text.substr(0, pos);

    // Add values based on type
    if constexpr (std::is_same_v<T, double>)
      node_values.push_back(std::stod(token));
    else if constexpr (std::is_same_v<T, int> || std::is_same_v<T, size_t>)
      node_values.push_back(std::stoi(token));
    else if constexpr (std::is_same_v<T, std::string>)
      node_values.push_back(token);

    // Remove the value in the text
    node_text.erase(0, pos + delimiter.length());
  }

  // Add the last value
  if (!node_text.empty()) {
    // Add values based on type
    if constexpr (std::is_same_v<T, double>)
      node_values.push_back(std::stod(node_text));
    else if constexpr (std::is_same_v<T, int> || std::is_same_v<T, size_t>)
      node_values.push_back(std::stoi(node_text));
    else if constexpr (std::is_same_v<T, std::string>)
      node_values.push_back(token);
  }

  return node_values;
}

} // namespace util

#endif // !_INPUT_
