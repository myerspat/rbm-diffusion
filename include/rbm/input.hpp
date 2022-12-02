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
/// loadFile() takes in xml file to run the code
///
/// @params file is a xml file object to open and read xml information
/// @params file_name is a string to the path of the xml file
/// @warnings if file not found then it throws a runtime error
void loadFile(pugi::xml_document& file, const std::string& file_name);

/// Get child of pugi node
/// @params node is the xml node found when loading the file
/// @params child_name is the name of the child to be found
///@warnings If the child node cant be found then it throws a runtime error
pugi::xml_node getNode(
  const pugi::xml_node& node, const std::string& child_name);

/// Parse materials node
/// @params root is the root xml node
///@retusn a vector of materials
std::vector<Material> parseMaterialsNode(const pugi::xml_node& root);

/// Parse mesh node and build Mesh
/// @params root is the root xml node
/// @returns mesh information
mesh::Mesh parseMeshNode(const pugi::xml_node& root);

/// Parse rbm node and build Perturb class
/// @params root takes in information to make Perturb object via the root xml
/// node
/// @params mesh the built mesh object is used to construct the Perturb Object
/// @returns a Purturb oject based off the information from in the xml file.
rbm::Perturb parseRBMNode(const pugi::xml_node& root, mesh::Mesh& mesh);

/// Get node attribute as type T
/// @params node is a xml node from input file
/// @params attribute is the attribute of the specific node being used
/// @warnings If node does not have an attribute then a runtime error is thrown.
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

/// Parse an XML child node
/// @params node is an xml from the input file
/// @value is the value being looked for in the node
/// @warnings If no value is found in the node then throws a runtime error.
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
