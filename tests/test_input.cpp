#include "rbm/input.hpp"
#include "rbm/mesh.hpp"
#include "rbm/rbm.hpp"
#include "unit_test_framework.hpp"
#include <iostream>
#include <pugixml.hpp>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor/xmath.hpp>

TEST(test_loadFile_1)
{
  // Initialize doc and supporting file
  pugi::xml_document doc;
  std::string supporting_file = "../../tests/supporting/test.xml";

  // Load doc
  util::loadFile(doc, supporting_file);

  // Assert top level node is rbm
  ASSERT_TRUE(doc.child("rbm"));
}

TEST(test_getNode_1)
{
  // Initialize doc and supporting file
  pugi::xml_document doc;
  std::string supporting_file = "../../tests/supporting/test.xml";

  // Load doc
  util::loadFile(doc, supporting_file);

  ASSERT_EQUAL(std::string(util::getNode(doc, "rbm").name()), "rbm")
}

TEST(test_getAttribute_1)
{
  // Initialize doc and supporting file
  pugi::xml_document doc;
  std::string supporting_file = "../../tests/supporting/test.xml";

  // Load doc
  util::loadFile(doc, supporting_file);
  auto material = doc.child("rbm").child("materials").child("material");

  // Assert attributes are expected values
  ASSERT_EQUAL(util::getAttribute<double>(material, "absorption"), 0.10);
  ASSERT_EQUAL(util::getAttribute<std::string>(material, "name"), "fuel");
}

TEST(test_parseString_1)
{
  // Fill node such that
  // <node>
  //    <values>0 1 2 3 4</values>
  // </node>
  pugi::xml_document doc;
  pugi::xml_node node = doc.append_child("node");
  pugi::xml_node values_node = node.append_child("values");
  values_node.append_child(pugi::node_pcdata).set_value("0 1 2 3 4");

  // Parse string to get a vector
  std::vector<int> values = util::parseString<int>(node, "values");

  // Asssert the data read in is 0 1 2 3 4
  for (int i = 0; i < 5; i++) {
    ASSERT_EQUAL(values[i], i);
  }
}

TEST(test_parseMaterialsNode_1)
{
  // Initialize doc and supporting file
  pugi::xml_document doc;
  std::string supporting_file = "../../tests/supporting/test.xml";

  // Load doc
  util::loadFile(doc, supporting_file);

  // Get rbm top level node
  auto root = util::getNode(doc, "rbm");

  // Parse materials node
  std::vector<material::Material> materials = util::parseMaterialsNode(root);

  // Assertions
  ASSERT_EQUAL(materials[0].getName(), "fuel");
  ASSERT_EQUAL(materials[0].getAbsorption(), 0.10);
  ASSERT_EQUAL(materials[0].getNuFission(), 0.11);
  ASSERT_EQUAL(materials[0].getD(), 2.0);
  ASSERT_EQUAL(materials[1].getName(), "reflector");
  ASSERT_EQUAL(materials[1].getAbsorption(), 0.01);
  ASSERT_EQUAL(materials[1].getNuFission(), 1e-7);
  ASSERT_EQUAL(materials[1].getD(), 1.5);
}

TEST(test_parseMeshNode_1)
{
  // Initialize doc and supporting file
  pugi::xml_document doc;
  std::string supporting_file = "../../tests/supporting/test.xml";

  // Load doc
  util::loadFile(doc, supporting_file);

  // Get rbm top level node
  auto root = util::getNode(doc, "rbm");

  // Parse and build mesh
  mesh::Mesh mesh = util::parseMeshNode(root);

  ASSERT_EQUAL(mesh.getXN(), 40);
  ASSERT_EQUAL(mesh.getYN(), 40);
  ASSERT_EQUAL(mesh.getSize(), 1600);
}

TEST(test_parseRBMNode_1)
{
  // Initialize doc and supporting file
  pugi::xml_document doc;
  std::string supporting_file = "../../tests/supporting/test.xml";

  // Load doc
  util::loadFile(doc, supporting_file);

  // Get rbm top level node
  auto root = util::getNode(doc, "rbm");

  // Parse and build mesh
  mesh::Mesh mesh = util::parseMeshNode(root);

  // Parse RBM node and initialize
  rbm::Perturb parameter = util::parseRBMNode(root, mesh);

  // No assertions included for this test as the RBM class does not have any
  // accessors at the moment
}

TEST_MAIN();
