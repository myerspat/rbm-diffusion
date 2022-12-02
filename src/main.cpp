#include "rbm/input.hpp"
#include "rbm/material.hpp"
#include "rbm/mesh.hpp"
#include "rbm/rbm.hpp"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <pugixml.hpp>
#include <string>
#include <xtensor/xio.hpp>

int main(int argc, char* argv[])
{
  // Top print statements
  std::cout
    << "\n rbm-diffusion\n"
    << "\n University of Michigan\n"
    << " Nuclear Engineering and Radiological Sciences 570\n"
    << " Authors: Patrick Myers, Connor Craig, Erin Burrell, Kyle Beyer\n";

  std::chrono::time_point<std::chrono::system_clock> start, end;

  // Start timing offline section
  start = std::chrono::system_clock::now();

  // Default file path
  std::string xml_file_name = "rbm.xml";
  if (argc > 1) {
    xml_file_name = argv[1];
  }

  // Load XML document
  pugi::xml_document xml_file;
  util::loadFile(xml_file, xml_file_name);

  // Print file path
  std::cout << "\n XML Input file found: " << xml_file_name
            << "\n Begin parsing input file";

  // Get rbm root node
  assert(xml_file.child("rbm"));
  pugi::xml_node rbm_root = util::getNode(xml_file, "rbm");

  // Build mesh from xml file
  mesh::Mesh mesh = util::parseMeshNode(rbm_root);

  // Build target pertubation parameter
  rbm::Perturb parameter = util::parseRBMNode(rbm_root, mesh);

  // Build subspace
  parameter.train();

  // End offline section timing
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> offline_time = end - start;

  // Calculate eigenvectors/values for target parameter values
  start = std::chrono::system_clock::now();
  parameter.calcTargets();
  end = std::chrono::system_clock::now();

  // Find online time
  std::chrono::duration<double> online_time = end - start;

  // Calculate error
  parameter.checkError();

  // Print offline timing results
  std::cout << std::setprecision(5)
            << "\n Offline Time = " << offline_time.count()
            << " s | Average Time per Training Point = "
            << offline_time.count() / parameter.getNumTraining() << " s"
            << std::endl;

  // Print online timing results
  std::cout << std::setprecision(5) << " Online Time = " << online_time.count()
            << " s | Average Time per Target Point = "
            << online_time.count() / parameter.getNumTarget() << " s"
            << std::endl;

  // Print finishing time
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << " Finished computation on " << std::ctime(&end_time)
            << std::endl;

  return 0;
}
