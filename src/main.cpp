#include "rbm/input.hpp"
#include "rbm/mesh.hpp"
#include <string>
#include <pugixml.hpp>

int main (int argc, char *argv[])
{
  // Default file path
  std::string xml_file_name = "rbm.xml";
  if (argc > 1) {
    xml_file_name = argv[1];
  }  

  pugi::xml_document xml_file;
  util::loadFile(xml_file, xml_file_name);

  assert(xml_file.child("rbm"));




  return 0;
}
