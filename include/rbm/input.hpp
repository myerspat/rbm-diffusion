#ifndef _INPUT_
#define _INPUT_

#include <pugixml.hpp>
#include <string>

namespace util {

// Load xml_document
void loadFile(pugi::xml_document& file, const std::string& file_name);

} // namespace util

#endif // !_INPUT_
