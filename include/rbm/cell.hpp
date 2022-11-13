#ifndef _CELL_
#define _CELL_

#include "rbm/material.hpp"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace mesh {

class Section {
public:
  //========================================================================
  // Constructor / Destructor
  Section() {};
  Section(int id, std::pair<double, double> p0, double l1, double l2)
    : _id(id), _p0(p0), _l1(l1), _l2(l2) {};

  //========================================================================
  // Methods
  void initialize(
    int id, std::pair<double, double> p0, double l1, double l2);

  //========================================================================
  // Getters
  int getID() { return _id; };
  std::pair<double, double> getP0() { return _p0; };
  double getL1() { return _l1; };
  double getL2() { return _l2; };

private:
  //========================================================================
  // Data
  int _id;                       // Section ID
  std::pair<double, double> _p0; // Bottom left point of the rectangle
  double _l1;                    // Length along the x-axis
  double _l2;                    // Length along the y-axis
};

class Cell {
public:
  //========================================================================
  // Constructor / Destructor
  Cell() {};
  Cell(int id) : _id(id) {};
  Cell(int id, std::vector<Section> sections, std::vector<Material> materials)
    : _id(id), _sections(sections), _materials(materials) {};

  //========================================================================
  // Methods
  void addSection(Section section, std::size_t x_bins, std::size_t y_bins);

  //========================================================================
  // Getters
  int getID() { return _id; };
  Section getSection(std::size_t idx) { return _sections[idx]; };
  Material getMaterial(std::size_t idx) { return _materials[idx]; };

private:
  //========================================================================
  // Data
  int _id;
  std::vector<Section> _sections;
  std::vector<Material> _materials;
};

} // namespace mesh

#endif // !_CELL_
