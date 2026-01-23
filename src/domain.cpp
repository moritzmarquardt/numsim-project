#include "domain.hpp"

#include <fstream>
#include <stdexcept>

CellType cellTypeFromChar(char c) {
  switch (c) {
  case '-':
    return CellType::Fluid;
  case 'x':
    return CellType::Obstacle;
  case 'i':
    return CellType::Inflow;
  case 'o':
    return CellType::Outflow;
  case 'n':
    return CellType::NoSlip;
  case 's':
    return CellType::Slip;
  default:
    throw std::runtime_error("Unknown domain character");
  }
}

char cellTypeToChar(CellType type) {
  switch (type) {
  case CellType::Fluid:
    return '-';
  case CellType::Obstacle:
    return 'x';
  case CellType::Inflow:
    return 'i';
  case CellType::Outflow:
    return 'o';
  case CellType::NoSlip:
    return 'n';
  case CellType::Slip:
    return 's';
  default:
    return '?';
  }
}

void Domain::loadFromFile(const std::string &filename, std::array<int, 2> nCells) {
  nCells_ = nCells;
  cells_.assign(static_cast<size_t>(nCells_[0] * nCells_[1]), CellType::Fluid);

  std::ifstream file(filename.c_str(), std::ios::in);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open domain file: " + filename);
  }

  std::string line;
  int row = 0;
  while (std::getline(file, line)) {
    // strip whitespace
    std::string compact;
    compact.reserve(line.size());
    for (char c : line) {
      if (c == '#') {
        break;
      }
      if (c != ' ' && c != '\t' && c != '\r') {
        compact.push_back(c);
      }
    }
    if (compact.empty()) {
      continue;
    }
    if (static_cast<int>(compact.size()) != nCells_[0]) {
      throw std::runtime_error("Domain row width does not match nCellsX");
    }
    if (row >= nCells_[1]) {
      throw std::runtime_error("Domain has more rows than nCellsY");
    }

    for (int i = 0; i < nCells_[0]; ++i) {
      cells_[static_cast<size_t>(row * nCells_[0] + i)] = cellTypeFromChar(compact[i]);
    }
    row++;
  }

  if (row != nCells_[1]) {
    throw std::runtime_error("Domain row count does not match nCellsY");
  }
}

CellType Domain::cellTypeGlobal(int i, int j) const {
  if (i < 0 || j < 0 || i >= nCells_[0] || j >= nCells_[1]) {
    return CellType::Obstacle;
  }
  return cells_[static_cast<size_t>(j * nCells_[0] + i)];
}

std::array<int, 2> Domain::nCells() const {
  return nCells_;
}
