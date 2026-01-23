#pragma once

#include <array>
#include <string>
#include <vector>

enum class CellType : unsigned char {
  Fluid,
  Obstacle,
  Inflow,
  Outflow,
  NoSlip,
  Slip
};

CellType cellTypeFromChar(char c);
char cellTypeToChar(CellType type);

class Domain {
public:
  void loadFromFile(const std::string &filename, std::array<int, 2> nCells);
  CellType cellTypeGlobal(int i, int j) const;
  std::array<int, 2> nCells() const;

private:
  std::array<int, 2> nCells_{0, 0};
  std::vector<CellType> cells_;
};
