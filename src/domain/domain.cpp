#include "domain.hpp"
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>

Domain::Domain(const Settings* settings, std::shared_ptr<Partitioning> partitioning)
    : settings_(settings), partitioning_(partitioning) {
    // Constructor implementation is empty, nothing has to be done here. the interesting stuff happens in readDomainFile
}

void Domain::readDomainFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open domain file: " << filename << std::endl;
        return;
    }

    const int nCellsX = settings_->nCells[0];
    const int nCellsY = settings_->nCells[1];

    obstacleMaskGlobal_ = std::make_unique<Array2D>(std::array<int, 2>{nCellsX, nCellsY});
    rightFacesBCGlobal_ = std::make_unique<Array2D>(std::array<int, 2>{nCellsX + 1, nCellsY});
    topFacesBCGlobal_ = std::make_unique<Array2D>(std::array<int, 2>{nCellsX, nCellsY + 1});

    obstacleMaskGlobal_->setToZero();
    rightFacesBCGlobal_->setToZero();
    topFacesBCGlobal_->setToZero();

    cellListAllLocal_ = std::make_unique<std::vector<CellInfo>>();
    cellListFluidLocal_ = std::make_unique<std::vector<CellInfo>>();
    redListFluidLocal_ = std::make_unique<std::vector<CellInfo>>();
    blackListFluidLocal_ = std::make_unique<std::vector<CellInfo>>();

    rightFaceBCInfo_.clear();
    topFaceBCInfo_.clear();
    rightFaceCodeToMarker_.clear();
    topFaceCodeToMarker_.clear();

    std::unordered_map<char, double> rightCharToCode{{'-', 0.0}};
    std::unordered_map<char, double> topCharToCode{{'-', 0.0}};
    double nextRightCode = 1.0;
    double nextTopCode = 1.0;

    auto parseValue = [](const std::string &line, const std::string &key) -> std::optional<double> {
        const std::size_t pos = line.find(key);
        if (pos == std::string::npos) {
            return std::nullopt;
        }
        const std::size_t valuePos = pos + key.size();
        std::size_t endPos = line.find_first_of(",; ", valuePos);
        const std::string valueStr = line.substr(valuePos, endPos - valuePos);
        try {
            return std::stod(valueStr);
        } catch (...) {
            return std::nullopt;
        }
    };

    auto ensureRightMapping = [&](char marker, const BoundaryInfo &info) {
        auto it = rightCharToCode.find(marker);
        if (it == rightCharToCode.end()) {
            const double code = nextRightCode++;
            rightCharToCode[marker] = code;
            rightFaceCodeToMarker_[code] = marker;
            rightFaceBCInfo_[code] = info;
        } else {
            const double code = it->second;
            rightFaceBCInfo_[code] = info;
            rightFaceCodeToMarker_[code] = marker;
        }
    };

    auto ensureTopMapping = [&](char marker, const BoundaryInfo &info) {
        auto it = topCharToCode.find(marker);
        if (it == topCharToCode.end()) {
            const double code = nextTopCode++;
            topCharToCode[marker] = code;
            topFaceCodeToMarker_[code] = marker;
            topFaceBCInfo_[code] = info;
        } else {
            const double code = it->second;
            topFaceBCInfo_[code] = info;
            topFaceCodeToMarker_[code] = marker;
        }
    };

    ensureRightMapping('-', BoundaryInfo{});
    ensureTopMapping('-', BoundaryInfo{});

    auto parseLegendLine = [&](const std::string &line, bool isRight) {
        if (line.size() < 2 || line[1] != ':') {
            return false;
        }
        const char marker = line[0];
        BoundaryInfo info;
        // accept both u_D= and u=
        if (auto val = parseValue(line, "u_D=")) info.dirichletU = val;
        if (auto val = parseValue(line, "v_D=")) info.dirichletV = val;
        if (auto val = parseValue(line, "u_N=")) info.neumannU = val;
        if (auto val = parseValue(line, "v_N=")) info.neumannV = val;
        if (auto val = parseValue(line, "u=")) info.dirichletU = val;
        if (auto val = parseValue(line, "v=")) info.dirichletV = val;

        if (isRight) {
            ensureRightMapping(marker, info);
        } else {
            ensureTopMapping(marker, info);
        }
        return true;
    };

    auto encodeRightMarker = [&](char marker) {
        auto it = rightCharToCode.find(marker);
        if (it == rightCharToCode.end()) {
            ensureRightMapping(marker, BoundaryInfo{});
            it = rightCharToCode.find(marker);
        }
        return it->second;
    };

    auto encodeTopMarker = [&](char marker) {
        auto it = topCharToCode.find(marker);
        if (it == topCharToCode.end()) {
            ensureTopMapping(marker, BoundaryInfo{});
            it = topCharToCode.find(marker);
        }
        return it->second;
    };

    enum class Section { None, Obstacles, FacesRight, FacesTop };
    Section currentSection = Section::None;
    int rowCount = 0;
    const int obstacleRows = nCellsY;
    const int rightRows = nCellsY;
    const int topRows = nCellsY + 1;

    auto trimTrailingCarriageReturn = [](std::string &value) {
        if (!value.empty() && value.back() == '\r') {
            value.pop_back();
        }
    };

    std::string line;
    while (std::getline(file, line)) {
        trimTrailingCarriageReturn(line);

        if (line.empty()) {
            continue;
        }

        if (line.rfind("ObstacleMarker", 0) == 0) {
            currentSection = Section::Obstacles;
            rowCount = 0;
            continue;
        }
        if (line.rfind("FacesRight", 0) == 0) {
            currentSection = Section::FacesRight;
            rowCount = 0;
            continue;
        }
        if (line.rfind("FacesTop", 0) == 0) {
            currentSection = Section::FacesTop;
            rowCount = 0;
            continue;
        }

        if (!line.empty() && line[0] == '#') {
            continue;
        }

        switch (currentSection) {
            case Section::Obstacles: {
                if (rowCount >= obstacleRows) {
                    continue;
                }
                const int j = nCellsY - 1 - rowCount;
                for (int i = 0; i < nCellsX; ++i) {
                    const char marker = (i < static_cast<int>(line.size())) ? line[i] : '-';
                    const double value = (marker == '#') ? 1.0 : 0.0;
                    (*obstacleMaskGlobal_)(i, j) = value;
                }
                rowCount++;
                break;
            }
            case Section::FacesRight: {
                if (parseLegendLine(line, true)) {
                    continue;
                }
                if (rowCount >= rightRows) {
                    continue;
                }
                const int j = nCellsY - 1 - rowCount;
                for (int i = 0; i < nCellsX + 1; ++i) {
                    const char marker = (i < static_cast<int>(line.size())) ? line[i] : '-';
                    const double encoded = encodeRightMarker(marker);
                    (*rightFacesBCGlobal_)(i, j) = encoded;
                }
                rowCount++;
                break;
            }
            case Section::FacesTop: {
                if (parseLegendLine(line, false)) {
                    continue;
                }
                if (rowCount >= topRows) {
                    continue;
                }
                const int j = nCellsY - rowCount;
                for (int i = 0; i < nCellsX; ++i) {
                    const char marker = (i < static_cast<int>(line.size())) ? line[i] : '-';
                    const double encoded = encodeTopMarker(marker);
                    (*topFacesBCGlobal_)(i, j) = encoded;
                }
                rowCount++;
                break;
            }
            case Section::None:
            default:
                break;
        }
    }

    file.close();


    // make local cell lists
    int nCellsXLocal = partitioning_->nCellsLocal()[0];
    int nCellsYLocal = partitioning_->nCellsLocal()[1];
    int xOffset = partitioning_->nodeOffset()[0];
    int yOffset = partitioning_->nodeOffset()[1];

    for (int jLocal = 0; jLocal < nCellsYLocal; ++jLocal) {
        int jGlobal = jLocal + yOffset;
        for (int iLocal = 0; iLocal < nCellsXLocal; ++iLocal) {
            int iGlobal = iLocal + xOffset;
            CellInfo cellInfo;
            // cellInfo.cellIndexGlobal = {iGlobal, jGlobal};
            cellInfo.cellIndexPartition = {iLocal, jLocal};
            cellInfo.fluidCell = ((*obstacleMaskGlobal_)(iGlobal, jGlobal) == 0.0);
            cellInfo.faceRight = rightFaceBCInfo_.at((*rightFacesBCGlobal_)(iGlobal+1, jGlobal));
            cellInfo.faceTop = topFaceBCInfo_.at((*topFacesBCGlobal_)(iGlobal, jGlobal+1));
            cellInfo.faceBottom = topFaceBCInfo_.at((*topFacesBCGlobal_)(iGlobal, jGlobal));
            cellInfo.faceLeft = rightFaceBCInfo_.at((*rightFacesBCGlobal_)(iGlobal, jGlobal));

            cellListAllLocal_->push_back(cellInfo);


            if (cellInfo.fluidCell) {
                cellListFluidLocal_->push_back(cellInfo);
                if ((iGlobal + jGlobal) % 2 == 0) {
                    redListFluidLocal_->push_back(cellInfo);
                } else {
                    blackListFluidLocal_->push_back(cellInfo);
                }
            }
        }
    }
}


