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
    ghostListLocal_ = std::make_unique<std::vector<CellInfo>>();

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

    // Create no-slip boundary condition codes if they don't exist yet
    double noSlipRightCode = -1.0;
    double noSlipTopCode = -1.0;
    
    // Check if a no-slip BC already exists in the maps
    for (const auto& [code, info] : rightFaceBCInfo_) {
        if (info.dirichletU.has_value() && info.dirichletU.value() == 0.0 &&
            info.dirichletV.has_value() && info.dirichletV.value() == 0.0 &&
            !info.neumannU.has_value() && !info.neumannV.has_value()) {
            noSlipRightCode = code;
            break;
        }
    }
    
    // If no no-slip code exists for right faces, create one
    if (noSlipRightCode < 0.0) {
        noSlipRightCode = nextRightCode++;
        BoundaryInfo noSlipInfo;
        noSlipInfo.dirichletU = 0.0;
        noSlipInfo.dirichletV = 0.0;
        rightFaceBCInfo_[noSlipRightCode] = noSlipInfo;
        rightFaceCodeToMarker_[noSlipRightCode] = 'O'; // 'O' for obstacle
    }
    
    // Same for top faces
    for (const auto& [code, info] : topFaceBCInfo_) {
        if (info.dirichletU.has_value() && info.dirichletU.value() == 0.0 &&
            info.dirichletV.has_value() && info.dirichletV.value() == 0.0 &&
            !info.neumannU.has_value() && !info.neumannV.has_value()) {
            noSlipTopCode = code;
            break;
        }
    }
    
    if (noSlipTopCode < 0.0) {
        noSlipTopCode = nextTopCode++;
        BoundaryInfo noSlipInfo;
        noSlipInfo.dirichletU = 0.0;
        noSlipInfo.dirichletV = 0.0;
        topFaceBCInfo_[noSlipTopCode] = noSlipInfo;
        topFaceCodeToMarker_[noSlipTopCode] = 'O'; // 'O' for obstacle
    }

    // Apply default no-slip boundary conditions for obstacle faces that have no defined BC
    for (int j = 0; j < nCellsY; ++j) {
        for (int i = 0; i < nCellsX; ++i) {
            bool isFluid = ((*obstacleMaskGlobal_)(i, j) < 0.5);
            
            if (isFluid) {
                // Check right face (at i+1, j)
                if (i + 1 < nCellsX) {
                    bool rightIsObstacle = ((*obstacleMaskGlobal_)(i + 1, j) > 0.5);
                    if (rightIsObstacle) {
                        double faceCode = (*rightFacesBCGlobal_)(i + 1, j);
                        const BoundaryInfo& faceInfo = rightFaceBCInfo_[faceCode];
                        // Only set no-slip if no boundary conditions are defined
                        if (!faceInfo.dirichletU.has_value() && !faceInfo.neumannU.has_value() &&
                            !faceInfo.dirichletV.has_value() && !faceInfo.neumannV.has_value()) {
                            (*rightFacesBCGlobal_)(i + 1, j) = noSlipRightCode;
                        }
                    }
                }
                
                // Check left face (at i, j)
                if (i - 1 >= 0) {
                    bool leftIsObstacle = ((*obstacleMaskGlobal_)(i - 1, j) > 0.5);
                    if (leftIsObstacle) {
                        double faceCode = (*rightFacesBCGlobal_)(i, j);
                        const BoundaryInfo& faceInfo = rightFaceBCInfo_[faceCode];
                        if (!faceInfo.dirichletU.has_value() && !faceInfo.neumannU.has_value() &&
                            !faceInfo.dirichletV.has_value() && !faceInfo.neumannV.has_value()) {
                            (*rightFacesBCGlobal_)(i, j) = noSlipRightCode;
                        }
                    }
                }
                
                // Check top face (at i, j+1)
                if (j + 1 < nCellsY) {
                    bool topIsObstacle = ((*obstacleMaskGlobal_)(i, j + 1) > 0.5);
                    if (topIsObstacle) {
                        double faceCode = (*topFacesBCGlobal_)(i, j + 1);
                        const BoundaryInfo& faceInfo = topFaceBCInfo_[faceCode];
                        if (!faceInfo.dirichletU.has_value() && !faceInfo.neumannU.has_value() &&
                            !faceInfo.dirichletV.has_value() && !faceInfo.neumannV.has_value()) {
                            (*topFacesBCGlobal_)(i, j + 1) = noSlipTopCode;
                        }
                    }
                }
                
                // Check bottom face (at i, j)
                if (j - 1 >= 0) {
                    bool bottomIsObstacle = ((*obstacleMaskGlobal_)(i, j - 1) > 0.5);
                    if (bottomIsObstacle) {
                        double faceCode = (*topFacesBCGlobal_)(i, j);
                        const BoundaryInfo& faceInfo = topFaceBCInfo_[faceCode];
                        if (!faceInfo.dirichletU.has_value() && !faceInfo.neumannU.has_value() &&
                            !faceInfo.dirichletV.has_value() && !faceInfo.neumannV.has_value()) {
                            (*topFacesBCGlobal_)(i, j) = noSlipTopCode;
                        }
                    }
                }
            }
        }
    }

    // make local cell lists
    int nCellsXLocal = partitioning_->nCellsLocal()[0];
    int nCellsYLocal = partitioning_->nCellsLocal()[1];
    int xOffset = partitioning_->nodeOffset()[0];
    int yOffset = partitioning_->nodeOffset()[1];

    for (int jLocal = 0; jLocal < nCellsYLocal; ++jLocal) {
        int jGlobal = jLocal + yOffset;
        for (int iLocal = 0; iLocal < nCellsXLocal; ++iLocal) {
            int iGlobal = iLocal + xOffset;

            CellInfo cellInfo = createCellInfo(iGlobal, jGlobal, iLocal, jLocal, nCellsXLocal, nCellsYLocal);
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

    // Fill ghostListLocal_ with cells adjacent to left and bottom partition boundaries
    // (only if those boundaries are NOT domain boundaries)
    
    // Add cells adjacent to left partition boundary (if left is not domain boundary)
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int jLocal = 0; jLocal < nCellsYLocal; ++jLocal) {
            int jGlobal = jLocal + yOffset;
            int iLocal = - 1;  // leftmost column in partition
            int iGlobal = iLocal + xOffset;
            std::cout << "Adding ghost cell at global (" << iGlobal << ", " << jGlobal << "), local (" << iLocal << ", " << jLocal << ")" << std::endl;
            
            CellInfo cellInfo = createCellInfo(iGlobal, jGlobal, iLocal, jLocal, nCellsXLocal, nCellsYLocal);
            std::cout << "  CellInfo: " << cellInfo.toString() << std::endl;
            ghostListLocal_->push_back(cellInfo);
        }
    }
    
    // Add cells adjacent to bottom partition boundary (if bottom is not domain boundary)
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int iLocal = 0; iLocal < nCellsXLocal; ++iLocal) {
            int iGlobal = iLocal + xOffset;
            int jLocal = - 1;  // bottommost row in partition
            int jGlobal = jLocal + yOffset;
            std::cout << "Adding ghost cell at global (" << iGlobal << ", " << jGlobal << "), local (" << iLocal << ", " << jLocal << ")" << std::endl;
            
            CellInfo cellInfo = createCellInfo(iGlobal, jGlobal, iLocal, jLocal, nCellsXLocal, nCellsYLocal);
            std::cout << "  CellInfo: " << cellInfo.toString() << std::endl;
            
            ghostListLocal_->push_back(cellInfo);
        }
    }
}


// helper function to fill cellinfo
CellInfo Domain::createCellInfo(int iGlobal, int jGlobal, int iLocal, int jLocal, int nCellsXLocal, int nCellsYLocal) {
    CellInfo cellInfo;
    // cellInfo.cellIndexGlobal = {iGlobal, jGlobal};
    cellInfo.cellIndexPartition = {iLocal + 2, jLocal + 2};
    cellInfo.fluidCell = ((*obstacleMaskGlobal_)(iGlobal, jGlobal) == 0.0);
    cellInfo.faceRight = rightFaceBCInfo_.at((*rightFacesBCGlobal_)(iGlobal+1, jGlobal));
    cellInfo.faceTop = topFaceBCInfo_.at((*topFacesBCGlobal_)(iGlobal, jGlobal+1));
    cellInfo.faceBottom = topFaceBCInfo_.at((*topFacesBCGlobal_)(iGlobal, jGlobal));
    cellInfo.faceLeft = rightFaceBCInfo_.at((*rightFacesBCGlobal_)(iGlobal, jGlobal));
    cellInfo.faceRight.isPartitionInnerFace = partitioning_->ownPartitionContainsRightBoundary() == false && (iLocal == nCellsXLocal - 1);
    cellInfo.faceTop.isPartitionInnerFace = partitioning_->ownPartitionContainsTopBoundary() == false && (jLocal == nCellsYLocal - 1);
    cellInfo.faceBottom.isPartitionInnerFace = partitioning_->ownPartitionContainsBottomBoundary() == false && (jLocal == 0);
    cellInfo.faceLeft.isPartitionInnerFace = partitioning_->ownPartitionContainsLeftBoundary() == false && (iLocal == 0);
    return cellInfo;
}