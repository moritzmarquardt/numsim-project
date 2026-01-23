#include "staggeredGrid.hpp"

// nCells has to be N + 2, where N is the number of cells in x or y direction
StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth, std::shared_ptr<Partitioning> partitioning)
        : nCells_(nCells), meshWidth_(meshWidth),
            u_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-1 * meshWidth[0], -1.5 * meshWidth[1]}, meshWidth)),
            f_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-1 * meshWidth[0], -1.5 * meshWidth[1]}, meshWidth)),
            v_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-1.5 * meshWidth[0], -1 * meshWidth[1]}, meshWidth)),
            g_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-1.5 * meshWidth[0], -1 * meshWidth[1]}, meshWidth)),
            p_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-1.5 * meshWidth[0], -1.5 * meshWidth[1]}, meshWidth)),
            rhs_(FieldVariable({nCells[0] + 3, nCells[1] + 3}, {-1.5 * meshWidth[0], -1.5 * meshWidth[1]}, meshWidth)),
            partitioning_ (partitioning),
            cellTypes_(static_cast<size_t>((nCells[0] + 3) * (nCells[1] + 3)), CellType::Obstacle){
                containsLeftBoundary_ = partitioning_->ownPartitionContainsLeftBoundary();
                containsRightBoundary_ = partitioning_->ownPartitionContainsRightBoundary();
                containsBottomBoundary_ = partitioning_->ownPartitionContainsBottomBoundary();
                containsTopBoundary_ = partitioning_->ownPartitionContainsTopBoundary();
            }

const std::array<double,2> StaggeredGrid::meshWidth() const {
    return meshWidth_;
}

/**
 * get the number of cells in x and y direction. Meant are only real cells, no ghost cells.
 */
const std::array<int,2> StaggeredGrid::nCells() const {
    return nCells_;
}

const FieldVariable& StaggeredGrid::u() const {
    return u_;
}

const FieldVariable& StaggeredGrid::v() const {
    return v_;
}

const FieldVariable& StaggeredGrid::p() const {
    return p_;
}

double StaggeredGrid::u(int i, int j) const {
    return u_(i, j);
}

double& StaggeredGrid::u(int i, int j) {
    return u_(i, j);
}

double StaggeredGrid::v(int i, int j) const {
    return v_(i, j);
}

double& StaggeredGrid::v(int i, int j) {
    return v_(i, j);
}

double StaggeredGrid::p(int i, int j) const {
    return p_(i, j);
}

double& StaggeredGrid::p(int i, int j) {
    return p_(i, j);
}

double& StaggeredGrid::rhs(int i, int j) {
    return rhs_(i, j);
}

double& StaggeredGrid::f(int i, int j) {
    return f_(i, j);
}

double& StaggeredGrid::g(int i, int j) {
    return g_(i, j);
}

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}

/**
 * 
 */
int StaggeredGrid::uIBegin() const {
    if (containsLeftBoundary_) {
        return 2;
    } else {
        return 1;
    }
}

/**
 * 
 */
int StaggeredGrid::uIEnd() const {
    if (containsRightBoundary_){
        return nCells_[0];
    } else {
        return nCells_[0] + 1;
    }
}

int StaggeredGrid::uJBegin() const {
    return 2;
}

int StaggeredGrid::uJEnd() const {
    return nCells_[1] + 1;
}

int StaggeredGrid::vIBegin() const {
    return 2;
}

int StaggeredGrid::vIEnd() const {
    return nCells_[0] + 1;
}

int StaggeredGrid::vJBegin() const { 
    if (containsBottomBoundary_) {
        return 2;
    } else {
        return 1;
    }
}

int StaggeredGrid::vJEnd() const {
    if (containsTopBoundary_) {
        return nCells_[1];
    } else {
        return nCells_[1] + 1;
    }
}

int StaggeredGrid::pIBegin() const {
    return 2;
}

int StaggeredGrid::pIEnd() const {
    return nCells_[0] + 1;
}

int StaggeredGrid::pJBegin() const {
    return 2;
}

int StaggeredGrid::pJEnd() const {
    return nCells_[1] + 1;
}

CellType StaggeredGrid::cellType(int i, int j) const {
    const int width = nCells_[0] + 3;
    const int height = nCells_[1] + 3;
    if (i < 0 || j < 0 || i >= width || j >= height) {
        return CellType::Obstacle;
    }
    return cellTypes_[static_cast<size_t>(j * width + i)];
}

void StaggeredGrid::setCellType(int i, int j, CellType type) {
    const int width = nCells_[0] + 3;
    const int height = nCells_[1] + 3;
    if (i < 0 || j < 0 || i >= width || j >= height) {
        return;
    }
    cellTypes_[static_cast<size_t>(j * width + i)] = type;
}

bool StaggeredGrid::isFluidCell(int i, int j) const {
    return cellType(i, j) == CellType::Fluid;
}

int StaggeredGrid::countFluidCellsLocal() const {
    int count = 0;
    for (int i = pIBegin(); i <= pIEnd(); ++i) {
        for (int j = pJBegin(); j <= pJEnd(); ++j) {
            if (isFluidCell(i, j)) {
                count++;
            }
        }
    }
    return count;
}

int StaggeredGrid::countActivePressureCellsLocal() const {
    int count = 0;
    for (int i = pIBegin(); i <= pIEnd(); ++i) {
        for (int j = pJBegin(); j <= pJEnd(); ++j) {
            if (isActivePressureCell(i, j)) {
                count++;
            }
        }
    }
    return count;
}

void StaggeredGrid::setPressureReferenceCell(int i, int j) {
    pressureReferenceI_ = i;
    pressureReferenceJ_ = j;
    hasPressureReference_ = (i >= 0 && j >= 0);
}

bool StaggeredGrid::isPressureReferenceCell(int i, int j) const {
    return hasPressureReference_ && i == pressureReferenceI_ && j == pressureReferenceJ_;
}

bool StaggeredGrid::isActivePressureCell(int i, int j) const {
    return isFluidCell(i, j) && !isPressureReferenceCell(i, j);
}

bool StaggeredGrid::hasPressureReferenceCell() const {
    return hasPressureReference_;
}

std::array<int,2> StaggeredGrid::pressureReferenceCell() const {
    return {pressureReferenceI_, pressureReferenceJ_};
}
