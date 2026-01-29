#pragma once

#include "domain/boundaryInfo.hpp"
#include <array>

struct CellInfo
{
    BoundaryInfo faceTop;
    BoundaryInfo faceRight;
    BoundaryInfo faceBottom;
    BoundaryInfo faceLeft;
    std::array<int, 2> cellIndexPartition;  //not the global one
    bool fluidCell = true; // true if fluid cell, false if obstacle cell

    bool hasAnyBoundaryFace() const {
        return faceTop.isBoundaryFace() || faceRight.isBoundaryFace() || faceBottom.isBoundaryFace() || faceLeft.isBoundaryFace();
    }

    /**
     * for this cell u and f have to be computed
     * do not calc u if there is a dirichlet condition on the right face
     */
    bool calcUF() const {
        return faceRight.dirichletU.has_value() == false;
    }

    // implement a toString method for pretty printing
    std::string toString() const {
        std::string result = "CellInfo(cellIndexPartition=[" + std::to_string(cellIndexPartition[0]) + ", " + std::to_string(cellIndexPartition[1]) +"], fluidCell=" + (fluidCell ? "true" : "false") + ", ";
        result += "faceTop=" + faceTop.toString() + ", ";
        result += "faceRight=" + faceRight.toString() + ", ";
        result += "faceBottom=" + faceBottom.toString() + ", ";
        result += "faceLeft=" + faceLeft.toString();
        result += ")";
        return result;
    }
};
