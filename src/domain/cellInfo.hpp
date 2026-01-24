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
        return topIsBoundaryFace() || rightIsBoundaryFace() || bottomIsBoundaryFace() || leftIsBoundaryFaceC();
    }

    bool topIsBoundaryFace() const {
        return topHasUBC() || topHasVBC();
    }
    bool rightIsBoundaryFace() const {
        return rightHasUBC() || rightHasVBC();
    }
    bool bottomIsBoundaryFace() const {
        return bottomHasUBC() || bottomHasVBC();
    }
    bool leftIsBoundaryFaceC() const {
        return leftHasUBC() || leftHasVBC();
    }

    bool topHasUBC() const {
        return faceTop.dirichletU.has_value() || faceTop.neumannU.has_value();
    }
    bool topHasVBC() const {
        return faceTop.dirichletV.has_value() || faceTop.neumannV.has_value();
    }
    bool rightHasUBC() const {
        return faceRight.dirichletU.has_value() || faceRight.neumannU.has_value();
    }
    bool rightHasVBC() const {
        return faceRight.dirichletV.has_value() || faceRight.neumannV.has_value();
    }
    bool bottomHasUBC() const {
        return faceBottom.dirichletU.has_value() || faceBottom.neumannU.has_value();
    }
    bool bottomHasVBC() const {
        return faceBottom.dirichletV.has_value() || faceBottom.neumannV.has_value();
    }
    bool leftHasUBC() const {
        return faceLeft.dirichletU.has_value() || faceLeft.neumannU.has_value();
    }
    bool leftHasVBC() const {
        return faceLeft.dirichletV.has_value() || faceLeft.neumannV.has_value();
    }
};
