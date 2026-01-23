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

    bool hasAnyBoundaryFace() const {
        return topIsBoundaryFace() || rightIsBoundaryFace() || bottomIsBoundaryFace() || leftIsBoundaryFaceC();
    }

    bool topIsBoundaryFace() const {
        return topHasXBC() || topHasYBC();
    }
    bool rightIsBoundaryFace() const {
        return rightHasXBC() || rightHasYBC();
    }
    bool bottomIsBoundaryFace() const {
        return bottomHasXBC() || bottomHasYBC();
    }
    bool leftIsBoundaryFaceC() const {
        return leftHasXBC() || leftHasYBC();
    }

    bool topHasXBC() const {
        return faceTop.dirichletX.has_value() || faceTop.neumannX.has_value();
    }
    bool topHasYBC() const {
        return faceTop.dirichletY.has_value() || faceTop.neumannY.has_value();
    }
    bool rightHasXBC() const {
        return faceRight.dirichletX.has_value() || faceRight.neumannX.has_value();
    }
    bool rightHasYBC() const {
        return faceRight.dirichletY.has_value() || faceRight.neumannY.has_value();
    }
    bool bottomHasXBC() const {
        return faceBottom.dirichletX.has_value() || faceBottom.neumannX.has_value();
    }
    bool bottomHasYBC() const {
        return faceBottom.dirichletY.has_value() || faceBottom.neumannY.has_value();
    }
    bool leftHasXBC() const {
        return faceLeft.dirichletX.has_value() || faceLeft.neumannX.has_value();
    }
    bool leftHasYBC() const {
        return faceLeft.dirichletY.has_value() || faceLeft.neumannY.has_value();
    }
};
