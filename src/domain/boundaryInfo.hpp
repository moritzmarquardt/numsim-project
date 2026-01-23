#pragma once

#include <optional>

struct BoundaryInfo {
    std::optional<double> dirichletU = std::nullopt;
    std::optional<double> dirichletV = std::nullopt;
    std::optional<double> neumannU = std::nullopt; // derivatives of u in normal direction to the boundary face (always outward normal)
    std::optional<double> neumannV = std::nullopt; // derivatives of v in normal direction to the boundary face
    // only neumar or dirichlet conditions can be set for u and v on each face
};