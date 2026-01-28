#pragma once

#include <optional>

struct BoundaryInfo {
    std::optional<double> dirichletU = std::nullopt;
    std::optional<double> dirichletV = std::nullopt;
    std::optional<double> neumannU = std::nullopt; // derivatives of u in normal direction to the boundary face (always outward normal)
    std::optional<double> neumannV = std::nullopt; // derivatives of v in normal direction to the boundary face
    // only neumar or dirichlet conditions can be set for u and v on each face
    bool isPartitionInnerFace = false; // true if this face is a partition boundary face inside the domain. false if it is a physical boundary face or inside the partition

    // implement a toString method for pretty printing
    std::string toString() const {
        std::string result = "BoundaryInfo(";
        if (dirichletU.has_value()) {
            result += "dirichletU=" + std::to_string(dirichletU.value()) + ", ";
        }
        if (dirichletV.has_value()) {
            result += "dirichletV=" + std::to_string(dirichletV.value()) + ", ";
        }
        if (neumannU.has_value()) {
            result += "neumannU=" + std::to_string(neumannU.value()) + ", ";
        }
        if (neumannV.has_value()) {
            result += "neumannV=" + std::to_string(neumannV.value()) + ", ";
        }
        if (result.size() > 13) { // length of "BoundaryInfo("
            result.pop_back(); // remove last comma and space
            result.pop_back();
        }
        result += ")";
        return result;
    }
};