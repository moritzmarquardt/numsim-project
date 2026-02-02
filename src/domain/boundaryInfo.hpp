#pragma once

#include <optional>

struct BoundaryInfo {
    std::optional<double> dirichletU = std::nullopt;
    std::optional<double> dirichletV = std::nullopt;
    std::optional<double> neumannU = std::nullopt; // derivatives of u in normal direction to the boundary face (always outward normal)
    std::optional<double> neumannV = std::nullopt; // derivatives of v in normal direction to the boundary face (always outward normal)

    /**
     * @brief Check if this boundary face has any boundary condition defined
     * @return true if any boundary condition is defined, false otherwise
     */
    bool isBoundaryFace() const {
        return dirichletU.has_value() || dirichletV.has_value() || neumannU.has_value() || neumannV.has_value();
    }

    /**
     * @brief Check if this boundary face has any boundary condition defined for the u component
     * @return true if any boundary condition for u is defined, false otherwise
     */
    bool hasUBC() const {
        return dirichletU.has_value() || neumannU.has_value();
    }

    /**
     * @brief Check if this boundary face has any boundary condition defined for the v component
     * @return true if any boundary condition for v is defined, false otherwise
     */
    bool hasVBC() const {
        return dirichletV.has_value() || neumannV.has_value();
    }


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