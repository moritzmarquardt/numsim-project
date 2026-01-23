#pragma once

#include <optional>

struct BoundaryInfo {
    std::optional<double> dirichletX = std::nullopt;
    std::optional<double> dirichletY = std::nullopt;
    std::optional<double> neumannX = std::nullopt;
    std::optional<double> neumannY = std::nullopt;
};