#pragma once

#include "pressureSolver/parallelPressureSolver.hpp"
#include "partitioning/partitioning.hpp"

class RedBlackGaussSeidel : public ParallelPressureSolver {
public:
    RedBlackGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);
    void solve() override;
};