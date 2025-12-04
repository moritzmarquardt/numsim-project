#pragma once

#include "pressureSolver/parallelPressureSolver.hpp"
#include "partitioning/partitioning.hpp"

class RedBlackGaussSeidel : public ParallelPressureSolver {
public:
    RedBlackGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);
    
    /**
     * solve the system of the Poisson equation for pressure using Red-Black Gauss-Seidel method
     * communicates values for ghost nodes after red and black updates
     */
    void solve() override;
};