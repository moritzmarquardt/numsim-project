#pragma once

#include "parallelPressureSolver.hpp"

class RedBlackSOR : public ParallelPressureSolver {
public:
    RedBlackSOR(std::shared_ptr<Discretization> discretization, double epsilon, 
        int maximumNumberOfIterations, double omega, std::shared_ptr<Partitioning> partitioning);
    
    /**
     * solve the system of the Poisson equation for pressure using Red-Black SOR method
     * communicates values for ghost nodes after red and black updates
     */
    void solve() override;

protected:
    double omega_;
};