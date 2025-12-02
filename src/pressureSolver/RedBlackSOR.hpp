#pragma once

#include "parallelPressureSolver.hpp"

class RedBlackSOR : public ParallelPressureSolver {
public:
    RedBlackSOR(std::shared_ptr<Discretization> discretization, double epsilon, 
        int maximumNumberOfIterations, double omega, std::shared_ptr<Partitioning> partitioning);
        
    void solve() override;

protected:
    double omega_;
};