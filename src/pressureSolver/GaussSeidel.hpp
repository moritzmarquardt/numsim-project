#pragma once

#include "pressureSolver/PressureSolver.hpp"

class GaussSeidel : public PressureSolver {
        
public:
    GaussSeidel::GaussSeidel(const std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);
    void GaussSeidel::solve() override;

};