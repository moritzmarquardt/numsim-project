#pragma once

#include "pressureSolver/PressureSolver.hpp"

/**
 * Standard Gauss-Seidel method to solve the Poisson equation for pressure.
 */
class GaussSeidel : public PressureSolver {
        
public:
    GaussSeidel(const std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    /**
     * solve the system of the Poisson equation for pressure
     * Implements PressureSolver.
     */
    void solve() override;

};