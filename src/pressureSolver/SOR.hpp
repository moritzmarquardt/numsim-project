#pragma once
#include "PressureSolver.hpp"

/**
 * Successive over-relaxation solver. 
 */
class SOR : public PressureSolver {
public :
    SOR(std::shared_ptr<Discretization> discretization, double epsilon, 
        int maximumNumberOfIterations, double omega);
    
        /**
         * solve the system of the Poisson equation for pressure
         * Implements PressureSolver.
         */
        void solve() override;

protected:
    // relaxation factor
    double omega_;
};