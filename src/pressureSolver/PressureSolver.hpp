#pragma once

#include "discretization/discretization.hpp"
#include <memory>

/**
 * Interface for the pressure solver. 
 * It computes the pressure field variable such that the continuity equation is fulfilled. 
 */
class PressureSolver {
    public:
        PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

        /**
         * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
         * Implemented in GaussSeidel, and SOR.
         */
        virtual void solve() = 0;

    protected:
        /**
         * set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration 
         */
        void setBoundaryValues();

        virtual void computeResidualNorm();

        double residualNorm();

        //object holding the needed field variables for rhs and p 
        std::shared_ptr<Discretization> discretization_; 
        double epsilon_;
        int maximumNumberOfIterations_;
        double residualNorm_;
        int numberOfIterations_;
};