#pragma once

#include "discretization/discretization.hpp"
#include <memory>

class PressureSolver {
    public:
        PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);
        virtual void solve() = 0;

    protected:
        void setBoundaryValues();

        virtual void computeResidualNorm();

        double residualNorm();

        std::shared_ptr<Discretization> discretization_;
        double epsilon_;
        int maximumNumberOfIterations_;
        double residualNorm_;
        int numberOfIterations_;
};