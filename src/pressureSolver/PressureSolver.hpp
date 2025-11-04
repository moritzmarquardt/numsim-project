#pragma once

#include "discretization.hpp"

class PressureSolver {
    public:
        PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);
        virtual void solve() = 0;

    protected:
        void setBoundaryValues();

        std::shared_ptr<Discretization> discretization_;
        double epsilon_;
        int maximumNumberOfIterations_;
};