#pragma once
#include "PressureSolver.hpp"

class SOR : public PressureSolver {
public :
    SOR(std::shared_ptr<Discretization> discretization, double epsilon, 
        int maximumNumberOfIterations, double omega);
    
        void solve() override;

protected:
    double omega_;
};