#pragma once
#include "pressureSolver/parallelPressureSolver.hpp"

class ParallelCG : public ParallelPressureSolver {
public:
    ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);

    void solve() override;

private:
    void communicateAndSetBoundaryValuesForDirection();

    FieldVariable direction_;
    FieldVariable residual_;
    FieldVariable w_;

    double residualOld2_;
    double alpha_;
    double residualNew2_;
    double dx2_;
    double dy2_;

};