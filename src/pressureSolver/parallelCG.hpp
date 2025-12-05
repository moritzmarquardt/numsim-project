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

    double rTr_;
    double alpha_;
    double beta_;
    double rTrNew_;
    double dx2_;
    double dy2_;

};