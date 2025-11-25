#pragma once

#include "discretization/discretization.hpp"
#include "pressureSolver/PressureSolver.hpp"
#include "partitioning/partitioning.hpp"
#include <mpi.h>

class ParallelPressureSolver : public PressureSolver {
public:
    ParallelPressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);
    virtual void solve() = 0;


protected:
    void computeResidualNorm() override;
    void communicateAndSetBoundaryValues();

    std::shared_ptr<Partitioning> partitioning_;
};