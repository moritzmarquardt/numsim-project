#include "pressureSolver/parallelPressureSolver.hpp"
#include "pressureSolver/PressureSolver.hpp"
#include "discretization/discretization.hpp"
#include "partitioning/partitioning.hpp"
#include <cmath>

ParallelPressureSolver::ParallelPressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning) {}

void ParallelPressureSolver::computeResidualNorm() {
    //TODO: implement parallel residual norm computation
}

void ParallelPressureSolver::communicateAndSetBoundaryValues() {
    //TODO: implement communication of boundary values between neighboring partitions
}