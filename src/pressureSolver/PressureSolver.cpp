#include "PressureSolver.hpp"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations) {}

void PressureSolver::setBoundaryValues() {
}
    
    