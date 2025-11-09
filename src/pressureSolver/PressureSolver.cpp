#include "pressureSolver/PressureSolver.hpp"
#include "discretization/discretization.hpp"
#include <cmath>
#include <iostream>


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization), 
    epsilon_(epsilon), 
    maximumNumberOfIterations_(maximumNumberOfIterations) {}

    /**
     * setze pressure im äußersten ring auf den wert des jeweils benachbarten inneren cells
     */
void PressureSolver::setBoundaryValues() {
    // links und rechts von unten nach oben
    for (int j = discretization_->pJBegin() - 1; j <= discretization_->pJEnd() + 1; j++) {
        discretization_->p(discretization_->pIBegin() - 1, j) = discretization_->p(discretization_->pIBegin(), j); // links
        discretization_->p(discretization_->pIEnd() + 1, j) = discretization_->p(discretization_->pIEnd(), j); // rechts
    }

    // unten und oben von links nach rechts
    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
        discretization_->p(i, discretization_->pJBegin() - 1) = discretization_->p(i, discretization_->pJBegin()); // unten
        discretization_->p(i, discretization_->pJEnd() + 1) = discretization_->p(i, discretization_->pJEnd()); // oben
    }
}

void PressureSolver::computeResidualNorm() {
    double residual_norm_squared = 0.0;
    const int N = discretization_->nCells()[0] * discretization_->nCells()[1]; // number of real cells
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();

    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
            const double rhs_ij = discretization_->rhs(i,j);
            const double p_xx = (discretization_->p(i+1,j) - 2.0 * discretization_->p(i,j) + discretization_->p(i-1,j)) / dx2;
            const double p_yy = (discretization_->p(i,j+1) - 2.0 * discretization_->p(i,j) + discretization_->p(i,j-1)) / dy2;
            const double residual_ij = rhs_ij - (p_xx + p_yy);
            residual_norm_squared += residual_ij * residual_ij;
        }
    }
    residualNorm_ = residual_norm_squared / N;

}

double PressureSolver::residualNorm() {
    return residualNorm_;
}
    
    