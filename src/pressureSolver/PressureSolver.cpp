#include "pressureSolver/PressureSolver.hpp"
#include "discretization/discretization.hpp"
#include <cmath>
#include <iostream>


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization), 
    epsilon_(epsilon), 
    maximumNumberOfIterations_(maximumNumberOfIterations) {}


void PressureSolver::setBoundaryValues() {
    //setze pressure im äußersten ring auf den wert des jeweils benachbarten inneren cells
    // links und rechts von unten nach oben
    // set the corner values in the j-iteration
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
    int N = 0;
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();

    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
            if (!discretization_->isActivePressureCell(i, j)) {
                continue;
            }
            N++;
            const double rhs_ij = discretization_->rhs(i,j);
            const double p_ij = discretization_->p(i,j);
            const double p_e = discretization_->isFluidCell(i+1,j) ? discretization_->p(i+1,j) : p_ij;
            const double p_w = discretization_->isFluidCell(i-1,j) ? discretization_->p(i-1,j) : p_ij;
            const double p_n = discretization_->isFluidCell(i,j+1) ? discretization_->p(i,j+1) : p_ij;
            const double p_s = discretization_->isFluidCell(i,j-1) ? discretization_->p(i,j-1) : p_ij;
            const double p_xx = (p_e - 2.0 * p_ij + p_w) / dx2;
            const double p_yy = (p_n - 2.0 * p_ij + p_s) / dy2;
            const double residual_ij = rhs_ij - (p_xx + p_yy);
            residual_norm_squared += residual_ij * residual_ij;
        }
    }
    if (N > 0) {
        residualNorm_ = residual_norm_squared / N;
    } else {
        residualNorm_ = 0.0;
    }

}

double PressureSolver::residualNorm() {
    return residualNorm_;
}
    
    
