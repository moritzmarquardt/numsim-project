#include "SOR.hpp"
#include <iostream>

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega)
        : PressureSolver(discretization, epsilon, maximumNumberOfIterations), 
        omega_(omega) {}

void SOR::solve() {
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();
    const double lek = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps_2 = epsilon_ * epsilon_;

    int iter = 0;

    computeResidualNorm();

    while (iter < maximumNumberOfIterations_ && residualNorm_ > eps_2) {
        iter++;
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
                if (!discretization_->isActivePressureCell(i, j)) {
                    continue;
                }
                const double p_ij = discretization_->p(i,j);
                const double p_e = discretization_->isFluidCell(i+1,j) ? discretization_->p(i+1,j) : p_ij;
                const double p_w = discretization_->isFluidCell(i-1,j) ? discretization_->p(i-1,j) : p_ij;
                const double p_n = discretization_->isFluidCell(i,j+1) ? discretization_->p(i,j+1) : p_ij;
                const double p_s = discretization_->isFluidCell(i,j-1) ? discretization_->p(i,j-1) : p_ij;
                const double ersterTerm = (p_e + p_w) / dx2;
                const double zweiterTerm = (p_n + p_s) / dy2;
                discretization_->p(i,j) = (1 - omega_)* p_ij + 
                omega_ * lek * (ersterTerm + zweiterTerm - discretization_->rhs(i,j));
            }
        }
        setBoundaryValues();
        computeResidualNorm();
    };
    this->numberOfIterations_ += iter;
};  
