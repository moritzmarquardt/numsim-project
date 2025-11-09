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
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd() + 1; i++) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd() + 1; j++) {
                double ersterTerm = (discretization_->p(i+1,j) + discretization_->p(i-1,j)) / dx2;
                double zweiterTerm = (discretization_->p(i,j+1) + discretization_->p(i,j-1)) / dy2;
                discretization_->p(i,j) = (1 - omega_)* discretization_->p(i,j) + 
                omega_ * lek * (ersterTerm + zweiterTerm - discretization_->rhs(i,j));
            }
        }
        setBoundaryValues();
        computeResidualNorm();
    };
    this->numberOfIterations_ += iter;
};  