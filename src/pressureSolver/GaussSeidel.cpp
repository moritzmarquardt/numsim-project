#include "pressureSolver/GaussSeidel.hpp"
#include <iostream>


GaussSeidel::GaussSeidel(const std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :

    PressureSolver(discretization, epsilon, maximumNumberOfIterations) {}

void GaussSeidel::solve() {
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
                const double ersterTerm = (discretization_->p(i+1,j) + discretization_->p(i-1,j)) / dx2;
                const double zweiterTerm = (discretization_->p(i,j+1) + discretization_->p(i,j-1)) / dy2;
                discretization_->p(i,j) = lek * (ersterTerm + zweiterTerm - discretization_->rhs(i,j));
            }
        }
        setBoundaryValues();
        computeResidualNorm();
    } ;
    this->numberOfIterations_ = iter;
};