#pragma once
#include "pressureSolver/domainPressureSolver.hpp"

class DomainCG : public DomainPressureSolver {
public:
    DomainCG(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Domain> domain);

    /**
     * solve the system of the Poisson equation for pressure using the Conjugate Gradient method with diagonal preconditioning (multiplying by the inverse of the diagonal elements of the matrix)
     */
    void solve() override;

private:

    /** 
     * Communicate direction_ field variable boundary values with neighboring processes similar to communicateAndSetBoundaryValues for p field variable
     */
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