#pragma once

#include "discretization/discretization.hpp"
#include "pressureSolver/PressureSolver.hpp"
#include "partitioning/partitioning.hpp"
#include "domain/domain.hpp"
#include <mpi.h>

class DomainPressureSolver : public PressureSolver {
public:
    DomainPressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Domain> domain);
    
    /**
     * solve the system of the Poisson equation for pressure
     * has to be implemented in derived classes
     */
    virtual void solve() = 0;


protected:
    void computeResidualNorm() override;
    void communicateAndSetBoundaryValues();

    std::shared_ptr<Partitioning> partitioning_;
    std::shared_ptr<Domain> domain_;
    MPI_Comm cartComm_;
};