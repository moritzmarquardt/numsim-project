#include "computation.hpp";
#include <cmath>
#include <iostream>


void Computation::initialize(int argc, char *argv[]) {
    settings_.loadFromFile(argv[1]);

    // calculate mesh width
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // create discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }

    // create pressure solver
    if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    } else if (settings_.pressureSolver == "SOR") {
        // TODO: calculate optimal omega (not in settings hardcoded)
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else {
        std::cerr << "Error: Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // create output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
}

void Computation::runSimulation() {
    applyInitialBoundaryValues();

    double currentTime = 0.0;
    int iterationCount = 0;

    
}

void Computation::applyInitialBoundaryValues() {
    // Apply initial boundary conditions for u, v, f, g
    for (int j = discretization_->uJBegin()-1; j < discretization_->uJEnd()+1; j++) {
        discretization_->u(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0]; // left
        discretization_->u(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0]; // right

        discretization_->f(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0]; // left
        discretization_->f(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0]; // right
    }
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->v(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1]; // bottom
        discretization_->v(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1]; // top

        discretization_->g(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1]; // bottom
        discretization_->g(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1]; // top
    }
}

void Computation::applyBoundaryValues() {
    // Apply Dirichlet boundary conditions for u and v
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->u(i,discretization_->uJBegin() - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i,discretization_->uJBegin());  // bottom
        discretization_->u(i,discretization_->uJEnd() + 1) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i,discretization_->uJEnd());  // top
    }
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        discretization_->v(discretization_->vIBegin() - 1,j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin(),j); // left
        discretization_->v(discretization_->vIEnd() + 1,j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd(),j); // right
    }
}

void Computation::computePreliminaryVelocities() {
    // calc F and G (leave out boundaries)
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd() - 1; i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            double A_ij = 1 / settings_.re * ( discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))
            - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j) + settings_.g[0];
            discretization_->f(i,j) = discretization_->u(i,j) + A_ij * dt_;
        }
    }

    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd() - 1; j++) {
            double B_ij = 1 / settings_.re * ( discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))
            - discretization_->computeDuvDx(i,j) - discretization_->computeDv2Dy(i,j) + settings_.g[1];
            discretization_->f(i,j) = discretization_->v(i,j) + B_ij * dt_;
        }
    }
}

void Computation::computeRightHandSide() {
    // compute rhs of pressure equation
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            discretization_->rhs(i,j) = ( (discretization_->f(i,j) - discretization_->f(i-1,j)) / discretization_->dx()
                                        + (discretization_->g(i,j) - discretization_->g(i,j-1)) / discretization_->dy() ) / dt_;
        }
    }

}

void Computation::computePressure() {
    // solve pressure equation
    pressureSolver_->solve();
}

void Computation::computeVelocities() {
    // update velocities u and v using F, G, and pressure p
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd() - 1; i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->u(i,j) = discretization_->f(i,j) - discretization_->computeDpDx(i,j) * dt_;
        }
    }

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd() - 1; j++) {
            discretization_->v(i,j) = discretization_->g(i,j) - discretization_->computeDpDy(i,j) * dt_; 
        }
    }
}   

void Computation::computeTimeStepWidth() {
    // compute time step width based on CFL condition

}