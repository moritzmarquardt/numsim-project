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




