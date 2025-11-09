#include "computation.hpp"
#include <cmath>
#include <iostream>


void Computation::initialize(int argc, char *argv[]) {
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();

    // calculate mesh width
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // create discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }

    std::cout << "Created discretization with mesh width dx: " << meshWidth_[0] << ", dy: " << meshWidth_[1] << std::endl;
    std::cout << "Number of cells in x direction: " << settings_.nCells[0] << ", y direction: " << settings_.nCells[1] << std::endl;

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

    while (currentTime < settings_.endTime) {
        applyBoundaryValues();
        std::cout << "Applied boundary values." << std::endl;
        // print current u and v for debugging
        // std::cout << "Current u field:" << std::endl;
        // discretization_->u().printAsArray();
        // std::cout << "Current v field:" << std::endl;
        // discretization_->v().printAsArray();

        computeTimeStepWidth();

        if (currentTime + dt_ > settings_.endTime) {
            dt_ = settings_.endTime - currentTime;
        }
        std::cout << "Computed time step width: " << dt_ << std::endl;


        computePreliminaryVelocities();
        std::cout << "Computed preliminary velocities." << std::endl;
        computeRightHandSide();
        std::cout << "Computed right hand side." << std::endl;
        computePressure();
        std::cout << "Computed pressure." << std::endl;
        computeVelocities();
        std::cout << "Computed velocities." << std::endl;

        currentTime += dt_;
        iterationCount++;
        std::cout << "Advanced to time: " << currentTime << std::endl;
        std::cout << "Completed iteration: " << iterationCount << std::endl;

        outputWriterParaview_->writeFile(currentTime);
        std::cout << "Wrote Paraview output." << std::endl;
        outputWriterText_->writeFile(currentTime);
        std::cout << "Wrote text output." << std::endl;

        std::cout << "Iteration: " << iterationCount << ", Time: " << currentTime << ", dt: " << dt_ << std::endl;
        
    }

    
}

void Computation::applyInitialBoundaryValues() {
    // Apply initial boundary conditions for u, v, f, g
    // u and f lay at the right side of a cell, so the left boundary is at index uIBegin() - 1 and the right boundary at uIEnd() + 1
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++) {
        // go through all j indices, that have real cells for u.
        const int i_left_bc = discretization_->uIBegin() - 1;
        const int i_right_bc = discretization_->uIEnd() + 1;
        discretization_->u(i_left_bc, j) = settings_.dirichletBcLeft[0]; // left
        discretization_->u(i_right_bc, j) = settings_.dirichletBcRight[0]; // right

        discretization_->f(i_left_bc, j) = settings_.dirichletBcLeft[0]; // left
        discretization_->f(i_right_bc, j) = settings_.dirichletBcRight[0]; // right
    }
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
        // go through all i indices, that have real cells for v.
        const int j_bottom_bc = discretization_->vJBegin() - 1;
        const int j_top_bc = discretization_->vJEnd() + 1;
        discretization_->v(i, j_bottom_bc) = settings_.dirichletBcBottom[1]; // bottom
        discretization_->v(i, j_top_bc) = settings_.dirichletBcTop[1]; // top

        discretization_->g(i, j_bottom_bc) = settings_.dirichletBcBottom[1]; // bottom
        discretization_->g(i, j_top_bc) = settings_.dirichletBcTop[1]; // top
    }
}

void Computation::applyBoundaryValues() {
    // Apply Dirichlet boundary conditions for u and v
    std::cout << "Applying boundary values for u and v." << std::endl;
    std::cout << "apply top bc from index i " << discretization_->uIBegin() << " to " << discretization_->uIEnd() << std::endl;
    for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++) { 
        // go through all i indices.
        const int j_bottom_bc = discretization_->uJBegin() - 1;
        const int j_bottom_inside = discretization_->uJBegin();
        const int j_top_bc = discretization_->uJEnd() + 1;
        const int j_top_inside = discretization_->uJEnd();
        discretization_->u(i,j_bottom_bc) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i,j_bottom_inside);  // bottom
        discretization_->u(i,j_top_bc) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i,j_top_inside);  // top
    }
    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++) {
        // go through all j indices.
        const int i_left_bc = discretization_->vIBegin() - 1;
        const int i_left_inside = discretization_->vIBegin();
        const int i_right_bc = discretization_->vIEnd() + 1;
        const int i_right_inside = discretization_->vIEnd();
        discretization_->v(i_left_bc,j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(i_left_inside,j); // left
        discretization_->v(i_right_bc,j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(i_right_inside,j); // right
    }
}

void Computation::computePreliminaryVelocities() {
    // calc F and G (leave out boundaries)
    for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++) {
            double A_ij = 1 / settings_.re * ( discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))
            - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j) + settings_.g[0];
            discretization_->f(i,j) = discretization_->u(i,j) + A_ij * dt_;
        }
    }

    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++) {
            double B_ij = 1 / settings_.re * ( discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))
            - discretization_->computeDuvDx(i,j) - discretization_->computeDv2Dy(i,j) + settings_.g[1];
            discretization_->g(i,j) = discretization_->v(i,j) + B_ij * dt_;
        }
    }
}

void Computation::computeRightHandSide() {
    // compute rhs of pressure equation
    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
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
    for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++) {
            discretization_->u(i,j) = discretization_->f(i,j) - discretization_->computeDpDx(i,j) * dt_;
        }
    }

    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++) {
            discretization_->v(i,j) = discretization_->g(i,j) - discretization_->computeDpDy(i,j) * dt_; 
        }
    }
}   

void Computation::computeTimeStepWidth() {
    // compute time step width based on CFL condition
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    double dx2 = dx * dx;
    double dy2 = dy * dy;  

    double dt_diff_cond = 0.5 * settings_.re * (dx2 * dy2) / (dx2 + dy2);
    double dt_conv_cond = std::min(dx / discretization_->u().computeMaxAbs(), dy / discretization_->v().computeMaxAbs());

    double dt_prelim = settings_.tau * std::min(dt_diff_cond, dt_conv_cond);

    if (dt_prelim < settings_.maximumDt) {
        dt_ = dt_prelim;
    } else {
        dt_ = settings_.maximumDt;
        std::cout << "Warning: Time step width limited by maximumDt!" << std::endl;
    }
}